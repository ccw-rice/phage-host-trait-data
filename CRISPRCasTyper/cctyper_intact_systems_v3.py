#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cctyper_intact_systems_v3.1.py  (CI/CA-free, optional --merge-near-as-orphan)

Counting rules:
1) An "intact system" is defined solely by presence in cas_operons.tab (certain operons); CI/CA checks are no longer performed.
2) System type statistics: Best_type (or Prediction if Best_type is missing) are mapped to Type I–VI, and subtype information is retained.
3) Cas gene count: calculated as the summed length of the Genes column in cas_operons.tab.
4) Spacer count: only Trusted arrays associated with intact operons are counted.
   Spacer counts are obtained preferentially from spacers/CRISPR_ID.fa record counts;
   if unavailable, N_repeats - 1 is used as a fallback.
   Each array is counted only once.
5) Orphan arrays: preferentially obtained from crisprs_orphan.tab (Trusted only).
   If unavailable, Trusted arrays from crisprs_all.tab that are absent from
   crisprs_near_cas.tab are treated as orphan arrays. Spacer counting follows
   the same procedure described above.
6) --merge-near-as-orphan: optionally merge Trusted near-Cas arrays that were
   not assigned to any intact system into the orphan count.

Output columns:
- sample, has_intact_system, intact_systems, intact_cas_genes, intact_spacers
- Type_I..Type_VI, intact_subtypes, intact_operons
- orphan_trusted_arrays_count, orphan_trusted_spacers
"""
from __future__ import annotations
import argparse, ast, re, sys, glob
from pathlib import Path
import pandas as pd

ROMAN_RE = re.compile(r"(?:Type[_-])?([IVX]+)", re.I)

# ───────── IO ─────────
def load_tab(p: Path) -> pd.DataFrame:
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(p, sep="\t", comment="#", low_memory=False)

def ensure_col(df: pd.DataFrame, want: str, fallback_index: int | None = 1) -> pd.DataFrame:
    """In some versions the target column is stored as the second column; rename it as a fallback."""
    if df.empty or want in df.columns:
        return df
    if fallback_index is not None and df.shape[1] > fallback_index:
        df = df.copy()
        df.rename(columns={df.columns[fallback_index]: want}, inplace=True)
    return df

# ───────── parsing ─────────
def parse_list(s):
    if pd.isna(s): return []
    if isinstance(s, list): return s
    t = str(s).strip()
    t = (t.replace("’","'").replace("‘","'")
           .replace("“",'"').replace("”",'"'))
    try:
        v = ast.literal_eval(t)
        return list(v) if isinstance(v,(list,tuple)) else [v]
    except Exception:
        t = t.strip("[]")
        if not t: return []
        return [x.strip().strip("'").strip('"') for x in re.split(r",\s*", t)]

def top_level_type(name: str) -> str|None:
    if not isinstance(name,str): return None
    m = ROMAN_RE.search(name)
    return f"Type_{m.group(1).upper()}" if m else None

def bool_from_str(x) -> bool:
    if isinstance(x,bool): return x
    return str(x).strip().lower() in ("true","t","1","yes")

# ───────── spacers ─────────
def n_spacers_from_repeats(n_repeats) -> int:
    try:
        n = int(float(str(n_repeats).strip()))
        return max(0, n-1)
    except Exception:
        return 0

def count_fasta_records(fa: Path) -> int:
    try:
        c = 0
        with fa.open("r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith(">"):
                    c += 1
        return c
    except Exception:
        return 0

def spacers_count_for_array(array_id: str, spacers_dir: Path, fallback_repeats=None) -> int:
    fa = spacers_dir / f"{array_id}.fa"
    if fa.exists():
        n = count_fasta_records(fa)
        if n > 0:
            return n
    return n_spacers_from_repeats(fallback_repeats)

# ───────── core ─────────
def summarise_dir(d: Path, merge_near_as_orphan: bool) -> dict:
    sample = d.name.replace("_cctyper","")

    oper = load_tab(d/"cas_operons.tab")
    cc   = load_tab(d/"CRISPR_Cas.tab")
    near = load_tab(d/"crisprs_near_cas.tab")
    allc = load_tab(d/"crisprs_all.tab")
    orph = load_tab(d/"crisprs_orphan.tab")  # Empty dataframe if file is absent or empty

    res = {
        "sample": sample,
        "has_intact_system": "NO",
        "intact_systems": 0,
        "intact_cas_genes": 0,
        "intact_spacers": 0,
        "Type_I": 0, "Type_II": 0, "Type_III": 0, "Type_IV": 0, "Type_V": 0, "Type_VI": 0,
        "intact_subtypes": "",
        "intact_operons": "",
        "orphan_trusted_arrays_count": 0,
        "orphan_trusted_spacers": 0,
    }

    if oper.empty:
        return res

    # Type / subtype / operon statistics and Cas gene count
    intact_operons = []
    intact_subtypes = []
    type_counts = {"Type_I":0,"Type_II":0,"Type_III":0,"Type_IV":0,"Type_V":0,"Type_VI":0}
    genes_total = 0

    oper = ensure_col(oper, "Operon", 0)
    for _, r in oper.iterrows():
        op = str(r["Operon"]).strip()
        intact_operons.append(op)

        subtype = None
        if "Best_type" in oper.columns and not pd.isna(r.get("Best_type")):
            s = parse_list(r["Best_type"])
            subtype = s[0] if s else None
        if not subtype:
            subtype = r.get("Prediction", None)

        if isinstance(subtype, str):
            intact_subtypes.append(subtype)
            tl = top_level_type(subtype)
            if tl in type_counts:
                type_counts[tl] += 1

        genes_total += len(parse_list(r.get("Genes", [])))

    res["has_intact_system"] = "YES"
    res["intact_systems"] = len(intact_operons)
    res["intact_cas_genes"] = int(genes_total)
    res["intact_operons"] = ",".join(intact_operons)
    res["intact_subtypes"] = ",".join(intact_subtypes)

    for k, v in type_counts.items():
        res[k] = v

    # Count spacers from Trusted arrays associated with intact operons
    spacers_dir = d / "spacers"
    counted_arrays = set()
    spacers_sum = 0

    if not cc.empty:
        cc = ensure_col(cc, "Operon", 0)

        if "CRISPRs" not in cc.columns:
            cc["CRISPRs"] = ""

        near = ensure_col(near, "CRISPR", 1) if not near.empty else pd.DataFrame(columns=["CRISPR","Trusted","N_repeats"])
        near_min = near[["CRISPR","Trusted","N_repeats"]].copy() if not near.empty else pd.DataFrame(columns=["CRISPR","Trusted","N_repeats"])

        if not near_min.empty:
            near_min["Trusted_bool"] = near_min["Trusted"].apply(bool_from_str)
            near_map = near_min.set_index("CRISPR")
        else:
            near_map = pd.DataFrame()

        valid_ops = set(intact_operons)

        for _, row in cc.iterrows():
            if str(row["Operon"]).strip() not in valid_ops:
                continue

            crisprs = parse_list(row.get("CRISPRs",""))

            for c in crisprs:
                if not near_map.empty and c in near_map.index and bool(near_map.loc[c, "Trusted_bool"]):
                    if c not in counted_arrays:
                        spacers_sum += spacers_count_for_array(
                            c,
                            spacers_dir,
                            fallback_repeats=near_map.loc[c, "N_repeats"]
                        )
                        counted_arrays.add(c)

    res["intact_spacers"] = int(spacers_sum)

    # Orphan arrays
    orphan_ids = set()
    orphan_sp = 0

    if not orph.empty:
        orph = ensure_col(orph, "CRISPR", 1)
        orph_min = orph[["CRISPR","Trusted","N_repeats"]].copy()
        orph_min["Trusted_bool"] = orph_min["Trusted"].apply(bool_from_str)
        orph_min = orph_min[orph_min["Trusted_bool"] == True].drop_duplicates(subset=["CRISPR"])

        for _, r in orph_min.iterrows():
            cid = str(r["CRISPR"])

            if cid in counted_arrays:
                continue

            orphan_ids.add(cid)
            orphan_sp += spacers_count_for_array(
                cid,
                spacers_dir,
                fallback_repeats=r["N_repeats"]
            )

    else:
        if not allc.empty:
            allc = ensure_col(allc, "CRISPR", 1)

            all_min = allc[["CRISPR","Trusted","N_repeats"]].copy()
            all_min["Trusted_bool"] = all_min["Trusted"].apply(bool_from_str)

            near = ensure_col(near, "CRISPR", 1) if not near.empty else pd.DataFrame(columns=["CRISPR"])
            near_set = set(near["CRISPR"]) if not near.empty else set()

            orphan_df = all_min[
                (all_min["Trusted_bool"] == True) &
                (~all_min["CRISPR"].isin(near_set))
            ].drop_duplicates(subset=["CRISPR"])

            for _, r in orphan_df.iterrows():
                cid = str(r["CRISPR"])

                if cid in counted_arrays:
                    continue

                orphan_ids.add(cid)
                orphan_sp += spacers_count_for_array(
                    cid,
                    spacers_dir,
                    fallback_repeats=r["N_repeats"]
                )

    # Optionally merge Trusted near-Cas arrays not assigned to any system
    if merge_near_as_orphan:

        if 'near_min' not in locals():
            near = ensure_col(near, "CRISPR", 1) if not near.empty else pd.DataFrame(columns=["CRISPR","Trusted","N_repeats"])
            near_min = near[["CRISPR","Trusted","N_repeats"]].copy() if not near.empty else pd.DataFrame(columns=["CRISPR","Trusted","N_repeats"])

            if not near_min.empty:
                near_min["Trusted_bool"] = near_min["Trusted"].apply(bool_from_str)

        if not near_min.empty:
            merge_ids = set()
            merge_sp = 0

            for _, r in near_min.iterrows():
                cid = str(r["CRISPR"])

                if not bool(r["Trusted_bool"]):
                    continue

                if cid in counted_arrays or cid in orphan_ids:
                    continue

                merge_ids.add(cid)
                merge_sp += spacers_count_for_array(
                    cid,
                    spacers_dir,
                    fallback_repeats=r["N_repeats"]
                )

            orphan_ids |= merge_ids
            orphan_sp += merge_sp

    res["orphan_trusted_arrays_count"] = int(len(orphan_ids))
    res["orphan_trusted_spacers"] = int(orphan_sp)

    return res

# ───────── CLI ─────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root", help="Root directory containing multiple *_cctyper result folders")
    ap.add_argument("--pattern", default="*_cctyper",
                    help="Glob pattern used to identify result directories under root")
    ap.add_argument("--merge-near-as-orphan", action="store_true",
                    help="Merge Trusted near-Cas arrays that are not assigned to any system into orphan counts")
    args = ap.parse_args()

    root = Path(args.root).expanduser().resolve()
    dirs = [Path(p).resolve() for p in glob.glob(str(root / args.pattern))]

    if not dirs:
        sys.exit(f"[ERR] no dirs matched {args.pattern} under {root}")

    rows = [
        summarise_dir(
            d,
            args["merge_near_as_orphan"] if isinstance(args, dict)
            else args.merge_near_as_orphan
        )
        for d in sorted(dirs)
    ]

    pd.DataFrame(rows).to_csv(sys.stdout, sep="\t", index=False)

if __name__ == "__main__":
    main()
