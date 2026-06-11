#!/usr/bin/env python3

import os
import sys
import csv
import glob
import gzip
from pathlib import Path

LIST_TSV   = os.environ.get("LIST_TSV",   "/scratch/cw106/phage_id_list_2col.ALL15k.tsv")
OUT_ROOT   = os.environ.get("OUT_ROOT",   "/scratch/cw106/phage_hmm_out_all15k")
DB_DIR     = os.environ.get("DB_DIR",     "/projects/alvarez/virus_hmm_db")
OUT_TSV    = os.environ.get("OUT_TSV",    "/scratch/cw106/phage_lifestyle_calls_all15k.tsv")

# Root directory containing original phage genome FASTA files
# (used for calculating genome length)
DATA_ROOT  = os.environ.get("DATA_ROOT",  "/scratch/cw106/phagescope_data_resource")

# Gene lists (prefer unified lists when available)
INTEG_LIST = os.environ.get("INTEG_LIST", f"{DB_DIR}/integrase_unified.list")
TERL_LIST  = os.environ.get("TERL_LIST",  f"{DB_DIR}/terL.list")
LYSIS_LIST = os.environ.get("LYSIS_LIST", f"{DB_DIR}/lysis_genes.list")

FA_EXTS = (".fasta", ".fa", ".fna")


def load_set(fp):
    s = set()
    with open(fp) as f:
        for ln in f:
            ln = ln.strip()
            if ln and not ln.startswith("#"):
                s.add(ln.split()[0])
    return s


def read_hmm_ids(tbl_path):
    """
    Read hmmsearch --tblout output.

    Verified format:
    Column 3 (query name) = HMM ID
    (e.g., phrog_xxx, VOG_xxx, etc.)
    """
    ids = set()

    with open(tbl_path) as f:
        for ln in f:
            if not ln or ln[0] == "#":
                continue

            parts = ln.rstrip("\n").split()

            if len(parts) >= 3:
                ids.add(parts[2])

    return ids


def lifestyle_call(has_int, has_terl, has_lysis):

    # Current relaxed classification:
    # integrase present -> temperate
    # no integrase -> initially considered virulent

    if has_int:
        return "temperate"

    # Structural evidence for dsDNA tailed phages
    # plus lysis support -> strong virulent

    if has_terl and has_lysis:
        return "strong_virulent"

    if has_terl:
        return "virulent"

    # No terL:
    # could be non-dsDNA-tailed,
    # incomplete contig,
    # or missing annotation

    return "uncertain"


def open_maybe_gz(path):

    if path.endswith(".gz"):
        return gzip.open(path, "rt")

    return open(path, "r")


def find_fasta(data_root, src, pid):
    """
    Expected layout:

        {DATA_ROOT}/{src}/{src}/{pid}.(fasta|fa|fna)[.gz]
    """

    base = Path(data_root) / src / src

    # Exact match first
    for ext in FA_EXTS:

        p = base / f"{pid}{ext}"

        if p.exists() and p.stat().st_size > 0:
            return str(p)

        pgz = base / f"{pid}{ext}.gz"

        if pgz.exists() and pgz.stat().st_size > 0:
            return str(pgz)

    # Fallback:
    # search pid.* for fasta/fa/fna or compressed equivalents

    for f in glob.glob(str(base / f"{pid}.*")):

        if f.endswith(FA_EXTS) or any(
            f.endswith(e + ".gz") for e in FA_EXTS
        ):
            try:
                if os.path.getsize(f) > 0:
                    return f
            except OSError:
                pass

    return None


def fasta_length_bp(fasta_path):
    """
    Sum sequence lengths across all contigs.

    Header lines are ignored.

    For multi-contig genomes,
    returns the total bp across all contigs.
    """

    total = 0

    with open_maybe_gz(fasta_path) as f:

        for ln in f:

            if not ln:
                continue

            if ln.startswith(">"):
                continue

            total += len(ln.strip())

    return total


def main():

    # Input validation

    for p in [LIST_TSV, INTEG_LIST, TERL_LIST, LYSIS_LIST]:

        if not os.path.isfile(p) or os.path.getsize(p) == 0:
            print(f"[ERR] missing/empty: {p}", file=sys.stderr)
            sys.exit(2)

    integ = load_set(INTEG_LIST)
    terl  = load_set(TERL_LIST)
    lysis = load_set(LYSIS_LIST)

    # Read master table

    master = []

    with open(LIST_TSV) as f:

        r = csv.reader(f, delimiter="\t")

        next(r, None)

        for row in r:

            if not row or len(row) < 2:
                continue

            src = row[0].strip()
            pid = row[1].strip()

            master.append((src, pid))

    out_dir = Path(OUT_ROOT)

    rows = []

    miss_tbl = 0
    miss_fa  = 0

    for src, pid in master:

        # Genome length

        fasta = find_fasta(DATA_ROOT, src, pid)

        if fasta is None:

            genome_len = ""
            miss_fa += 1

        else:

            try:
                genome_len = fasta_length_bp(fasta)

            except Exception:
                genome_len = ""
                miss_fa += 1

        # HMM results

        tbl = out_dir / src / f"{pid}.hmm.tbl"

        if not tbl.exists() or tbl.stat().st_size == 0:

            miss_tbl += 1

            rows.append([
                src,
                pid,
                genome_len,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                "missing_tbl"
            ])

            continue

        hmm_ids = read_hmm_ids(str(tbl))

        n_total = len(hmm_ids)

        hit_int  = hmm_ids & integ
        hit_terl = hmm_ids & terl
        hit_lys  = hmm_ids & lysis

        has_int  = 1 if hit_int else 0
        has_terl = 1 if hit_terl else 0
        has_lys  = 1 if hit_lys else 0

        call = lifestyle_call(
            has_int,
            has_terl,
            has_lys
        )

        rows.append([
            src,
            pid,
            genome_len,
            n_total,
            len(hit_int),
            len(hit_terl),
            len(hit_lys),
            has_int,
            has_terl,
            has_lys,
            call
        ])

    Path(OUT_TSV).parent.mkdir(
        parents=True,
        exist_ok=True
    )

    with open(OUT_TSV, "w", newline="") as fo:

        w = csv.writer(
            fo,
            delimiter="\t"
        )

        w.writerow([
            "Phage_source",
            "Phage_ID",
            "Genome_length_bp",
            "n_unique_hmm_hits",
            "n_integrase_hits",
            "n_terL_hits",
            "n_lysis_hits",
            "Has_integrase",
            "Has_terL",
            "Has_lysis",
            "Lifestyle_call"
        ])

        w.writerows(rows)

    print("[DONE]", OUT_TSV)

    print(
        "[INFO] master:",
        len(master),
        "missing_tbl:",
        miss_tbl,
        "missing_fasta:",
        miss_fa
    )


if __name__ == "__main__":
    main()
