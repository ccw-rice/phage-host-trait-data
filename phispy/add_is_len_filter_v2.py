#!/usr/bin/env python3
"""
add_is_len_filter_v2.py
=======================

Batch example
-------------
python add_is_len_filter_v2.py \
  --root         /scratch/cw106/phispy_output \
  --isescan-root /scratch/cw106/isescan_output
"""

import argparse
import csv
import glob
import os
from collections import defaultdict
from pathlib import Path


# ---------- Basic utilities ----------
def merge_intervals(intervals):
    if not intervals:
        return []

    intervals = sorted(intervals, key=lambda x: x[0])

    out = []
    s, e = intervals[0]

    for a, b in intervals[1:]:
        if a <= e + 1:  # overlapping / adjacent
            e = max(e, b)
        else:
            out.append((s, e))
            s, e = a, b

    out.append((s, e))
    return out


def read_is(csv_path):

    with open(csv_path, "r") as fh:
        first = fh.readline()

    delim = "," if first.count(",") > first.count("\t") else "\t"

    hits = defaultdict(list)

    with open(csv_path, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter=delim)

        for row in rdr:
            try:
                a, b = int(row["isBegin"]), int(row["isEnd"])
            except (KeyError, ValueError):
                continue

            contig = row["seqID"].strip("|").split("|")[-1]
            hits[contig].append((a, b))

    return hits


def calc_is_len(hit_list, start, end):
    """
    Calculate the total overlap length between IS elements
    and the prophage region [start, end].

    Any overlap is counted, but only the intersecting
    coordinates contribute to the final length.
    """

    segs = []

    for a, b in hit_list:

        if b < start or a > end:  # completely left/right
            continue

        segs.append(
            (max(a, start), min(b, end))
        )  # overlap coordinates

    if not segs:
        return 0

    merged = merge_intervals(segs)

    return sum(e - s + 1 for s, e in merged)


# ---------- Process one genome ----------
def run_one(filter_tsv, csv_path, zero_is=False):

    raw_db = {} if zero_is else read_is(csv_path)

    out_tsv = os.path.splitext(filter_tsv)[0] + "_islen.tsv"

    with open(filter_tsv, newline="") as fin, \
         open(out_tsv, "w", newline="") as fout:

        rd = csv.DictReader(fin, delimiter="\t")

        hdr = rd.fieldnames + [
            "contig",
            "length",
            "IS_bp",
            "IS_pct"
        ]

        wr = csv.DictWriter(fout, hdr, delimiter="\t")
        wr.writeheader()

        for row in rd:

            contig, st, ed = row["prophage"].rsplit("_", 2)

            s, e = int(st), int(ed)
            length = e - s + 1

            bp = calc_is_len(
                raw_db.get(contig, []),
                s,
                e
            )

            pct = bp / length if length else 0.0

            note_list = [
                n.strip()
                for n in row.get("note", "").split(";")
                if n.strip()
            ]

            status = row.get("status", "kept")

            if length < 18_000:
                status = "discard"
                note_list.append("len<18kb")

            if pct > 0.25:
                status = "discard"
                note_list.append("IS>25%")

            row.update(
                status=status,
                note="; ".join(sorted(set(note_list))),
                contig=contig,
                length=length,
                IS_bp=bp,
                IS_pct=f"{pct:.3f}",
            )

            wr.writerow(row)

    print(f"✔ Output written: {out_tsv}")


# ---------- Arguments ----------
def parse_args():

    ap = argparse.ArgumentParser()

    mx = ap.add_mutually_exclusive_group(required=True)

    mx.add_argument(
        "--filter",
        help="Single-sample filter_summary.tsv"
    )

    mx.add_argument(
        "--root",
        help="Batch mode: PhiSpy output root directory"
    )

    ap.add_argument(
        "--csv",
        help="Single-sample ISEScan CSV"
    )

    ap.add_argument(
        "--gff",
        help="Single-sample ISEScan GFF (CSV inferred automatically)"
    )

    ap.add_argument(
        "--isescan-root",
        help="Batch mode: ISEScan output root directory"
    )

    return ap.parse_args()


# ---------- Main ----------
def main():

    args = parse_args()

    # -------- Batch mode --------
    if args.root:

        filt_files = glob.glob(
            os.path.join(args.root, "*", "filter_summary.tsv")
        )

        if not filt_files:
            raise FileNotFoundError(
                f"No filter_summary.tsv found under {args.root}"
            )

        for f in sorted(filt_files):

            gcf_full = os.path.basename(
                os.path.dirname(f)
            )  # GCF_xxx_genomic

            gcf_base = gcf_full.removesuffix(
                "_genomic"
            )  # GCF_xxx

            rootdir = Path(
                args.isescan_root or args.root
            )

            cand_dirs = [
                rootdir / gcf_full / "phispy_input",
                rootdir / f"{gcf_base}_genomic" / "phispy_input",
                rootdir / gcf_base / "phispy_input",
                rootdir / gcf_full,
                rootdir / gcf_base,
                rootdir / gcf_base / gcf_base,
                rootdir / gcf_full / gcf_base,
            ]

            csv_p = None

            for d in cand_dirs:

                if not d.is_dir():
                    continue

                exact = d / f"{gcf_base}_genomic.fna.csv"

                if exact.is_file():
                    csv_p = str(exact)
                    break

                hits = list(d.glob("*.fna.csv"))

                if hits:
                    csv_p = str(hits[0])
                    break

            if csv_p:
                run_one(f, csv_p)

            else:
                print(
                    f"⚠ Missing ISEScan CSV: {gcf_full} → treated as 0 IS"
                )
                run_one(
                    f,
                    csv_path=None,
                    zero_is=True
                )

    # -------- Single-sample mode --------
    else:

        if args.csv:

            run_one(
                args.filter,
                args.csv
            )

        elif args.gff:

            csv_p = (
                args.gff
                .replace(".gff.gz", ".csv")
                .replace(".gff", ".csv")
            )

            run_one(
                args.filter,
                csv_p
            )

        else:
            raise ValueError(
                "Single-sample mode requires --csv or --gff"
            )


if __name__ == "__main__":
    main()
