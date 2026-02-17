#!/usr/bin/env python3
"""
add_is_len_filter_v2.py
=======================

批量示例
--------
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


# ---------- 基础工具 ----------
def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: x[0])

    out = []
    s, e = intervals[0]
    for a, b in intervals[1:]:
        if a <= e + 1:  # 重叠 / 相邻
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
    计算 IS 与 prophage 区段 [start, end] 的交集总长度
    （任何相交都会计入，只取交集部分）
    """
    segs = []
    for a, b in hit_list:
        if b < start or a > end:  # 完全在左 / 右侧
            continue
        segs.append((max(a, start), min(b, end)))  # 交集坐标

    if not segs:
        return 0
    merged = merge_intervals(segs)
    return sum(e - s + 1 for s, e in merged)


# ---------- 处理单 genome ----------
def run_one(filter_tsv, csv_path, zero_is=False):
    raw_db = {} if zero_is else read_is(csv_path)

    out_tsv = os.path.splitext(filter_tsv)[0] + "_islen.tsv"
    with open(filter_tsv, newline="") as fin, open(out_tsv, "w", newline="") as fout:
        rd = csv.DictReader(fin, delimiter="\t")
        hdr = rd.fieldnames + ["contig", "length", "IS_bp", "IS_pct"]
        wr = csv.DictWriter(fout, hdr, delimiter="\t")
        wr.writeheader()

        for row in rd:
            contig, st, ed = row["prophage"].rsplit("_", 2)
            s, e = int(st), int(ed)
            length = e - s + 1

            bp = calc_is_len(raw_db.get(contig, []), s, e)
            pct = bp / length if length else 0.0

            note_list = [n.strip() for n in row.get("note", "").split(";") if n.strip()]
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

    print(f"✔ 已输出: {out_tsv}")


# ---------- 参数 ----------
def parse_args():
    ap = argparse.ArgumentParser()
    mx = ap.add_mutually_exclusive_group(required=True)
    mx.add_argument("--filter", help="单样本 filter_summary.tsv")
    mx.add_argument("--root", help="批量：PhiSpy 输出根目录")
    ap.add_argument("--csv", help="单样本 ISEScan CSV")
    ap.add_argument("--gff", help="单样本 ISEScan GFF（推断 csv）")
    ap.add_argument("--isescan-root", help="批量：ISEScan 输出根目录")
    return ap.parse_args()


# ---------- 主入口 ----------
def main():
    args = parse_args()

    # -------- 批量模式 --------
    if args.root:
        filt_files = glob.glob(os.path.join(args.root, "*", "filter_summary.tsv"))
        if not filt_files:
            raise FileNotFoundError(f"{args.root} 下未找到 filter_summary.tsv")

        for f in sorted(filt_files):
            gcf_full = os.path.basename(os.path.dirname(f))       # GCF_..._genomic
            gcf_base = gcf_full.removesuffix("_genomic")          # GCF_...
            rootdir = Path(args.isescan_root or args.root)

            cand_dirs = [
                rootdir / gcf_full / "phispy_input",
                rootdir / f"{gcf_base}_genomic" / "phispy_input",
                rootdir / gcf_base / "phispy_input",
                rootdir / gcf_full,                # GCF_..._genomic
                rootdir / gcf_base,                # GCF_...
                rootdir / gcf_base / gcf_base,     # GCF_xxx/GCF_xxx
                rootdir / gcf_full / gcf_base,     # GCF_xxx_genomic/GCF_xxx
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
                print(f"⚠ 跳过 ISEScan CSV: {gcf_full} → 视为 0 IS")
                run_one(f, csv_path=None, zero_is=True)

    # -------- 单样本模式 --------
    else:
        if args.csv:
            run_one(args.filter, args.csv)
        elif args.gff:
            csv_p = args.gff.replace(".gff.gz", ".csv").replace(".gff", ".csv")
            run_one(args.filter, csv_p)
        else:
            raise ValueError("单样本需 --csv 或 --gff")


if __name__ == "__main__":
    main()
