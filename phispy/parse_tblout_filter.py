#!/usr/bin/env python
import sys
import csv
from pathlib import Path


# ---------------- 参数检查 ----------------
if len(sys.argv) != 3:
    sys.exit("Usage: python parse_tblout_filter.py <tblout> <summary.tsv>")
tblout, summary = map(Path, sys.argv[1:3])

# ---------------- 数据库加载 ----------------
db = Path("/projects/alvarez/virus_hmm_db")
INT = set(open(db / "integrase_unified.list").read().split())
STR = set(open(db / "structural_unified.list").read().split())

# {pid: {'int': bool, 'str': bool}}
hits = {}

# ---------------- 解析 tblout ----------------
with tblout.open() as f:
    for line in f:
        if line.startswith('#'):
            continue
        cols = line.split()
        qname, profile = cols[0], cols[2]

        # ---------- 关键修改 ----------
        parts = qname.split('_')
        pid = '_'.join(parts[:-1])          # 同一 prophage 共用前缀
        # --------------------------------

        rec = hits.setdefault(pid, {'int': False, 'str': False})
        if profile in INT:
            rec['int'] = True
        if profile in STR:
            rec['str'] = True

# ---------------- 写入汇总 ----------------
with summary.open('w', newline='') as out:
    w = csv.writer(out, delimiter='\t')
    w.writerow(['prophage', 'status', 'note'])
    for pid, flag in hits.items():
        keep = flag['int'] and flag['str']
        note = []
        if not flag['int']:
            note.append('no_integrase')
        if not flag['str']:
            note.append('no_structural')
        w.writerow([pid, 'kept' if keep else 'removed', ','.join(note)])

print(f"[DONE] summary written to {summary}")
