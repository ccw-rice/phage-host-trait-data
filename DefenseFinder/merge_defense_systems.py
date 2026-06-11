#!/usr/bin/env python3
"""
merge_defense_systems.py   (for NOTSX)

Merge and summarize:
/scratch/cw106/df_out/**/**/*_defense_finder_systems.tsv

Outputs:
  all_genomes_defense_systems.tsv
  systems_count_by_type.tsv
"""

from pathlib import Path
import pandas as pd
import sys, os

ROOT = Path("/scratch/cw106/df_out2")          # Root directory

KEEP = {
    "GCF_ID", "accession", "replicon",
    "sys_id", "type", "subtype", "activity",
    "sys_beg", "sys_end", "proteins_in_system",
    "genes_count", "list_of_profiles",
}

print(">>> Root :", ROOT)

tsv_files = sorted(ROOT.glob("**/*_defense_finder_systems.tsv"))
print(">>> Found:", len(tsv_files), "systems.tsv files")

if not tsv_files:
    sys.exit("[ERROR] No systems.tsv files found")

frames = []

for f in tsv_files:
    df = pd.read_csv(f, sep="\t")

    # Standardize genome identifier column:
    # accession / replicon → GCF_ID
    if "GCF_ID" not in df.columns:

        if "accession" in df.columns:
            df = df.rename(columns={"accession": "GCF_ID"})

        elif "replicon" in df.columns:
            df = df.rename(columns={"replicon": "GCF_ID"})

        else:
            # Fallback: extract GCF accession from filename
            gcf = f.name.split("_")[0]
            df["GCF_ID"] = gcf

    df = df[[c for c in df.columns if c in KEEP]]

    df["source_file"] = f.name

    frames.append(df)

all_df = pd.concat(frames, ignore_index=True)

out1 = ROOT / "all_genomes_defense_systems.tsv"
all_df.to_csv(out1, sep="\t", index=False)

print(">>> Merged table written to:", out1)

pivot = (
    all_df.groupby(["GCF_ID", "type"], as_index=False)
          .size()
          .pivot(index="GCF_ID", columns="type", values="size")
          .fillna(0)
          .astype(int)
)

out2 = ROOT / "systems_count_by_type.tsv"
pivot.to_csv(out2, sep="\t")

print(">>> Summary table written to:", out2)
print(">>> Done! ")
