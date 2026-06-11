#!/bin/bash
#SBATCH -J prophage_filter
#SBATCH --array=5753
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 01:00:00
#SBATCH -o /projects/alvarez/logs/filter_%A_%a.out
#SBATCH -e /projects/alvarez/logs/filter_%A_%a.err

# ────────── 1. Environment ──────────
module purge
module use /opt/apps/modules/Core          # Core directory contains the Miniforge3 modulefile
module --ignore_cache load Miniforge3/24.11.3-0

eval "$(conda shell.bash hook)"
conda activate phispy_env

set -euo pipefail

# ────────── 2. Genome ID ──────────
GENOME=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" /projects/alvarez/genomes.list)

# ────────── 3. Source directory under /scratch ──────────
dir="/scratch/cw106/phispy_output/${GENOME}"

[[ -d "${dir}_genomic" ]] && dir="${dir}_genomic"

TBL_SRC="${dir}/hmm.tbl"
SUM_DST="${dir}/filter_summary.tsv"

[[ ! -s "$TBL_SRC" ]] && {
    echo "[SKIP] Missing $TBL_SRC"
    exit 0
}

# ────────── 4. Working directory on /scratch ──────────
SCR=/scratch/$USER/$SLURM_JOB_ID/$GENOME

mkdir -p "$SCR"

cp "$TBL_SRC" "$SCR/"

# ────────── 5. Run filtering script (on /scratch) ──────────
python /projects/alvarez/scripts/parse_tblout_filter.py \
       "$SCR/hmm.tbl" \
       "$SCR/filter_summary.tsv"

# ────────── 6. Sync results back to /scratch and clean up ──────────
rsync -av "$SCR/filter_summary.tsv" "$SUM_DST"

rm -rf "$SCR"
