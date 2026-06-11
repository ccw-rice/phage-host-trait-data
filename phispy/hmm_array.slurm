#!/bin/bash
#SBATCH --job-name=hmmsearch_phage
#SBATCH --partition=commons
#SBATCH --array=5753
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=02:00:00
#SBATCH -o /projects/alvarez/logs/hmm_%A_%a.out
#SBATCH -e /projects/alvarez/logs/hmm_%A_%a.err

# ────────── 1. Environment ──────────
module purge
module use /opt/apps/modules/Core
module --ignore_cache load Miniforge3/24.11.3-0

eval "$(conda shell.bash hook)"
conda activate phispy_env

set -euo pipefail

# ────────── 2. Genome ID ──────────
GENOME=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" /projects/alvarez/genomes.list)

[[ -d /scratch/cw106/phispy_output/${GENOME}_genomic ]] && \
    GENOME="${GENOME}_genomic"

# ────────── 3. Source and destination paths ──────────
SRC_DIR=/scratch/cw106/phispy_output/$GENOME

FAA_SRC=$SRC_DIR/phage.faa
TBL_DST=$SRC_DIR/hmm.tbl

[[ ! -s "$FAA_SRC" ]] && {
    echo "[SKIP] $FAA_SRC is missing or empty"
    exit 0
}

# ────────── 4. Temporary working directory ──────────
SCR=/scratch/$USER/$SLURM_JOB_ID/$GENOME

mkdir -p "$SCR"

cp "$FAA_SRC" "$SCR/phage.faa"

# ────────── 5. Run HMMER ──────────
echo "[RUN ] hmmsearch on $GENOME"

hmmsearch \
    --cpu 4 \
    -E 1e-3 \
    --tblout "$SCR/hmm.tbl" \
    /projects/alvarez/virus_hmm_db/merged_phage.hmm \
    "$SCR/phage.faa" > /dev/null

echo "[DONE] Generated hmm.tbl"

# ────────── 6. Sync results and clean up ──────────
rsync -av "$SCR/hmm.tbl" "$TBL_DST"

rm -rf "$SCR"
