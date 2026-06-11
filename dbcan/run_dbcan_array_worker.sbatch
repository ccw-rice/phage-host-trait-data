#!/bin/bash -l
#SBATCH -J dbcan_array
#SBATCH -p commons
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH -o /projects/alvarez/scripts/logs/dbcan.%A_%a.out
#SBATCH -e /projects/alvarez/scripts/logs/dbcan.%A_%a.err

set -Eeuo pipefail
set -x
trap 'echo "[ERR] line:$LINENO cmd:$BASH_COMMAND" >&2' ERR

# --- Configuration (modify as needed) ---
export MAMBA_ROOT_PREFIX=/scratch/cw106/software
ENV=/scratch/cw106/software/envs/dbcan_env
DBROOT=/scratch/cw106/dbcan_db
INDIR=/scratch/cw106/eggnog_input/ncbi_dataset/data
OUTBASE=/scratch/cw106/run_dbcan_out
LIST=/projects/alvarez/scripts/gcf_list.txt
OUT_PRE=run1
CPU=${SLURM_CPUS_PER_TASK:-16}
# ----------------------------------------

# Get the GCF corresponding to the current array task
GCF=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$LIST" || true)
if [ -z "${GCF:-}" ]; then
  echo "[WARN] No GCF for task $SLURM_ARRAY_TASK_ID"
  exit 0
fi

INFA="$INDIR/$GCF/protein.faa"
OUTD="$OUTBASE/$GCF"
mkdir -p "$OUTD"

# Activate environment
module purge || true
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook --shell bash)"
micromamba activate "$ENV"
export PATH="$ENV/bin:$PATH"
export PYTHONNOUSERSITE=1

# Sanity checks
if [ ! -s "$INFA" ]; then
  echo "[WARN] Missing input: $INFA ; skipping"
  exit 0
fi

test -s "$DBROOT/dbCAN-HMMdb-V14.txt" || {
  echo "[FATAL] Missing: $DBROOT/dbCAN-HMMdb-V14.txt"
  exit 2
}

test -s "$DBROOT/CAZy.dmnd" || {
  echo "[FATAL] Missing: $DBROOT/CAZy.dmnd"
  exit 2
}

which python
python --version

echo "[INFO] GCF=$GCF"

diamond dbinfo -d "$DBROOT/CAZy.dmnd" | head -n 5 || true

# Run dbCAN (run-dbcan==2.0.11 has no global --cpu option;
# CPU threads must be specified separately for HMMER and DIAMOND)
python "$ENV/bin/run_dbcan.py" \
  "$INFA" protein \
  --out_dir "$OUTD" \
  --db_dir "$DBROOT" \
  --dbCANFile dbCAN-HMMdb-V14.txt \
  --out_pre "$OUT_PRE" \
  --tools hmmer diamond \
  --hmm_cpu "$CPU" --hmm_eval 1e-15 --hmm_cov 0.35 \
  --dia_cpu "$CPU" --dia_eval 1e-102

# Standardize output naming and perform validation
OV1="$OUTD/${OUT_PRE}overview.txt"

[ -s "$OV1" ] && ln -sf "$OV1" "$OUTD/overview.tsv"

if [ ! -s "$OUTD/overview.tsv" ] && [ ! -s "$OV1" ]; then
  if [ ! -s "$OUTD/${OUT_PRE}hmmer.out" ] && \
     [ ! -s "$OUTD/${OUT_PRE}diamond.out" ]; then

    echo "[FATAL] No outputs found in $OUTD"
    ls -lah "$OUTD" || true
    exit 4
  fi
fi

echo "[OK] Done $GCF"
