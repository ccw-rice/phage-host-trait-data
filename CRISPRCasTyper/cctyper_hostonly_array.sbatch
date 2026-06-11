#!/bin/bash
#SBATCH --job-name=cctyper_host
#SBATCH --partition=commons
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH -o /scratch/cw106/logs/cctyper_host_%A_%a.out
#SBATCH -e /scratch/cw106/logs/cctyper_host_%A_%a.err

set -euo pipefail

echo "[INFO] Hostname: $(hostname)"
echo "[INFO] SLURM_JOB_ID=${SLURM_JOB_ID}, TASK_ID=${SLURM_ARRAY_TASK_ID}"

##################################
# 1) Activate CCTyper environment
##################################
module purge
module use /opt/apps/modules/Core
module --ignore_cache load Miniforge3/24.1.2-0

eval "$(conda shell.bash hook)"
conda activate /scratch/cw106/cctyper_env_new

echo "[INFO] Conda env: $CONDA_DEFAULT_ENV"
echo "[INFO] which cctyper:"
which cctyper || { echo "[ERROR] cctyper not found in PATH"; exit 1; }

##################################
# 2) Set CCTYPER_DB
##################################
export CCTYPER_DB=/scratch/cw106/cctyper_db/data
if [[ ! -d "$CCTYPER_DB" ]]; then
    echo "[ERROR] CCTYPER_DB not found at $CCTYPER_DB"
    exit 1
fi
echo "[INFO] Using CCTYPER_DB=$CCTYPER_DB"

##################################
# 3) Manifest and output directory (host-only version)
##################################
MANIFEST=/scratch/cw106/displace_prophage/crispr_defense_manifest_hostonly.txt
OUTDIR_ROOT=/scratch/cw106/displace_prophage/crisprcas_out_hostonly
mkdir -p "$OUTDIR_ROOT"

##################################
# 4) Retrieve the input .fna file for this task
##################################
fa=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST" || true)

if [[ -z "${fa}" ]]; then
    echo "[WARN] No input file found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}; exiting"
    exit 0
fi

if [[ ! -s "$fa" ]]; then
    echo "[ERROR] Input file not found or empty: $fa"
    exit 1
fi

base=$(basename "$fa" .fna)   # e.g. GCF_000005845.2_host_only
outdir="${OUTDIR_ROOT}/${base}_cctyper"

##################################
# 5) Handle existing output directory
##################################
if [[ -d "$outdir" ]]; then
    # If the directory already contains files, assume the job completed successfully
    if [[ $(ls -A "$outdir" 2>/dev/null | wc -l) -gt 0 ]]; then
        echo "[SKIP] $base already has CCTyper output in $outdir"
        exit 0
    else
        echo "[WARN] $outdir exists but is empty; removing and rerunning"
        rm -rf "$outdir"
    fi
fi

echo "[INFO] [$SLURM_ARRAY_TASK_ID] Running CCTyper on: $fa"
echo "[INFO] Output dir: $outdir"

##################################
# 6) Run CCTyper
##################################
set -x
cctyper -t 4 "$fa" "$outdir"
set +x

echo "[DONE] [$SLURM_ARRAY_TASK_ID] Completed $base"
