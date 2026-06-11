#!/bin/bash
#SBATCH -J barrnap_batch
#SBATCH -p commons
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH -o /scratch/%u/logs/barrnap.%A_%a.out
#SBATCH -e /scratch/%u/logs/barrnap.%A_%a.err

set -Eeuo pipefail

# Usage:
#   sbatch --array=1-N /projects/alvarez/scripts/barrnap_batch.v3.sbatch MANIFEST [OUTROOT] [KINGDOM]
#
# Arguments:
#   MANIFEST : File containing one absolute path to a .fna/.fa/.fasta file per line
#   OUTROOT  : Output root directory (default: /scratch/$USER/barrnap_out)
#   KINGDOM  : bac/arc/euk (default: bac)

MANIFEST="${1:?Need MANIFEST file}"
OUTROOT="${2:-/scratch/$USER/barrnap_out}"
KINGDOM="${3:-bac}"

# Environment path
ENV="/scratch/$USER/envs/barrnap_env"

# Absolute path to barrnap (avoid dependency on PATH lookup)
BARRNAP="$ENV/bin/barrnap"
if [[ ! -x "$BARRNAP" ]]; then
  echo "[ERROR] $BARRNAP not found or not executable." >&2
  exit 127
fi

mkdir -p "/scratch/$USER/logs" "$OUTROOT"

# Retrieve the input FASTA file for this array task (remove potential CR characters)
INFA="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST" | tr -d '\r' || true)"

if [[ -z "${INFA:-}" ]]; then
  echo "[WARN] Empty line for task ${SLURM_ARRAY_TASK_ID}"
  exit 0
fi

if [[ ! -f "$INFA" ]]; then
  echo "[ERROR] Input file not found: $INFA"
  exit 2
fi

# Sample name: parent directory name (e.g. 2617270923);
# if unavailable, use the input filename without extension
SAMPLE="$(basename "$(dirname "$INFA")")"

if [[ -z "$SAMPLE" || "$SAMPLE" == "/" ]]; then
  SAMPLE="$(basename "$INFA")"
  SAMPLE="${SAMPLE%%.*}"
fi

OUTDIR="$OUTROOT/$SAMPLE"
mkdir -p "$OUTDIR"

# Output GFF file name: retain only the basename without extension
BASENAME="$(basename "$INFA")"
OUTGFF="$OUTDIR/${BASENAME%%.*}_rRNA.gff"

# Ensure hmmer and bedtools are taken from this environment
export PATH="$ENV/bin:$PATH"

# Skip existing non-empty output files (supports restart/resume)
if [[ -s "$OUTGFF" ]]; then
  echo "[SKIP] $OUTGFF"
else
  echo "[INFO] Host: $(hostname)"
  echo "[INFO] Input: $INFA"
  echo "[INFO] Out  : $OUTGFF"
  echo "[INFO] KINGDOM=$KINGDOM  Threads=${SLURM_CPUS_PER_TASK}"

  # Use the environment-specific Perl interpreter to avoid system Perl
  "$ENV/bin/perl" "$BARRNAP" \
    --kingdom "$KINGDOM" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    "$INFA" > "$OUTGFF"
fi

# === Count 16S / 23S / 5S rRNA copy numbers and write summary CSV ===

SUMMARY="$OUTROOT/summary_16S_copy_number.csv"
LOCK="$SUMMARY.lock"

# Create summary table header if needed (with file lock)
{
  flock 200

  if [[ ! -f "$SUMMARY" ]]; then
    echo "sample,input,copy_16S,copy_23S,copy_5S" > "$SUMMARY"
  fi

} 200>"$LOCK"

# Parse GFF and count rRNA genes
count_16S=$(awk '$3=="rRNA" && $9 ~ /16S_rRNA/' "$OUTGFF" | wc -l | awk '{print $1}')
count_23S=$(awk '$3=="rRNA" && $9 ~ /23S_rRNA/' "$OUTGFF" | wc -l | awk '{print $1}')
count_5S=$(awk  '$3=="rRNA" && $9 ~ /5S_rRNA/'  "$OUTGFF" | wc -l | awk '{print $1}')

# Append results to summary CSV (with file lock to avoid race conditions)
{
  flock 200

  printf "%s,%s,%s,%s,%s\n" \
    "$SAMPLE" \
    "$INFA" \
    "$count_16S" \
    "$count_23S" \
    "$count_5S" >> "$SUMMARY"

} 200>"$LOCK"

echo "[DONE] $SAMPLE | 16S=$count_16S 23S=$count_23S 5S=$count_5S"
