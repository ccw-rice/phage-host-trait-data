#!/bin/bash
#SBATCH -J grodon2_batch_v2
#SBATCH -p commons
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 03:00:00
#SBATCH -o /scratch/%u/logs/grodon2.v2.%A_%a.out
#SBATCH -e /scratch/%u/logs/grodon2.v2.%A_%a.err

set -Eeuo pipefail

# Usage:
#   sbatch --array=1-N /projects/alvarez/scripts/grodon2_batch.v2.sbatch MANIFEST [OUTROOT]
#
#   MANIFEST: one /abs/path/to/<SAMPLE>/cds_from_genomic.fna per line
#   OUTROOT : default = /scratch/$USER/grodon2_out

MANIFEST="${1:?Need manifest file}"
OUTROOT="${2:-/scratch/$USER/grodon2_out}"

mkdir -p "$OUTROOT" "/scratch/$USER/logs"

# Explicitly specify Rscript and library path
# (ensures the environment with installed packages is used)
RSCRIPT="/scratch/$USER/envs/grodon2_env/bin/Rscript"
export R_LIBS_USER="/scratch/$USER/envs/grodon2_env/lib/R/library"

if [[ ! -x "$RSCRIPT" ]]; then
  echo "[ERROR] $RSCRIPT not found or not executable" >&2
  exit 127
fi

# Retrieve input for this array task (line N of the manifest)
IN_FFN="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST" || true)"

[[ -z "${IN_FFN:-}" ]] && {
  echo "[WARN] Empty MANIFEST line for task ${SLURM_ARRAY_TASK_ID}"
  exit 0
}

[[ ! -f "$IN_FFN" ]] && {
  echo "[ERROR] Input file not found: $IN_FFN"
  exit 2
}

# Sample name = parent directory name
SAMPLE="$(basename "$(dirname "$IN_FFN")")"
IN_DIR="$(dirname "$IN_FFN")"
OUT_DIR="$OUTROOT/$SAMPLE"
OUT_CSV="$OUT_DIR/${SAMPLE}_gRodon.csv"

mkdir -p "$OUT_DIR"

echo "[INFO] Host=$(hostname)"
echo "[INFO] SAMPLE=$SAMPLE"
echo "[INFO] ffn=$IN_FFN"
echo "[INFO] out=$OUT_CSV"

# === Critical: export variables for use in R ===
export SAMPLE
export IN_FFN
export OUT_CSV

# Skip if output already exists and is non-empty
# (supports checkpoint/restart workflows)
if [[ -s "$OUT_CSV" ]]; then

  echo "[SKIP] $OUT_CSV"

else

  "$RSCRIPT" - <<'EOF'
suppressPackageStartupMessages({
  library(Biostrings)
  library(gRodon)  # Package name is gRodon (installed version: 2.5.1)
})

sample  <- Sys.getenv("SAMPLE")
ffn     <- Sys.getenv("IN_FFN")
out_csv <- Sys.getenv("OUT_CSV")

if (!nzchar(sample) || !file.exists(ffn) || !nzchar(out_csv)) {
  stop("Missing SAMPLE/IN_FFN/OUT_CSV or file not found.")
}

# Read CDS sequences (DNAStringSet)
genes <- Biostrings::readDNAStringSet(ffn)

if (length(genes) == 0L)
  stop("No CDS found in: ", ffn)

# Identify highly expressed genes (robust matching)
he <- grepl(
  "ribosomal protein|\\brps\\b|\\brpl\\b|elongation factor",
  names(genes),
  ignore.case = TRUE
)

# If no HE genes are detected, issue a warning but continue
# (all FALSE values are still acceptable)
if (!any(he)) {
  warning("No HE genes matched in headers; proceeding with he=FALSE for all CDS.")
  he <- rep(FALSE, length(genes))
}

# New interface: predictGrowth(genes, highly_expressed)
res <- gRodon::predictGrowth(genes, he)

# Convert results to a one-row data.frame and add custom fields
df <- as.data.frame(t(unlist(res)), stringsAsFactors = FALSE)

df$organism <- sample

# If doubling time ("d", hours) is returned, calculate mu_max
if ("d" %in% names(df)) {
  df$mu_max <- log(2) / as.numeric(df$d)
}

# Write output
utils::write.csv(df, out_csv, row.names = FALSE)
EOF

fi

# === Thread-safe summary aggregation ===
SUMMARY="$OUTROOT/summary_grodon2.csv"
LOCK="$SUMMARY.lock"

{
  flock 200

  if [[ ! -f "$SUMMARY" ]]; then
    head -n 1 "$OUT_CSV" > "$SUMMARY"
  fi

  tail -n +2 "$OUT_CSV" >> "$SUMMARY"

} 200>"$LOCK"

echo "[DONE] $SAMPLE -> $OUT_CSV"
