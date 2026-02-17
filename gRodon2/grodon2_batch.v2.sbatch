#!/bin/bash
#SBATCH -J grodon2_batch_v2
#SBATCH -p commons
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 03:00:00
#SBATCH -o /scratch/%u/logs/grodon2.v2.%A_%a.out
#SBATCH -e /scratch/%u/logs/grodon2.v2.%A_%a.err

set -Eeuo pipefail

# 用法：
#   sbatch --array=1-N /projects/alvarez/scripts/grodon2_batch.v2.sbatch MANIFEST [OUTROOT]
#     MANIFEST：每行一个 /abs/path/to/<SAMPLE>/cds_from_genomic.fna
#     OUTROOT ：默认 /scratch/$USER/grodon2_out

MANIFEST="${1:?Need manifest file}"
OUTROOT="${2:-/scratch/$USER/grodon2_out}"

mkdir -p "$OUTROOT" "/scratch/$USER/logs"

# 显式指定 Rscript 与库路径（保证用到你安装好包的环境）
RSCRIPT="/scratch/$USER/envs/grodon2_env/bin/Rscript"
export R_LIBS_USER="/scratch/$USER/envs/grodon2_env/lib/R/library"

if [[ ! -x "$RSCRIPT" ]]; then
  echo "[ERROR] $RSCRIPT not found or not executable" >&2
  exit 127
fi

# 取该数组任务的输入（清单第 N 行）
IN_FFN="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST" || true)"
[[ -z "${IN_FFN:-}" ]] && { echo "[WARN] Empty MANIFEST line for task ${SLURM_ARRAY_TASK_ID}"; exit 0; }
[[ ! -f "$IN_FFN" ]] && { echo "[ERROR] Input file not found: $IN_FFN"; exit 2; }

# 样本名 = 上一级目录名
SAMPLE="$(basename "$(dirname "$IN_FFN")")"
IN_DIR="$(dirname "$IN_FFN")"
OUT_DIR="$OUTROOT/$SAMPLE"
OUT_CSV="$OUT_DIR/${SAMPLE}_gRodon.csv"
mkdir -p "$OUT_DIR"

echo "[INFO] Host=$(hostname)"
echo "[INFO] SAMPLE=$SAMPLE"
echo "[INFO] ffn=$IN_FFN"
echo "[INFO] out=$OUT_CSV"

# === 关键：导出变量给 R 使用 ===
export SAMPLE
export IN_FFN
export OUT_CSV

# 已有且非空则跳过（支持断点续跑）
if [[ -s "$OUT_CSV" ]]; then
  echo "[SKIP] $OUT_CSV"
else
  "$RSCRIPT" - <<'EOF'
suppressPackageStartupMessages({
  library(Biostrings)
  library(gRodon)  # 包名 gRodon（你安装的是 2.5.1）
})

sample  <- Sys.getenv("SAMPLE")
ffn     <- Sys.getenv("IN_FFN")
out_csv <- Sys.getenv("OUT_CSV")

if (!nzchar(sample) || !file.exists(ffn) || !nzchar(out_csv)) {
  stop("Missing SAMPLE/IN_FFN/OUT_CSV or file not found.")
}

# 读取 CDS（DNAStringSet）
genes <- Biostrings::readDNAStringSet(ffn)
if (length(genes) == 0L) stop("No CDS found in: ", ffn)

# 识别高表达基因（尽量稳健）
he <- grepl("ribosomal protein|\\brps\\b|\\brpl\\b|elongation factor",
            names(genes), ignore.case = TRUE)

# 若一个都没匹配到，给出警告但允许继续（全 FALSE 也可运行）
if (!any(he)) {
  warning("No HE genes matched in headers; proceeding with he=FALSE for all CDS.")
  he <- rep(FALSE, length(genes))
}

# 新接口：predictGrowth(genes, highly_expressed)
res <- gRodon::predictGrowth(genes, he)

# 整理为一行 data.frame，并追加自定义字段
df <- as.data.frame(t(unlist(res)), stringsAsFactors = FALSE)
df$organism <- sample
# 若返回名为 "d" 的倍增时间(h)，计算 mu_max
if ("d" %in% names(df)) {
  df$mu_max <- log(2) / as.numeric(df$d)
}

# 输出
utils::write.csv(df, out_csv, row.names = FALSE)
EOF
fi

# === 并发安全地汇总 ===
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
