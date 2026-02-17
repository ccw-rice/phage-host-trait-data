#!/bin/bash
#SBATCH -J barrnap_batch
#SBATCH -p commons
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH -o /scratch/%u/logs/barrnap.%A_%a.out
#SBATCH -e /scratch/%u/logs/barrnap.%A_%a.err

set -Eeuo pipefail

# 用法：
#   sbatch --array=1-N /projects/alvarez/scripts/barrnap_batch.v3.sbatch MANIFEST [OUTROOT] [KINGDOM]
# 参数说明：
#   MANIFEST : 每行一个绝对路径的 .fna/.fa/.fasta
#   OUTROOT  : 输出根目录，默认 /scratch/$USER/barrnap_out
#   KINGDOM  : bac/arc/euk（默认 bac）

MANIFEST="${1:?Need MANIFEST file}"
OUTROOT="${2:-/scratch/$USER/barrnap_out}"
KINGDOM="${3:-bac}"

# 你的环境路径
ENV="/scratch/$USER/envs/barrnap_env"

# 绝对路径 barrnap（避免依赖 PATH 查找）
BARRNAP="$ENV/bin/barrnap"
if [[ ! -x "$BARRNAP" ]]; then
  echo "[ERROR] $BARRNAP not found or not executable." >&2
  exit 127
fi

mkdir -p "/scratch/$USER/logs" "$OUTROOT"

# 取该数组任务的输入 fasta（去掉潜在的 \r）
INFA="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST" | tr -d '\r' || true)"
if [[ -z "${INFA:-}" ]]; then
  echo "[WARN] Empty line for task ${SLURM_ARRAY_TASK_ID}"
  exit 0
fi
if [[ ! -f "$INFA" ]]; then
  echo "[ERROR] Input file not found: $INFA"
  exit 2
fi

# 样本名：上一级目录名（例如 2617270923），若为空则用去后缀的文件名
SAMPLE="$(basename "$(dirname "$INFA")")"
if [[ -z "$SAMPLE" || "$SAMPLE" == "/" ]]; then
  SAMPLE="$(basename "$INFA")"; SAMPLE="${SAMPLE%%.*}"
fi

OUTDIR="$OUTROOT/$SAMPLE"
mkdir -p "$OUTDIR"

# 输出 GFF 文件名：去掉扩展只保留基名
BASENAME="$(basename "$INFA")"
OUTGFF="$OUTDIR/${BASENAME%%.*}_rRNA.gff"

# 保证 hmmer/bedtools 来自该环境
export PATH="$ENV/bin:$PATH"

# 已有且非空则跳过（方便断点续跑）
if [[ -s "$OUTGFF" ]]; then
  echo "[SKIP] $OUTGFF"
else
  echo "[INFO] Host: $(hostname)"
  echo "[INFO] Input: $INFA"
  echo "[INFO] Out  : $OUTGFF"
  echo "[INFO] KINGDOM=$KINGDOM  Threads=${SLURM_CPUS_PER_TASK}"
  # 用环境里的 perl 调用 barrnap，避免系统 perl
  "$ENV/bin/perl" "$BARRNAP" --kingdom "$KINGDOM" --threads "${SLURM_CPUS_PER_TASK}" "$INFA" > "$OUTGFF"
fi

# === 统计 16S/23S/5S copy numbers 并写入汇总 CSV ===
SUMMARY="$OUTROOT/summary_16S_copy_number.csv"
LOCK="$SUMMARY.lock"

# 确保汇总表头（带锁）
{
  flock 200
  if [[ ! -f "$SUMMARY" ]]; then
    echo "sample,input,copy_16S,copy_23S,copy_5S" > "$SUMMARY"
  fi
} 200>"$LOCK"

# 解析 GFF 统计
count_16S=$(awk '$3=="rRNA" && $9 ~ /16S_rRNA/' "$OUTGFF" | wc -l | awk '{print $1}')
count_23S=$(awk '$3=="rRNA" && $9 ~ /23S_rRNA/' "$OUTGFF" | wc -l | awk '{print $1}')
count_5S=$(awk  '$3=="rRNA" && $9 ~ /5S_rRNA/'  "$OUTGFF" | wc -l | awk '{print $1}')

# 追加到汇总 CSV（带锁，避免并发写冲突）
{
  flock 200
  printf "%s,%s,%s,%s,%s\n" "$SAMPLE" "$INFA" "$count_16S" "$count_23S" "$count_5S" >> "$SUMMARY"
} 200>"$LOCK"

echo "[DONE] $SAMPLE | 16S=$count_16S 23S=$count_23S 5S=$count_5S"
