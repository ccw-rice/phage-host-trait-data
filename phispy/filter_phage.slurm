#!/bin/bash
#SBATCH -J prophage_filter
#SBATCH --array=5753
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 01:00:00
#SBATCH -o /projects/alvarez/logs/filter_%A_%a.out
#SBATCH -e /projects/alvarez/logs/filter_%A_%a.err

# ────────── 1. 环境 ──────────
module purge
module use /opt/apps/modules/Core          # Core 目录里有 Miniforge3 的 modulefile
module --ignore_cache load Miniforge3/24.11.3-0
eval "$(conda shell.bash hook)"
conda activate phispy_env
set -euo pipefail

# —— 2. 样本 ID ——
GENOME=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" /projects/alvarez/genomes.list)

# —— 3. 源目录指向 /scratch —— 
dir="/scratch/cw106/phispy_output/${GENOME}"
[[ -d "${dir}_genomic" ]] && dir="${dir}_genomic"

TBL_SRC="${dir}/hmm.tbl"
SUM_DST="${dir}/filter_summary.tsv"
[[ ! -s "$TBL_SRC" ]] && { echo "[SKIP] 缺少 $TBL_SRC"; exit 0; }

# ========= 4. /scratch 工作目录 =========
SCR=/scratch/$USER/$SLURM_JOB_ID/$GENOME
mkdir -p "$SCR"
cp "$TBL_SRC" "$SCR/"

# —— 5. 运行过滤脚本（在 /scratch） ——
python /projects/alvarez/scripts/parse_tblout_filter.py \
       "$SCR/hmm.tbl" "$SCR/filter_summary.tsv"

# ========= 6. 同步结果回 /scratch 并清理 =========
rsync -av "$SCR/filter_summary.tsv" "$SUM_DST"
rm -rf "$SCR"
