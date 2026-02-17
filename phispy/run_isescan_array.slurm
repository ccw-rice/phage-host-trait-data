#!/bin/bash
#SBATCH -J isescan_default
#SBATCH --partition=commons
#SBATCH --array=585,4358,5556,5668
#SBATCH --cpus-per-task=8          # ← 用长选项
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH -o /projects/alvarez/logs/isescan_%A_%a.out
#SBATCH -e /projects/alvarez/logs/isescan_%A_%a.err

#################### 1. Conda 环境 ####################
module purge
module use /opt/apps/modules/Core
module --ignore_cache load Miniforge3/24.11.3-0
eval "$(conda shell.bash hook)"
conda activate isescan_env
export HMMER_NCPU=$SLURM_CPUS_PER_TASK
export PYTHONNOUSERSITE=1
set -euo pipefail
#################### 2. 样本 ID ####################
GCF=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" /projects/alvarez/genomes.list)

#################### 3. /scratch 工作目录 ####################
SCR=/scratch/$USER/$SLURM_JOB_ID/$GCF
mkdir -p "$SCR"

# ★ 输入 fasta 从 /scratch/cw106/phispy_input 复制
cp /scratch/cw106/phispy_input/${GCF}_genomic.fna "$SCR/"

#################### 4. 运行 ISEScan ####################
isescan.py \
    --seqfile  "$SCR/${GCF}_genomic.fna" \
    --output   "$SCR/isescan_out" \
    --nthread  "$SLURM_CPUS_PER_TASK"

#################### 5. 同步结果回 /scratch 并清理 ####################
# ★ 把结果放进你准备好的持久目录
DEST=/scratch/cw106/isescan_output/$GCF
mkdir -p "$DEST"
rsync -av "$SCR/isescan_out/" "$DEST/"

rm -rf "$SCR"
echo "[DONE] ISEScan finished for $GCF   (结果位于 $DEST)"
