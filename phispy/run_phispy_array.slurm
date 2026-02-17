#!/bin/bash
#SBATCH --job-name=phispy_array
#SBATCH --partition=commons
#SBATCH --account=commons
#SBATCH --array=0-1529%30                  # ★ 先跑 8 条；正式改 0-5752%30
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=03:30:00
#SBATCH -o /projects/alvarez/logs/phispy_%A_%a.out
#SBATCH -e /projects/alvarez/logs/phispy_%A_%a.err

############################ 1. 环境 ############################
module purge
module use /opt/apps/modules/Core                      # ★ 关键：让 compute 节点找得到 Core 模块目录
module --ignore_cache load Miniforge3/24.11.3-0        # ★ 指定可用版本
eval "$(conda shell.bash hook)"
conda activate phispy_env                              # 预先 conda install -c bioconda phispy
set -euo pipefail

# --- 可选调试输出（跑通后可删） ---
echo "[DEBUG] host: $(hostname)"
echo "[DEBUG] loaded modules:"; module list
echo "-----------------------------------------"

########################### 2. scratch 目录 ###########################
SCR=/scratch/$USER/$SLURM_JOB_ID/$SLURM_ARRAY_TASK_ID
mkdir -p "$SCR"

INPUT_DIR=/scratch/cw106/phispy_input
OUTPUT_DIR=$SCR/phispy_out
mkdir -p "$OUTPUT_DIR"

DEST_BASE=/scratch/cw106/phispy_output
mkdir -p "$DEST_BASE"

##################### 3. 选取本子任务的 gbff 文件 #####################
files=("$INPUT_DIR"/*.gbff)
gbff=${files[$SLURM_ARRAY_TASK_ID]}
id=$(basename "${gbff%.*}")                   # 例 GCF_..._genomic

LIST_FILE=${LIST_FILE:-/projects/alvarez/genomes.list}
GCF=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$LIST_FILE")
[[ "$id" != "$GCF" ]] && echo "[WARN] list 与文件名不符: $id vs $GCF"

################## 4. 复制到本地盘并运行 PhiSpy ##################
cp "$gbff" "$SCR/"
local_gbff=$SCR/$(basename "$gbff")

################## 4. 运行 PhiSpy ##################
PhiSpy.py \
        -o "$OUTPUT_DIR/$id" \
        --threads "$SLURM_CPUS_PER_TASK" \
        --output_choice 7 \
        --keep_dropped_predictions \
        "$local_gbff"

################### 5. 同步结果并清理 ###################
rsync -av "$OUTPUT_DIR/$id/" "$DEST_BASE/$id/"
rm -rf "$SCR"
echo "[DONE] $id 完成 (结果存于 $DEST_BASE/$id)"
