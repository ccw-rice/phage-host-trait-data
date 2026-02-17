#!/bin/bash
#SBATCH -J phage_lifestyle_call
#SBATCH --partition=commons
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 00:45:00
#SBATCH -o /projects/alvarez/logs/lifestyle_call.%j.out
#SBATCH -e /projects/alvarez/logs/lifestyle_call.%j.err

set -euo pipefail
echo "[START] $(date)"

# ────────── 1. 环境（与你之前所有成功任务一致：micromamba） ──────────
module purge
module use /opt/apps/modules/Core
module --ignore_cache load Miniforge3/24.11.3-0

export MAMBA_ROOT_PREFIX=/scratch/cw106/software
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook --shell bash)"
micromamba activate phatyp_v1

echo "[INFO] Python:"
which python
python -V

# ────────── 2. 输入 / 输出路径 ──────────
export LIST_TSV=/scratch/cw106/phage_id_list_2col.ALL15k.tsv

# prodigal + hmmsearch 的输出根目录
export OUT_ROOT=/scratch/cw106/phage_hmm_out_all15k

# 原始 phage genome 所在目录（用于算 genome length）
export DATA_ROOT=/scratch/cw106/phagescope_data_resource

# HMM / list 数据库
export DB_DIR=/projects/alvarez/virus_hmm_db
export INTEG_LIST=$DB_DIR/integrase_unified.list
export TERL_LIST=$DB_DIR/terL.list
export LYSIS_LIST=$DB_DIR/lysis_genes.list

# 最终输出表
export OUT_TSV=/scratch/cw106/phage_lifestyle_calls_all15k.tsv

# ────────── 3. 基本存在性检查（快速 fail，避免 silent error） ──────────
for f in \
  "$LIST_TSV" \
  "$INTEG_LIST" \
  "$TERL_LIST" \
  "$LYSIS_LIST"
do
  if [[ ! -s "$f" ]]; then
    echo "[ERR] missing or empty: $f" >&2
    exit 2
  fi
done

# ────────── 4. 运行判定脚本 ──────────
python /projects/alvarez/scripts/call_lifestyle_from_tbl.py

echo "[DONE] $(date)"
