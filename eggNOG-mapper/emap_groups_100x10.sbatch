#!/bin/bash -l
#SBATCH -J emap_groups_restx10
#SBATCH -p commons
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH -o emap_groups.%A_%a.out
#SBATCH -e emap_groups.%A_%a.err
#SBATCH --array=0-5991%30     # 并发 30；索引写到 5991（若超出实际组数会优雅退出）

set -Eeuo pipefail
trap 'echo "[ERR] line:$LINENO cmd:$BASH_COMMAND" >&2' ERR

################ 0) 环境 ################
module purge
export MAMBA_ROOT_PREFIX=/scratch/cw106/software
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook --shell bash)"
micromamba activate /projects/alvarez/conda_envs/eggnog310_env

EMAPPER=/projects/alvarez/conda_envs/eggnog310_env/bin/emapper.py
export PATH="$(dirname "$EMAPPER"):$PATH"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
echo "[`date`] $(emapper.py --version | head -n1)"

################ 1) 参数与输入收集（加锁，最小化修改） ################
DATA=/scratch/cw106/eggnog_input/ncbi_dataset/data
OUTDIR=/scratch/cw106/eggnog_input/new_results
mkdir -p "$OUTDIR"

LIST=/scratch/cw106/eggnog_test/genomes_rest.list
LOCK="${LIST}.lock"

# 只生成一次“剩余清单”，其余任务复用，避免并发写入竞态
(
  flock -w 300 9 || { echo "[ERR] lock timeout on $LOCK"; exit 1; }
  if [ ! -s "$LIST" ]; then
    echo "[INFO] building rest list -> $LIST"
    TMP_ALL=$(mktemp)
    find "$DATA" -maxdepth 1 -type d -name 'GCF_*' | sort > "$TMP_ALL"
    : > "$LIST"
    while IFS= read -r d; do
      g=$(basename "$d")
      if [ ! -s "$OUTDIR/${g}.emapper.annotations" ]; then
        echo "$d" >> "$LIST"
      fi
    done < "$TMP_ALL"
    rm -f "$TMP_ALL"
  else
    echo "[INFO] using existing $LIST"
  fi
) 9>"$LOCK"

GROUP_SIZE=10

# 可见性：打印剩余条目与组数（便于你核对但不改变流程）
M=$(wc -l < "$LIST" || echo 0)
G=$(( (M + GROUP_SIZE - 1) / GROUP_SIZE ))
echo "[INFO] REST M=$M  =>  groups G=$G  (size=$GROUP_SIZE)"

# 计算本任务负责的区间（1-based）
start=$(( SLURM_ARRAY_TASK_ID * GROUP_SIZE + 1 ))
end=$(( start + GROUP_SIZE - 1 ))

# 取该组目录列表
mapfile -t GROUP_DIRS < <(sed -n "${start},${end}p" "$LIST" | grep -v '^\s*$' || true)

# 没有输入则优雅退出（当数组索引超出真实组数时会发生）
if [ ${#GROUP_DIRS[@]} -eq 0 ]; then
  echo "[INFO] No inputs for task $SLURM_ARRAY_TASK_ID (probably beyond group count). Exit."
  exit 0
fi

echo "[`date`] Task $SLURM_ARRAY_TASK_ID -> ${#GROUP_DIRS[@]} genomes"

################ 2) 本地工作目录与数据库缓存 ################
workdir=${SLURM_TMPDIR:-/scratch/$USER/tmp_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}}
mkdir -p "$workdir"

EGGNOG_DATA_SRC=/scratch/cw106/eggnog_data
for f in eggnog_proteins.dmnd eggnog.db eggnog.taxa.db eggnog.taxa.db.traverse.pkl; do
  cp -n "$EGGNOG_DATA_SRC/$f" "$workdir/"   # 已存在则不覆盖
done

################ 3) 合并本组 10 个基因组蛋白并加前缀 ################
merged="$workdir/group_${SLURM_ARRAY_TASK_ID}.faa"
: > "$merged"
for d in "${GROUP_DIRS[@]}"; do
  gcf=$(basename "$d")
  faa="$d/protein.faa"
  if [ ! -s "$faa" ]; then
    echo "[WARN] missing $faa"
    continue
  fi
  # 给每条标题加 GCF 前缀，便于后续按基因组拆回
  awk -v pre="${gcf}|" '/^>/{sub(/^>/,">"pre)}{print}' "$faa" >> "$merged"
done

nseq=$(grep -c '^>' "$merged" || true)
if [ "$nseq" -le 0 ]; then
  echo "[WARN] empty merged for task $SLURM_ARRAY_TASK_ID"
  exit 0
fi
echo "[`date`] merged ${#GROUP_DIRS[@]} genomes, $nseq proteins"

################ 4) 运行 eggNOG-mapper ################
# 速度优先：DIAMOND + fast + iterate + dbmem
"$EMAPPER" \
  -i "$merged" \
  -o "group_${SLURM_ARRAY_TASK_ID}" \
  -m diamond \
  --cpu "$SLURM_CPUS_PER_TASK" \
  --sensmode fast \
  --dmnd_iterate yes \
  --dbmem \
  --target_taxa 2 \
  --no_file_comments \
  --data_dir "$workdir" \
  --temp_dir "$workdir" \
  --output_dir "$OUTDIR" \
  --override

ann="$OUTDIR/group_${SLURM_ARRAY_TASK_ID}.emapper.annotations"
[ -s "$ann" ] || { echo "[ERR] missing annotations $ann"; exit 1; }

################ 5) 拆分注释表回各自 GCF ################
awk -F'\t' '
  NR==1 { hdr=$0; next }
  {
    split($1, a, "|"); g=a[1]         # query 形如 GCF_xxx|protein_id
    if (!(g in seen)) { print hdr > (dir "/" g ".emapper.annotations"); seen[g]=1 }
    print >> (dir "/" g ".emapper.annotations")
  }
' dir="$OUTDIR" "$ann"

echo "[`date`] ✓ task $SLURM_ARRAY_TASK_ID done → $OUTDIR"
