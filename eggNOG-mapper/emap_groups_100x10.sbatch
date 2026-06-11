#!/bin/bash -l
#SBATCH -J emap_groups_restx10
#SBATCH -p commons
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH -o emap_groups.%A_%a.out
#SBATCH -e emap_groups.%A_%a.err
#SBATCH --array=0-5991%30     # Concurrency = 30; array index set to 5991 (tasks beyond the actual group count will exit gracefully)

set -Eeuo pipefail
trap 'echo "[ERR] line:$LINENO cmd:$BASH_COMMAND" >&2' ERR

################ 0) Environment ################
module purge
export MAMBA_ROOT_PREFIX=/scratch/cw106/software
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook --shell bash)"
micromamba activate /projects/alvarez/conda_envs/eggnog310_env

EMAPPER=/projects/alvarez/conda_envs/eggnog310_env/bin/emapper.py
export PATH="$(dirname "$EMAPPER"):$PATH"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
echo "[`date`] $(emapper.py --version | head -n1)"

################ 1) Parameters and input collection (locked, minimal modification) ################
DATA=/scratch/cw106/eggnog_input/ncbi_dataset/data
OUTDIR=/scratch/cw106/eggnog_input/new_results
mkdir -p "$OUTDIR"

LIST=/scratch/cw106/eggnog_test/genomes_rest.list
LOCK="${LIST}.lock"

# Generate the remaining-genome list only once;
# all other tasks reuse it to avoid race conditions during concurrent writes
(
  flock -w 300 9 || { echo "[ERR] lock timeout on $LOCK"; exit 1; }

  if [ ! -s "$LIST" ]; then
    echo "[INFO] Building remaining-genome list -> $LIST"

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
    echo "[INFO] Using existing $LIST"
  fi
) 9>"$LOCK"

GROUP_SIZE=10

# Report remaining entries and group count for verification
M=$(wc -l < "$LIST" || echo 0)
G=$(( (M + GROUP_SIZE - 1) / GROUP_SIZE ))
echo "[INFO] REMAINING M=$M  =>  groups G=$G  (size=$GROUP_SIZE)"

# Compute the range assigned to this task (1-based indexing)
start=$(( SLURM_ARRAY_TASK_ID * GROUP_SIZE + 1 ))
end=$(( start + GROUP_SIZE - 1 ))

# Retrieve directories assigned to this group
mapfile -t GROUP_DIRS < <(
  sed -n "${start},${end}p" "$LIST" | grep -v '^\s*$' || true
)

# Exit gracefully if no inputs are assigned
# (this occurs when the array index exceeds the actual group count)
if [ ${#GROUP_DIRS[@]} -eq 0 ]; then
  echo "[INFO] No inputs for task $SLURM_ARRAY_TASK_ID (probably beyond group count). Exit."
  exit 0
fi

echo "[`date`] Task $SLURM_ARRAY_TASK_ID -> ${#GROUP_DIRS[@]} genomes"

################ 2) Local working directory and database cache ################
workdir=${SLURM_TMPDIR:-/scratch/$USER/tmp_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}}
mkdir -p "$workdir"

EGGNOG_DATA_SRC=/scratch/cw106/eggnog_data

for f in eggnog_proteins.dmnd eggnog.db eggnog.taxa.db eggnog.taxa.db.traverse.pkl; do
  cp -n "$EGGNOG_DATA_SRC/$f" "$workdir/"   # Do not overwrite if already present
done

################ 3) Merge proteins from 10 genomes and add genome prefixes ################
merged="$workdir/group_${SLURM_ARRAY_TASK_ID}.faa"

: > "$merged"

for d in "${GROUP_DIRS[@]}"; do

  gcf=$(basename "$d")
  faa="$d/protein.faa"

  if [ ! -s "$faa" ]; then
    echo "[WARN] Missing $faa"
    continue
  fi

  # Add GCF prefix to each sequence header
  # to facilitate splitting annotations back by genome later
  awk -v pre="${gcf}|" '/^>/{sub(/^>/,">"pre)}{print}' "$faa" >> "$merged"

done

nseq=$(grep -c '^>' "$merged" || true)

if [ "$nseq" -le 0 ]; then
  echo "[WARN] Empty merged file for task $SLURM_ARRAY_TASK_ID"
  exit 0
fi

echo "[`date`] Merged ${#GROUP_DIRS[@]} genomes, $nseq proteins"

################ 4) Run eggNOG-mapper ################
# Speed-optimized settings: DIAMOND + fast + iterate + dbmem
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

[ -s "$ann" ] || {
  echo "[ERR] Missing annotations file: $ann"
  exit 1
}

################ 5) Split annotation table back into individual GCF files ################
awk -F'\t' '
  NR==1 { hdr=$0; next }
  {
    split($1, a, "|")
    g=a[1]         # query format: GCF_xxx|protein_id

    if (!(g in seen)) {
      print hdr > (dir "/" g ".emapper.annotations")
      seen[g]=1
    }

    print >> (dir "/" g ".emapper.annotations")
  }
' dir="$OUTDIR" "$ann"

echo "[`date`] Task $SLURM_ARRAY_TASK_ID completed successfully → $OUTDIR"
