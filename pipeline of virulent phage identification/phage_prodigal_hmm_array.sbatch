#!/bin/bash
#SBATCH -J phage_pHMM
#SBATCH -o /scratch/cw106/logs/phage_pHMM.%A_%a.out
#SBATCH -e /scratch/cw106/logs/phage_pHMM.%A_%a.err
#SBATCH --partition=commons
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=24:00:00
#SBATCH --array=0-0

set -Eeuo pipefail
umask 022
trap 'echo "[ERR] line:$LINENO cmd:$BASH_COMMAND" >&2' ERR

# ====== inputs you care about ======
LIST_TSV="${LIST_TSV:-/scratch/cw106/phage_id_list_2col.ALL15k.tsv}"
DATA_ROOT="${DATA_ROOT:-/scratch/cw106/phagescope_data_resource}"
HMM_DB="${HMM_DB:-/projects/alvarez/virus_hmm_db/merged_phage.hmm}"

OUT_ROOT="${OUT_ROOT:-/scratch/cw106/phage_hmm_out_all15k}"
WORK_ROOT="${WORK_ROOT:-/scratch/cw106/phage_hmm_work_all15k}"

CHUNK="${CHUNK:-40}"   # rows per array task (excluding header)

mkdir -p "$OUT_ROOT" "$WORK_ROOT" /scratch/cw106/logs

# ====== env (reuse your working style) ======
module purge || true
module use /opt/apps/modules/Core
module --ignore_cache load Miniforge3/24.11.3-0
eval "$(conda shell.bash hook)"
conda activate phispy_env   # has prodigal + hmmsearch in your earlier pipeline

command -v prodigal
command -v hmmsearch

# ====== compute slice ======
N=$(($(wc -l < "$LIST_TSV") - 1))
start=$(( SLURM_ARRAY_TASK_ID * CHUNK + 2 ))
end=$(( start + CHUNK - 1 ))
(( end > N+1 )) && end=$((N+1))

echo "[INFO] LIST=$LIST_TSV N=$N task=$SLURM_ARRAY_TASK_ID rows=$start..$end"

# ====== process each row ======
for line_no in $(seq "$start" "$end"); do
  row=$(sed -n "${line_no}p" "$LIST_TSV" || true)
  [[ -z "${row:-}" ]] && continue
  src=$(echo "$row" | cut -f1)
  pid=$(echo "$row" | cut -f2)

  # resolve fasta path: DATA_ROOT/<src>/<src>/<pid>.fasta
  fa="$DATA_ROOT/$src/$src/${pid}.fasta"
  [[ -s "$fa" ]] || fa="$DATA_ROOT/$src/$src/${pid}.fa"
  [[ -s "$fa" ]] || fa="$DATA_ROOT/$src/$src/${pid}.fna"
  [[ -s "$fa" ]] || fa="$DATA_ROOT/$src/$src/${pid}.fasta.gz"
  [[ -s "$fa" ]] || fa="$DATA_ROOT/$src/$src/${pid}.fa.gz"
  [[ -s "$fa" ]] || fa="$DATA_ROOT/$src/$src/${pid}.fna.gz"

  if [[ ! -s "$fa" ]]; then
    echo "[MISS] $src $pid  (no fasta found under $DATA_ROOT/$src/$src/)" >&2
    continue
  fi

  # output locations (one dir per source)
  out_dir="$OUT_ROOT/$src"
  mkdir -p "$out_dir"

  faa="$out_dir/$pid.faa"
  tbl="$out_dir/$pid.hmm.tbl"

  # skip if already done
  if [[ -s "$faa" && -s "$tbl" ]]; then
    echo "[SKIP] $src $pid"
    continue
  fi

  # scratch temp
  scr="/scratch/$USER/$SLURM_JOB_ID/${src}_${pid}"
  mkdir -p "$scr"
  cd "$scr"

  # copy/decompress fasta
  if [[ "$fa" == *.gz ]]; then
    zcat "$fa" > phage.fna
  else
    cp "$fa" phage.fna
  fi

  echo "[RUN ] $src $pid prodigal"
  prodigal -i phage.fna -a phage.faa -p meta -q

  echo "[RUN ] $src $pid hmmsearch"
  hmmsearch --cpu "$SLURM_CPUS_PER_TASK" -E 1e-3 \
    --tblout hmm.tbl \
    "$HMM_DB" phage.faa > /dev/null

  # sync back
  rsync -av phage.faa "$faa"
  rsync -av hmm.tbl   "$tbl"

  cd /
  rm -rf "$scr"
  echo "[DONE] $src $pid"
done

echo "[ALL DONE]"
