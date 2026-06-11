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

# ────────── 1. Environment (same micromamba setup used by previous successful jobs) ──────────
module purge
module use /opt/apps/modules/Core
module --ignore_cache load Miniforge3/24.11.3-0

export MAMBA_ROOT_PREFIX=/scratch/cw106/software

eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook --shell bash)"

micromamba activate phatyp_v1

echo "[INFO] Python:"
which python
python -V

# ────────── 2. Input and output paths ──────────
export LIST_TSV=/scratch/cw106/phage_id_list_2col.ALL15k.tsv

# Output root directory from prodigal + hmmsearch
export OUT_ROOT=/scratch/cw106/phage_hmm_out_all15k

# Original phage genome directory
# (used to calculate genome length)
export DATA_ROOT=/scratch/cw106/phagescope_data_resource

# HMM and gene-list databases
export DB_DIR=/projects/alvarez/virus_hmm_db
export INTEG_LIST=$DB_DIR/integrase_unified.list
export TERL_LIST=$DB_DIR/terL.list
export LYSIS_LIST=$DB_DIR/lysis_genes.list

# Final output table
export OUT_TSV=/scratch/cw106/phage_lifestyle_calls_all15k.tsv

# ────────── 3. Basic file existence checks
#             (fail fast to avoid silent errors) ──────────
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

# ────────── 4. Run lifestyle classification script ──────────
python /projects/alvarez/scripts/call_lifestyle_from_tbl.py

echo "[DONE] $(date)"
