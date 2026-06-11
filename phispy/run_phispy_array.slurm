#!/bin/bash
#SBATCH --job-name=phispy_array
#SBATCH --partition=commons
#SBATCH --account=commons
#SBATCH --array=0-1529%30                  # Run first 8 jobs for testing; change to 0-5752%30 for full production run
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=03:30:00
#SBATCH -o /projects/alvarez/logs/phispy_%A_%a.out
#SBATCH -e /projects/alvarez/logs/phispy_%A_%a.err

############################ 1. Environment ############################
module purge

# Critical: allow compute nodes to locate the Core module directory
module use /opt/apps/modules/Core

# Load the available Miniforge3 version
module --ignore_cache load Miniforge3/24.11.3-0

eval "$(conda shell.bash hook)"

# PhiSpy should already be installed via:
# conda install -c bioconda phispy
conda activate phispy_env

set -euo pipefail

# --- Optional debug output (can be removed after validation) ---
echo "[DEBUG] host: $(hostname)"
echo "[DEBUG] loaded modules:"
module list
echo "-----------------------------------------"

########################### 2. Scratch directories ###########################
SCR=/scratch/$USER/$SLURM_JOB_ID/$SLURM_ARRAY_TASK_ID

mkdir -p "$SCR"

INPUT_DIR=/scratch/cw106/phispy_input

OUTPUT_DIR=$SCR/phispy_out
mkdir -p "$OUTPUT_DIR"

DEST_BASE=/scratch/cw106/phispy_output
mkdir -p "$DEST_BASE"

##################### 3. Select the GBFF file for this array task #####################
files=("$INPUT_DIR"/*.gbff)

gbff=${files[$SLURM_ARRAY_TASK_ID]}

id=$(basename "${gbff%.*}")   # Example: GCF_xxx_genomic

LIST_FILE=${LIST_FILE:-/projects/alvarez/genomes.list}

GCF=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$LIST_FILE")

[[ "$id" != "$GCF" ]] && \
    echo "[WARN] Genome list and filename do not match: $id vs $GCF"

################## 4. Copy to local scratch and run PhiSpy ##################
cp "$gbff" "$SCR/"

local_gbff=$SCR/$(basename "$gbff")

################## 4. Run PhiSpy ##################
PhiSpy.py \
    -o "$OUTPUT_DIR/$id" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --output_choice 7 \
    --keep_dropped_predictions \
    "$local_gbff"

################### 5. Sync results and clean up ###################
rsync -av "$OUTPUT_DIR/$id/" "$DEST_BASE/$id/"

rm -rf "$SCR"

echo "[DONE] $id completed (results stored in $DEST_BASE/$id)"
