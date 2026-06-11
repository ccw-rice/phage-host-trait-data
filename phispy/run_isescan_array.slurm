#!/bin/bash
#SBATCH -J isescan_default
#SBATCH --partition=commons
#SBATCH --array=585,4358,5556,5668
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH -o /projects/alvarez/logs/isescan_%A_%a.out
#SBATCH -e /projects/alvarez/logs/isescan_%A_%a.err

#################### 1. Conda environment ####################
module purge
module use /opt/apps/modules/Core
module --ignore_cache load Miniforge3/24.11.3-0

eval "$(conda shell.bash hook)"
conda activate isescan_env

export HMMER_NCPU=$SLURM_CPUS_PER_TASK
export PYTHONNOUSERSITE=1

set -euo pipefail

#################### 2. Genome ID ####################
GCF=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" /projects/alvarez/genomes.list)

#################### 3. Working directory on /scratch ####################
SCR=/scratch/$USER/$SLURM_JOB_ID/$GCF

mkdir -p "$SCR"

# Copy input FASTA from /scratch/cw106/phispy_input
cp /scratch/cw106/phispy_input/${GCF}_genomic.fna "$SCR/"

#################### 4. Run ISEScan ####################
isescan.py \
    --seqfile "$SCR/${GCF}_genomic.fna" \
    --output "$SCR/isescan_out" \
    --nthread "$SLURM_CPUS_PER_TASK"

#################### 5. Sync results back to /scratch and clean up ####################
# Store results in the designated persistent output directory
DEST=/scratch/cw106/isescan_output/$GCF

mkdir -p "$DEST"

rsync -av "$SCR/isescan_out/" "$DEST/"

rm -rf "$SCR"

echo "[DONE] ISEScan finished for $GCF (results located in $DEST)"
