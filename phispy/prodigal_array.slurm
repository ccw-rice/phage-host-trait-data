#!/bin/bash
#SBATCH --job-name=prodigal_faa
#SBATCH --partition=commons
#SBATCH --array=5753
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH -o /projects/alvarez/logs/prodigal_%A_%a.out
#SBATCH -e /projects/alvarez/logs/prodigal_%A_%a.err

# ─────────── 1. Conda environment ───────────
module purge
module use /opt/apps/modules/Core          # Parent directory containing the Miniforge3 modulefile
module --ignore_cache load Miniforge3/24.11.3-0

eval "$(conda shell.bash hook)"
conda activate phispy_env

set -euo pipefail

# ─────────── 2. Genome ID ───────────
GENOME_ID=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" /projects/alvarez/genomes.list)

[[ -d /scratch/cw106/phispy_output/${GENOME_ID}_genomic ]] && \
    GENOME_ID="${GENOME_ID}_genomic"

# ─────────── 3. Source and destination directories under /scratch ───────────
SRC_DIR=/scratch/cw106/phispy_output/${GENOME_ID}

fasta="$SRC_DIR/phage.fasta"
faa="$SRC_DIR/phage.faa"
gff="$SRC_DIR/phage.gff"

# Skip if results already exist
if [[ -s "$faa" && -s "$gff" ]]; then
    echo "[SKIP] ${GENOME_ID} already has phage.faa and phage.gff"
    exit 0
fi

# ─────────── 4. Temporary working directory on /scratch ───────────
SCR=/scratch/$USER/$SLURM_JOB_ID/$GENOME_ID

mkdir -p "$SCR"

cp "$fasta" "$SCR/"

cd "$SCR"

# ─────────── 5. Run Prodigal ───────────
echo "[RUN ] Prodigal on $GENOME_ID"

prodigal \
    -i phage.fasta \
    -a phage.faa \
    -o phage.gff \
    -p meta \
    -q

echo "[DONE] Generated phage.faa and phage.gff"

# ─────────── 6. Sync results back to /scratch/phispy_output and clean up ───────────
rsync -av phage.faa phage.gff "$SRC_DIR/"

cd /

rm -rf "$SCR"
