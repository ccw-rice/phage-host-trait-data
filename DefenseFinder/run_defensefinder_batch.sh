#!/bin/bash

module purge
module load Miniforge3/24.1.2-0
conda activate /scratch/$USER/defensefinder_env

export MACSY_DATA_DIR=/scratch/$USER/macsydata
export TMPDIR=/scratch/$USER/tmp_defensefinder
export HOME=/scratch/$USER

mkdir -p "$TMPDIR" df_out2

for f in /scratch/$USER/df_data/*.fna; do
    base=$(basename "$f" .fna)
    OUT=/scratch/$USER/df_out2/"${base}_$(date +%s)"
    
    defense-finder run "$f" \
        --models-dir "$MACSY_DATA_DIR" \
        --workers 8 \
        -o "$OUT"
done
