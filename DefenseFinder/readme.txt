# DefenseFinder pipeline for anti-phage defense system detection

## Overview

This pipeline identifies anti-phage defense systems in bacterial genomes using DefenseFinder and summarizes the results into genome-level tables.

DefenseFinder detects known anti-phage defense systems such as:

* Restriction-modification (RM)
* CBASS
* BREX
* DISARM
* CRISPR-associated systems
* and other defense modules

A custom script (`merge_defense_systems.py`) is used to merge and summarize DefenseFinder output across all genomes.

---

## Software requirements

DefenseFinder
https://github.com/mdmparis/defense-finder

Python ≥ 3.8

Required Python package:

```
pandas
```

---

## Input data

Input genomes must be nucleotide FASTA files:

```
*.fna
```

Example directory:

```
/scratch/cw106/df_data/

    GCF_000008625.1.fna
    GCF_000005845.2.fna
    ...
```

Each file represents one genome.

---

## Running DefenseFinder

DefenseFinder was executed using the following command:

```
defense-finder run <genome.fna> \
    --models-dir <MACSY_DATA_DIR> \
    --workers 8 \
    -o <output_directory>
```

Example batch execution:

```
module purge
module load Miniforge3/24.1.2-0
conda activate /scratch/$USER/defensefinder_env

export MACSY_DATA_DIR=/scratch/$USER/macsydata

for f in /scratch/$USER/df_data/*.fna; do

    base=$(basename "$f" .fna)

    OUT=/scratch/$USER/df_out2/${base}_$(date +%s)

    defense-finder run "$f" \
        --models-dir "$MACSY_DATA_DIR" \
        --workers 8 \
        -o "$OUT"

done
```

---

## DefenseFinder output structure

Each genome produces an output directory containing multiple files, including:

```
*_defense_finder_systems.tsv
*_defense_finder_genes.tsv
*_defense_finder_summary.tsv
```

This pipeline uses:

```
*_defense_finder_systems.tsv
```

which contains detected defense system annotations.

---

## Merging results across genomes

Script used:

```
merge_defense_systems.py
```

This script searches recursively under:

```
/scratch/cw106/df_out2/
```

for all files matching:

```
*_defense_finder_systems.tsv
```

---

## Output files

### 1. all_genomes_defense_systems.tsv

Contains all detected defense systems across all genomes.

Example columns:

```
GCF_ID
type
subtype
genes_count
proteins_in_system
sys_beg
sys_end
activity
source_file
```

Each row represents one detected defense system.

---

### 2. systems_count_by_type.tsv

Matrix summarizing defense system counts per genome.

Example:

```
GCF_ID    RM    CBASS    BREX    DISARM
GCF_001   3     0        1       2
```

Each value represents the number of detected systems of each type in that genome.

---

## Running the merge script

Execute:

```
python merge_defense_systems.py
```

The script automatically:

* Finds all DefenseFinder result files
* Merges annotations
* Generates summary tables

Output location:

```
/scratch/cw106/df_out2/

    all_genomes_defense_systems.tsv
    systems_count_by_type.tsv
```

---

## Script description

merge_defense_systems.py performs the following steps:

1. Recursively locates all DefenseFinder systems files
2. Standardizes genome identifiers (GCF_ID)
3. Merges all detected systems into one table
4. Generates genome-level count matrices

No filtering or modification of DefenseFinder results is performed.

---

## Files included in this repository

```
scripts/
    merge_defense_systems.py

README_defensefinder.md
```

---

## Reproducibility

All defense systems were detected using DefenseFinder with default models.

This script aggregates DefenseFinder outputs without modifying individual system annotations.

---

## Contact

For questions regarding this pipeline, please contact the authors.
