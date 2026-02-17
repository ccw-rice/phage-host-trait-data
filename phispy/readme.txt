# Prophage identification pipeline (PhiSpy + secondary screening)

This directory contains scripts and marker lists used for prophage identification and secondary quality filtering in:

Wu et al. Host life-history strategies structure phage lifestyle associations across bacterial species.

All scripts provided here were used to generate the prophage dataset analyzed in the manuscript.

---

# Overview

Prophage regions were first predicted using PhiSpy and subsequently subjected to secondary screening based on hallmark phage genes and insertion sequence (IS) content.

The pipeline was executed on the Rice University NOTSX high-performance computing cluster using Slurm job arrays. Each genome was processed independently, and filtering criteria were applied in a reproducible and automated manner.

---

# Pipeline steps

## Step 1. Prophage prediction using PhiSpy

Script:

run_phispy_array.slurm


Input:

- GenBank genome files (.gbff)

Output:

- Candidate prophage sequences:
phage.fasta


Software:

- PhiSpy v4.2.21

---

## Step 2. Gene prediction using Prodigal

Script:

prodigal_array.slurm


Input:

phage.fasta


Output:

phage.faa
phage.gff


Software:

- Prodigal v2.6.3

---

## Step 3. Functional annotation using HMMER

Script:

hmm_array.slurm


Input:

phage.faa


Output:

hmm.tbl


Software:

- HMMER v3.4

Database:

- merged PHROG and VOG HMM profile database

---

## Step 4. Marker-based prophage filtering

Scripts:

parse_tblout_filter.py
filter_phage.slurm


Marker lists:

integrase_unified.list
structural_unified.list


Filtering criteria:

Candidate prophage regions were evaluated based on the presence of hallmark phage genes:

- integrase
- structural proteins (capsid, tail, portal, terminase, etc.)

Output:

filter_summary.tsv


---

## Step 5. Insertion sequence detection using ISEScan

Script:

run_isescan_array.slurm


Software:

- ISEScan v1.7.3

Output:

ISEScan annotation files used to calculate insertion sequence overlap.

---

## Step 6. IS content calculation and filtering

Script:

add_is_len_filter_v2.py


This script computes the following metrics for each candidate prophage region:

- prophage length
- IS overlap length (IS_bp)
- IS percentage (IS_pct)

The script includes filtering logic based on:

- prophage length ≥ 18 kb
- IS content < 25%

IMPORTANT NOTE:

Although the script includes a minimum length threshold of 18 kb, this length filter was NOT applied in the final dataset used in the manuscript. All prophage regions were retained regardless of length, provided they satisfied structural gene and IS-content criteria.

The length filter was included for conservative screening during pipeline development but was not enforced in the final analysis.

Output:

*_islen.tsv


---

## Step 7. Final summary generation

Script:

merge_filter_summary.sh


This script aggregates prophage filtering results across all genomes and produces:

all_filter_summary_islen.tsv


This file contains the final set of high-confidence prophage regions used in downstream analyses.

---

# Marker lists

The following marker lists define hallmark phage genes used for prophage identification:

integrase_unified.list
structural_unified.list


These lists were constructed based on annotations from PHROG and VOG databases.

---

# Computational environment

All analyses were performed on:

Rice University NOTSX high-performance computing cluster

Software used:

- PhiSpy v4.2.21
- Prodigal v2.6.3
- HMMER v3.4
- ISEScan v1.7.3
- Python 3
- Bash
- Slurm workload manager

---

# Reproducibility

All scripts required to reproduce prophage identification and filtering are provided in this directory.

External databases required:

- PHROG database
- VOG database
- merged HMM database used for annotation

---

# Contact

Chuncheng Wu  
Rice University  