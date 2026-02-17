# 16S rRNA copy number prediction pipeline using Barrnap

This directory contains the script used to identify rRNA genes and calculate rRNA copy numbers from bacterial genomes using Barrnap.

This pipeline was used in:

Wu et al. Host life-history strategies structure phage lifestyle associations across bacterial species.

---

# Overview

Barrnap (v0.9) was used to predict rRNA genes in bacterial genomes. The number of 16S, 23S, and 5S rRNA genes was extracted from Barrnap output and summarized for each genome.

These values were used as genome-level traits in downstream host life-history strategy analyses.

All computations were performed on the Rice University NOTSX high-performance computing cluster using Slurm job arrays.

---

# Input data

Genome FASTA files in one of the following formats:

.fna
.fa
.fasta


Manifest file format:

Each line contains the absolute path to one genome FASTA file:

/path/to/genome1.fna
/path/to/genome2.fna
/path/to/genome3.fna


---

# Pipeline execution

Script:

barrnap_batch.v3.sbatch


Example usage:

MANIFEST=/scratch/cw106/barrnap_manifest.txt

N=$(wc -l < "$MANIFEST")

sbatch --array=1-$N barrnap_batch.v3.sbatch
"$MANIFEST"
/scratch/cw106/barrnap_out
bac


Arguments:

Argument 1: manifest file containing genome paths
Argument 2: output directory
Argument 3: kingdom (bac, arc, or euk)


Default kingdom used in this study:

bac


---

# Processing steps

For each genome, the script performs the following steps:

1. Run Barrnap to identify rRNA genes:

barrnap --kingdom bac genome.fna


2. Output rRNA annotation in GFF format:

genome_rRNA.gff


3. Count number of rRNA genes:

- 16S rRNA copy number
- 23S rRNA copy number
- 5S rRNA copy number

4. Append results to summary table.

File locking is used to ensure safe parallel writing during Slurm array execution.

---

# Output structure

Output directory structure:

output_directory/

sample1/
    genome_rRNA.gff

sample2/
    genome_rRNA.gff

summary_16S_copy_number.csv

---

# Summary output file

summary_16S_copy_number.csv


Columns:

sample
input
copy_16S
copy_23S
copy_5S


These values represent the predicted number of rRNA genes in each genome.

The 16S rRNA copy number was used in downstream analyses.

---

# Software used

Barrnap v0.9

Barrnap repository:

https://github.com/tseemann/barrnap

---

# Computational environment

Rice University NOTSX high-performance computing cluster

Slurm workload manager

Micromamba environment containing Barrnap and dependencies.

---

# Reproducibility

This repository provides the complete script required to reproduce rRNA copy number prediction.

Users must install Barrnap and provide genome FASTA files.

---

# Contact

Chuncheng Wu  
Rice University  