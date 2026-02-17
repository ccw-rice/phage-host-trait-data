# gRodon2 Growth Rate Prediction Pipeline

This pipeline predicts bacterial maximal growth potential and doubling time using gRodon v2 based on coding sequence (CDS) nucleotide data.

All analyses were performed on the Rice University NOTSX high-performance computing cluster using Slurm.

---

## Software Requirements

The pipeline requires the following software:

R (tested with gRodon v2.5.1)

Required R packages:

gRodon  
Biostrings  

Slurm workload manager

The script uses an explicit R environment:

/scratch/$USER/envs/grodon2_env/bin/Rscript

and library path:

/scratch/$USER/envs/grodon2_env/lib/R/library

No manual environment activation is required.

---

## Input Data

Each genome must provide:

cds_from_genomic.fna

This file contains CDS nucleotide sequences.

Example directory structure:

/scratch/cw106/gRodon2_test/

    GCF_xxxxxxxx.x/
        cds_from_genomic.fna

---

## Manifest File

A manifest file is required. It must contain one absolute path per line:

Example:

/scratch/cw106/gRodon2_test/GCF_000001405.40/cds_from_genomic.fna  
/scratch/cw106/gRodon2_test/GCF_000005845.2/cds_from_genomic.fna  

Example command to generate manifest:

find /scratch/cw106/gRodon2_test \
-type f -name "cds_from_genomic.fna" \
| sort > /scratch/cw106/grodon2_manifest.txt

---

## Running the Pipeline

Script:

grodon2_batch.v2.sbatch

Submit using Slurm array:

N=$(wc -l < /scratch/cw106/grodon2_manifest.txt)

sbatch --array=1-$N \
/projects/alvarez/scripts/grodon2_batch.v2.sbatch \
/scratch/cw106/grodon2_manifest.txt \
/scratch/cw106/grodon2_out

Arguments:

Argument 1: manifest file (required)

Argument 2: output directory (optional, default: /scratch/$USER/grodon2_out)

---

## Pipeline Description

For each genome, the script performs the following steps:

1. Reads CDS nucleotide sequences using Biostrings.

2. Identifies highly expressed genes using header pattern matching:

ribosomal protein  
rps  
rpl  
elongation factor  

If no highly expressed genes are detected, prediction proceeds using default behavior.

3. Calls growth prediction using:

gRodon::predictGrowth(genes, highly_expressed)

4. Writes per-genome output file:

<OUTROOT>/<SAMPLE>/<SAMPLE>_gRodon.csv

5. Appends results to summary file using file locking to ensure concurrency safety:

<OUTROOT>/summary_grodon2.csv

---

## Output Files

Per-genome output:

<OUTROOT>/<SAMPLE>/<SAMPLE>_gRodon.csv

Combined summary output:

<OUTROOT>/summary_grodon2.csv

---

## Output Fields

Typical output variables include:

d  
Predicted doubling time (hours)

mu_max  
Maximum growth rate (calculated as ln(2)/d when available)

CUBHE  
Codon usage bias of highly expressed genes

ConsistencyHE  
Consistency of codon usage among highly expressed genes

GC  
GC content

nHE  
Number of highly expressed genes detected

---

## Reproducibility

Growth rate predictions were performed using gRodon v2.5.1 with CDS nucleotide sequences.

Execution parameters and workflow are fully defined in:

grodon2_batch.v2.sbatch