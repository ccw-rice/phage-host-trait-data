# Virulent Phage identification pipeline (HMM-based)

This directory contains scripts and marker lists used to infer phage lifestyle (virulent) based on genome-wide HMM annotation.

This pipeline was used in:

Wu et al. Host life-history strategies structure phage lifestyle associations across bacterial species.

---

# Overview

Phage genomes were annotated using Prodigal and HMMER against a merged PHROG + VOG HMM database. Lifestyle was inferred based on the presence of hallmark genes including integrase, terminase large subunit (terL), and lysis-related genes.

All analyses were performed on the Rice University NOTSX high-performance computing cluster.

---

# Input data

Phage genome FASTA files organized as:

DATA_ROOT/<Phage_source>/<Phage_source>/<Phage_ID>.fasta(.gz)


Master index file:

phage_id_list_2col.ALL15k.tsv


Format:

Phage_source Phage_ID


---

# Pipeline steps

## Step 1. Gene prediction and HMM annotation

Script:

phage_prodigal_hmm_array.sbatch


Steps performed:

- Gene prediction using Prodigal v2.6.3
- Functional annotation using HMMER v3.4
- Database: merged PHROG + VOG HMM database

Outputs:

<Phage_ID>.faa
<Phage_ID>.hmm.tbl


---

## Step 2. Lifestyle classification

Scripts:

call_lifestyle_from_tbl.py
call_lifestyle_from_tbl.sbatch


This script:

- Reads HMM annotation results (.hmm.tbl)
- Counts hits to marker genes
- Calculates genome length
- Assigns lifestyle classification

---

# Marker lists

The following marker lists define hallmark phage genes:

integrase_unified.list
terl.list
lysis_genes.list


These lists were constructed from PHROG and VOG databases.

---

# Lifestyle classification rules

Lifestyle is assigned using the following logic:

if Has_integrase == 1:
Lifestyle_call = "temperate"

elif Has_terL == 1 and Has_lysis == 1:
Lifestyle_call = "strong_virulent"

elif Has_terL == 1:
Lifestyle_call = "virulent"

else:
Lifestyle_call = "uncertain"


Definitions:

- integrase: marker of lysogeny capability
- terL: marker of dsDNA tailed phages
- lysis genes: holin, endolysin, spanin, and related genes

---

# Output

Final output table:

phage_lifestyle_calls_all15k.tsv


Columns:

Phage_source
Phage_ID
Genome_length_bp
n_unique_hmm_hits
n_integrase_hits
n_terL_hits
n_lysis_hits
Has_integrase
Has_terL
Has_lysis
Lifestyle_call


---

# Software used

- Prodigal v2.6.3
- HMMER v3.4
- Python 3
- Slurm workload manager

---

# Computational environment

All analyses were performed on:

Rice University NOTSX high-performance computing cluster

---

# Reproducibility

All scripts and marker lists required to reproduce lifestyle classification are provided in this directory.

External databases required:

- PHROG database
- VOG database
- merged HMM database

---

# Contact

Chuncheng Wu  
Rice University  