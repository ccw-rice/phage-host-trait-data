# CRISPR–Cas System Identification Using CRISPRCasTyper

## Overview

CRISPR–Cas systems were identified from bacterial genome assemblies using CRISPRCasTyper. Batch analysis was performed on a high-performance computing (HPC) cluster using SLURM array jobs. A custom summary script was used to extract intact CRISPR–Cas system information and generate genome-level summary tables.

This pipeline identifies:

* Presence or absence of intact CRISPR–Cas systems
* Number of intact operons
* CRISPR–Cas system types and subtypes
* Number of Cas genes
* Number of spacers associated with intact systems
* Number of orphan CRISPR arrays

---

## Software and Environment

CRISPRCasTyper was installed in a dedicated Conda environment:

```
/scratch/cw106/cctyper_env_new
```

The CRISPRCasTyper database was installed at:

```
/scratch/cw106/cctyper_db/data
```

Environment activation:

```
module purge
module load Miniforge3/24.1.2-0

eval "$(conda shell.bash hook)"
conda activate /scratch/cw106/cctyper_env_new

export CCTYPER_DB=/scratch/cw106/cctyper_db/data
```

---

## Input Data

Input files consisted of bacterial genome assemblies in FASTA format:

```
*.fna
```

Example input:

```
/scratch/cw106/displace_prophage/crispr_defense_input_hostonly/GCF_000005845.2.fna
```

A manifest file was created listing all genome paths:

```
/scratch/cw106/displace_prophage/crispr_defense_manifest_hostonly.txt
```

Format:

```
/absolute/path/to/genome1.fna
/absolute/path/to/genome2.fna
...
```

Each line corresponds to one genome.

---

## Batch Execution

Batch processing was performed using SLURM array jobs with the script:

```
scripts/cctyper_hostonly_array.sbatch
```

Submission example:

```
sbatch --array=1-N%20 scripts/cctyper_hostonly_array.sbatch
```

This script automatically:

* Loads the Conda environment
* Reads genome paths from the manifest file
* Executes CRISPRCasTyper
* Writes output to per-genome directories

Output directory structure:

```
/scratch/cw106/displace_prophage/crisprcas_out_hostonly/

├── GCF_000005845.2_cctyper/
│   ├── cas_operons.tab
│   ├── CRISPR_Cas.tab
│   ├── crisprs_all.tab
│   ├── crisprs_orphan.tab
│   ├── crisprs_near_cas.tab
│   ├── hmmer.tab
│   ├── spacers/
│   └── plot.svg
```

---

## Summary Table Generation

Genome-level summary tables were generated using:

```
scripts/cctyper_intact_systems_v3.py
```

Execution example:

```
python scripts/cctyper_intact_systems_v3.py \
    /scratch/cw106/displace_prophage/crisprcas_out_hostonly \
    --pattern "*_cctyper" \
    --merge-near-as-orphan \
    > intact_summary_hostonly.tsv
```

---

## Summary Output

The output file contains one row per genome:

```
intact_summary_hostonly.tsv
```

Columns include:

| Column                      | Description                         |
| --------------------------- | ----------------------------------- |
| sample                      | Genome identifier                   |
| has_intact_system           | YES / NO                            |
| intact_systems              | Number of intact CRISPR–Cas systems |
| intact_cas_genes            | Total Cas genes in intact systems   |
| intact_spacers              | Number of spacers in intact systems |
| Type_I–Type_VI              | Counts per CRISPR–Cas type          |
| intact_subtypes             | Subtype list                        |
| intact_operons              | Operon identifiers                  |
| orphan_trusted_arrays_count | Number of orphan CRISPR arrays      |
| orphan_trusted_spacers      | Number of orphan spacers            |

---

## Criteria for Intact Systems

A CRISPR–Cas system was considered intact if it was listed in:

```
cas_operons.tab
```

Spacer counts were derived from:

```
spacers/*.fa
```

Orphan CRISPR arrays were identified from:

```
crisprs_orphan.tab
```

and optionally merged with nearby arrays using the `--merge-near-as-orphan` option.

---

## Scripts Provided

The following scripts are included in this repository:

```
scripts/cctyper_hostonly_array.sbatch
scripts/cctyper_intact_systems_v3.py
```

These scripts fully reproduce the CRISPR–Cas detection and summary pipeline.

---

## Reproducibility Notes

* All analyses were performed on complete bacterial genomes
* Default CRISPRCasTyper detection parameters were used
* Genome assemblies were obtained from NCBI RefSeq

---

## Citation

If using CRISPRCasTyper, please cite:

Russel J, et al. CRISPRCasTyper: automated identification and classification of CRISPR–Cas loci.
Nucleic Acids Research.

---

## Contact

For questions regarding this pipeline, please contact:

Chuncheng Wu
Rice University
Department of Civil and Environmental Engineering
