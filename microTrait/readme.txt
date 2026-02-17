# microTrait batch pipeline (hmmscan + microtrait mapping)

This pipeline runs microTrait functional trait inference in batch mode using:
1) HMMER `hmmscan` against the microTrait HMM database (and optionally a dbCAN HMM database),
2) normalization of `--domtblout` into a 23-column TSV,
3) trait mapping using the R package `microtrait`,
4) per-sample outputs and trait count tables at two granularities (g1 and g2).

All analyses were performed on the Rice University NOTSX high-performance computing cluster (Slurm).

---

## 1. Requirements

### Software
- HMMER (hmmscan)
- R (via the same environment)
- R packages:
  - `data.table`
  - `microtrait`

### Environment layout (as used in the script)
The Slurm script expects an environment directory containing both hmmscan and Rscript:

ENV=/scratch/$USER/envs/microtrait_env

The script uses:
- $ENV/bin/hmmscan
- $ENV/bin/Rscript

---

## 2. Databases

The pipeline requires the microTrait HMM database:

HMM_DIR=/scratch/$USER/microtrait-db

Expected file:
- $HMM_DIR/hmm/microtrait-hmmdb/microtrait.hmmdb

If the microTrait HMM file contains trusted cutoffs (lines starting with `TC`),
the script will automatically enable `hmmscan --cut_tc`.
Otherwise, it runs without `--cut_tc`.

### Optional dbCAN HMM database
If present, dbCAN hmmscan will also be executed. The script searches for one of:

- $HMM_DIR/hmm/dbcan/dbcan.noGA.hmmdb
- $HMM_DIR/hmm/dbcan/dbcan.select.v8.hmmdb

If neither exists, dbCAN is skipped and an empty dbCAN table is used as a fallback
for trait mapping.

---

## 3. Input: manifest TSV (two columns)

A manifest file is required. It must be a **two-column TSV**:

<sample_id>    /absolute/path/to/protein.faa

Example:
GCF_000005845.2    /scratch/cw106/eggnog_input/ncbi_dataset/data/GCF_000005845.2/protein.faa
GCF_000195955.2    /scratch/cw106/eggnog_input/ncbi_dataset/data/GCF_000195955.2/protein.faa

Notes:
- The first column (`sample_id`) is used as the output folder name.
- The second column must be a readable, non-empty `.faa` file.

IMPORTANT:
The Slurm array index is used as the manifest line number.
Therefore, the manifest should NOT include a header line unless you intentionally
account for it in the array range.

---

## 4. Running the pipeline

Slurm script:
- run_microtrait_batch.v3i.sbatch

Example submission (6000 samples, concurrency 30):

sbatch \
  --export=ALL,ENV=/scratch/cw106/envs/microtrait_env,HMM_DIR=/scratch/cw106/microtrait-db,MANIFEST=/scratch/cw106/microtrait_manifest.gcf.tsv,OUTROOT=/scratch/cw106/microtrait_batch_gcf \
  --array=1-6000%30 \
  /projects/alvarez/scripts/run_microtrait_batch.v3i.sbatch

---

## 5. Parameters and thresholds

These parameters can be overridden via `--export`:

ENV        : environment path (default: /scratch/$USER/envs/microtrait_env)
HMM_DIR    : database root (default: /scratch/$USER/microtrait-db)
MANIFEST   : manifest TSV (default: /scratch/$USER/microtrait_manifest.gcf.tsv)
OUTROOT    : output root (default: /scratch/$USER/microtrait_batch_gcf)

Filtering thresholds (defaults from the script):
- COV_GENE = 50
- COV_HMM  = 55
- IE_OR    = 1e-5
- PROT_G   = 45
- PROT_H   = 50

Filtering logic (implemented in R):
- For most hits: keep if (union gene coverage ≥ COV_GENE AND union HMM coverage ≥ COV_HMM) OR (i_evalue ≤ IE_OR)
- For hits matching protease-like HMM names (regex: `^(pep|lon|clp|hsl|ftsH|degP|prt|pepT)`):
  use PROT_G / PROT_H thresholds instead of COV_GENE / COV_HMM (still OR by i_evalue)

Coverage is computed as **union coverage** across all aligned segments per (gene, hmm) pair.

---

## 6. What the script produces

For each sample, outputs are written to:

$OUTROOT/<sample_id>/

Main files:
- microtrait.domtblout
- microtrait.hmmscan.log

Optional dbCAN files (if dbCAN DB is found):
- dbcan.domtblout
- dbcan.hmmscan.log

Sanitized domtblout (always generated):
- sanitized/microtrait.domtblout.tsv
- sanitized/dbcan.domtblout.tsv (empty file if dbCAN is skipped)

Trait outputs (core products):
- traits/trait_counts_g1.tsv
- traits/trait_counts_g2.tsv
- traits/genes_detected.tsv (if produced by microtrait mapping)

The main downstream inputs are:
- traits/trait_counts_g1.tsv
- traits/trait_counts_g2.tsv

---

## 7. Notes on reproducibility

- The pipeline uses HMMER hmmscan output normalized into a 23-column TSV format.
- Trait mapping is performed using the internal microtrait mapping function:
  `microtrait:::map.traits.fromdomtblout()`.
- dbCAN is optional; if absent, the pipeline still runs using an empty dbCAN table.
- All execution details (paths, thresholds, outputs) are defined in:
  run_microtrait_batch.v3i.sbatch
