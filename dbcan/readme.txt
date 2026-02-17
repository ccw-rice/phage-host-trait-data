# dbCAN CAZyme Annotation Pipeline

This pipeline performs genome-wide carbohydrate-active enzyme (CAZyme) annotation using run_dbcan (v2.0.11) with HMMER and DIAMOND, applied to bacterial protein sequences.

This analysis was used to quantify host resource-acquisition traits for downstream life-history strategy analysis.

------------------------------------------------------------
Software requirements
------------------------------------------------------------

run_dbcan v2.0.11
HMMER v3.3+
DIAMOND v2.0+
Micromamba / Conda environment

dbCAN database v14:

dbCAN-HMMdb-V14.txt
CAZy.dmnd

------------------------------------------------------------
Input data structure
------------------------------------------------------------

Each genome must contain a protein FASTA file:

DATA_ROOT/
  GCF_xxxxxxxx.x/
      protein.faa

Example:

/scratch/cw106/eggnog_input/ncbi_dataset/data/GCF_000005845.2/protein.faa

A genome list file is required:

gcf_list.txt

Format:

GCF_000005845.2
GCF_000006765.1
...

------------------------------------------------------------
Pipeline overview
------------------------------------------------------------

Step 1. Generate genome list

Example:

find DATA_ROOT -name protein.faa \
| sed 's|DATA_ROOT/||; s|/protein.faa||' \
> gcf_list.txt

------------------------------------------------------------

Step 2. Submit recursive batch jobs

Example:

sbatch \
--export=ALL,LIST=gcf_list.txt,WORKER=run_dbcan_array_worker.sbatch,START=0,N=<TOTAL>,CHUNK=200,MAXPAR=30 \
submit_dbcan_recursive.sbatch

------------------------------------------------------------

Step 3. Worker job execution

Each worker performs:

run_dbcan.py protein.faa protein \
--tools hmmer diamond \
--hmm_eval 1e-15 \
--hmm_cov 0.35 \
--dia_eval 1e-102

------------------------------------------------------------

Output structure
------------------------------------------------------------

run_dbcan_out/

  GCF_xxxxxxxx.x/

      overview.tsv

      hmmer.out

      diamond.out

overview.tsv contains genome-level CAZyme annotations.

------------------------------------------------------------

Output interpretation
------------------------------------------------------------

overview.tsv provides:

CAZyme family assignments
Gene-level annotations
Functional trait information

These annotations were used for downstream functional trait and life-history analyses.

------------------------------------------------------------

Reproducibility
------------------------------------------------------------

All genome annotations were performed using run_dbcan v2.0.11 with dbCAN v14 database and default thresholds as specified in run_dbcan_array_worker.sbatch.

------------------------------------------------------------
