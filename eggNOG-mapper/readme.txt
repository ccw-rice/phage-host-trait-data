# eggNOG-mapper Functional Annotation Pipeline

This pipeline performs genome-wide functional annotation using eggNOG-mapper (emapper.py) with DIAMOND search.

This analysis was used to assign functional traits for downstream host life-history strategy analysis.

------------------------------------------------------------
Software requirements
------------------------------------------------------------

eggNOG-mapper v2.x
DIAMOND v2.x
Micromamba / Conda

eggNOG database files:

eggnog_proteins.dmnd
eggnog.db
eggnog.taxa.db
eggnog.taxa.db.traverse.pkl

------------------------------------------------------------
Input data structure
------------------------------------------------------------

Each genome must contain:

protein.faa

Directory structure:

DATA_ROOT/
    GCF_xxxxxxxx.x/
        protein.faa

Example:

/scratch/cw106/eggnog_input/ncbi_dataset/data/GCF_000005845.2/protein.faa

------------------------------------------------------------
Pipeline overview
------------------------------------------------------------

Step 1. Generate genome list

Example:

find DATA_ROOT -maxdepth 1 -type d -name 'GCF_*' \
> genomes.list

------------------------------------------------------------

Step 2. Submit batch job

Example:

sbatch emap_groups_100x10.sbatch

------------------------------------------------------------

Step 3. Protein merging and annotation

For each array task:

10 genomes are merged into a single protein FASTA file

Each protein ID is prefixed with genome ID:

Example:

GCF_xxxxx|protein_id

Annotation performed using:

emapper.py \
-i merged.faa \
-m diamond \
--cpu 16 \
--sensmode fast \
--dmnd_iterate yes \
--dbmem \
--target_taxa 2

------------------------------------------------------------

Step 4. Output splitting

Combined annotation results are split back into individual genome files:

Output structure:

OUTPUT_ROOT/

    GCF_xxxxxxxx.x.emapper.annotations

------------------------------------------------------------

Output file format
------------------------------------------------------------

Each output file contains eggNOG functional annotations:

COG category
Functional description
KEGG ortholog
Gene ontology
etc.

------------------------------------------------------------

Reproducibility
------------------------------------------------------------

Annotations were performed using eggNOG-mapper with DIAMOND search mode and bacterial taxonomic scope (target_taxa=2).

Full parameters are provided in emap_groups_100x10.sbatch.
