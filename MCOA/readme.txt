Overview

This repository contains the scripts used to construct the host life-history strategy space using multitable co-inertia analysis (MCOA), identify key contributing variables, and evaluate their associations with prophage burden.

All analyses were performed on the Rice University NOTSX high-performance computing cluster using GTDB species-level aggregated genomic trait matrices.

This pipeline consists of three main steps:

Construction of host strategy space using MCOA

Identification of key variables contributing to MCOA axes

Correlation analysis between key variables and prophage burden

Input Data

All input matrices must be species-level aggregated feature tables in TSV format.

Required files:

Ko.merged_mean.tsv or Ko.merged_mean.ln1p.tsv
Cog.merged_mean.tsv or Cog.merged_mean.ln1p.tsv
Cazy.merged_mean.tsv or Cazy.merged_mean.ln1p.tsv
Traits.merged_mean.tsv or Traits.merged_mean.ln1p.tsv
Prophage.merged_mean.ln1p.tsv


Each file must have:

Column 1: species ID
Column 2+: feature values


Rows must correspond to the same GTDB species identifiers across all tables.

Software Requirements

The following R packages are required:

ade4
vegan
data.table
dplyr
ggplot2
ggpubr
viridis
tibble


Install example:

micromamba activate r_ade4
micromamba install r-ade4 r-vegan r-data.table r-dplyr r-ggplot2 r-ggpubr r-viridis r-tibble

Step 1 — Construct Host Life-History Strategy Space

Script:

mcoa_author_style.sbatch


Run:

sbatch mcoa_author_style.sbatch


This script performs:

• Standardization of feature matrices
• Principal component analysis (PCA) on each table
• Multitable co-inertia analysis (MCOA) integration
• Projection of species into shared strategy space

Outputs:

mcoa_author_style/

samples_SynVar123.tsv
variables_contribution_Tco.tsv
variables_vs_axes_lm.tsv
axis_variation_share.tsv
MCOA_samples_SynVar1_vs_2.pdf


Key output:

samples_SynVar123.tsv


Contains species coordinates in host strategy space (MCOA1, MCOA2, MCOA3).

Step 2 — Identify Key Variables Contributing to MCOA Axes

Script:

filter_best_vars.R


Run inside MCOA output directory:

Rscript filter_best_vars.R


Selection criteria:

p < 0.001
R² > 0.2
Top 50 variables ranked by absolute contribution (|SV|)


Outputs:

best_vars_filtered/

best_vars_axis1_filtered.tsv
best_vars_axis2_filtered.tsv
best_vars_axis3_filtered.tsv
filter_summary.tsv


These represent variables strongly associated with host strategy axes.