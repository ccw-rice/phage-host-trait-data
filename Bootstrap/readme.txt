MCOA Bootstrap Pipeline for Host Strategy Space Analysis
🧠 Overview

This repository contains a reproducible pipeline to construct and evaluate host life-history strategy space using multitable co-inertia analysis (MCOA), with a focus on:

Robustness of MCOA structure under subsampling
Association between host traits and MCOA axes
Spatial patterns of prophage burden and virulent phage hosts

The pipeline integrates multiple genome-derived feature tables (KO, COG, CAZy, Traits), applies dynamic transformation, and performs bootstrap resampling across phylogenetic groups (phylum level) to assess stability.

⚙️ Pipeline Structure

The workflow consists of three main components:

1. Bootstrap MCOA Construction

Core script:
👉

Key steps:
Phylogenetically stratified subsampling
Up to 100 genomes per phylum
Reduces taxonomic bias
Dynamic data transformation (ln(x+1))
Applied only when needed
Based on distribution criteria:
non-negative
high dynamic range or strong skew
Matrix cleaning
Remove NA-heavy rows/columns
MCOA construction
PCA per dataset (KO, COG, CAZy, Traits)
Integrated using ade4::mcoa
Outputs per run (×50 runs):
samples.tsv → genome coordinates (MCOA space)
variables_contribution_Tco.tsv → variable contributions
variables_vs_axes_lm.tsv → linear models (trait vs axis)
axis1.tsv, axis2.tsv → top contributing variables
2. Stability Analysis of MCOA Space

Core script:
👉

What it evaluates:
Axis consistency
Correlation of SynVar1 / SynVar2 between runs
Global structure similarity
Procrustes analysis
Genome-level displacement
Euclidean shift between bootstrap and full dataset
Outputs:
sample_stability_summary.tsv
cor_x, cor_y
Procrustes R
mean / median shift
sample_shift_all.tsv
per-genome displacement across runs
3. Trait Projection in MCOA Space
(A) Prophage trajectory

Script:
👉

Calculates centroid positions by prophage count group
Aligns all bootstrap runs to a reference coordinate system
Produces:

Outputs:

centroids_all.tsv
centroid_mean.tsv
trajectory_mean.png
(B) Virulent phage distribution

Script:
👉

Focuses on genomes with virulent phage presence
Computes centroid per bootstrap run
Produces spatial distribution (“cloud”)

Outputs:

virulent_centroids.tsv
virulent_centroid.png
4. Optional: Subsampling Only

Script:
👉

Performs phylum-balanced subsampling
Outputs subset tables for all input matrices
📂 Input Data

All input matrices must be species-level aggregated TSV tables.

Required files:

Ko.merged_mean.tsv
Cog.merged_mean.tsv
Cazy.merged_mean.tsv
Traits.merged_mean.tsv
Prophage.merged_mean.tsv
Virulent.merged_mean.tsv
genome_and_phylum.csv

Each table format:

ID    feature1    feature2    ...
📁 Directory Structure
mcoa_project/
│
├── data_raw/                 # input matrices
├── data_raw_subset/          # subsampled data (optional)
├── mcoa_bootstrap/           # bootstrap runs
│   ├── run_01/
│   ├── run_02/
│   └── ...
│
├── check/
│   └── samples.tsv           # full dataset reference
│
├── bootstrap_sample_comparison/
├── prophage_trend_aligned/
└── virulent_trend_aligned/
🔁 Workflow Summary
Step 1: Subsample genomes by phylum
Step 2: Apply dynamic ln(x+1) transformation
Step 3: Construct MCOA (50 bootstrap runs)
Step 4: Align coordinates to full dataset
Step 5: Evaluate stability (correlation + Procrustes)
Step 6: Map biological signals:
        - prophage burden trajectory
        - virulent host distribution
🔬 Key Design Features
✅ Dynamic ln(x+1) (important difference vs old pipeline)
Not all features are transformed
Avoids over-normalization of already well-behaved variables
Improves interpretability of trait contributions
✅ Phylum-balanced bootstrap
Controls for taxonomic overrepresentation
Ensures robustness is not driven by dominant clades (e.g., Pseudomonadota)
✅ Coordinate alignment (critical)

Across all scripts:

if (cor(x_boot, x_full) < 0) flip
Ensures axis direction consistency
Prevents artificial inversion across runs
✅ Multi-level validation
Axis correlation (rank consistency)
Procrustes (global geometry)
Centroid shifts (biological interpretation)
📊 Outputs Interpretation
High cor_x / cor_y (~0.95+)
→ stable axis structure
Low mean shift
→ consistent genome positioning
Centroid trajectory
→ biological gradients (e.g., prophage accumulation)
Virulent cloud
→ spatial enrichment pattern in strategy space
🚀 How to Run

Example:

Rscript run_bootstrap_pipeline_dynamic_ln1p.R
Rscript compare_coordinates.R
Rscript prophage_coordinates.R
Rscript virulent_phage_coordinates.R
📌 Notes
All analyses were performed on HPC (Linux environment)
Requires R packages:
data.table
ade4
vegan
ggplot2
