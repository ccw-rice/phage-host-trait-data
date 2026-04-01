Overview

This repository provides an end-to-end pipeline for constructing the host life-history strategy space using multitable co-inertia analysis (MCOA) and identifying key variables associated with major axes.

Unlike earlier modular workflows, this pipeline integrates:

Data preprocessing
Adaptive transformation (ln(x+1))
Quality control (NA filtering)
Multitable co-inertia analysis (MCOA)
Variable–axis association analysis
Feature selection

All steps are performed within a single script, ensuring reproducibility and consistency.

📂 Input Data

All input files must be species-level aggregated feature matrices in TSV format.

Required files (placed in data_raw/):
Ko.merged_mean.tsv
Cog.merged_mean.tsv
Cazy.merged_mean.tsv
Traits.merged_mean.tsv
Format requirements:
Column 1: species ID
Column 2+: numeric feature values
All files must share consistent species IDs
Notes:
Non-feature files (e.g., genome_and_phylum) are automatically ignored
All values are coerced to numeric during preprocessing
⚙️ Pipeline Workflow

The entire workflow is implemented in:

all_in_one_mcoa_pipeline.R
Step 1 — Data loading and type conversion
All input tables are read into memory
All feature columns are converted to numeric
Non-numeric values are coerced (invalid values become NA)
Step 2 — Adaptive ln(x+1) transformation

Log transformation is selectively applied based on feature distribution.

A feature is transformed if:

All values are non-negative
AND one of the following is true:
Large dynamic range (max > 20)
Strong right-skewness (q90 / max < 0.1 and max > 1)

This approach targets:

Count-based features (e.g., KO, COG, CAZy)
Highly skewed abundance traits

and avoids transforming:

Bounded variables (e.g., GC content)
Approximately symmetric variables
Outputs:

For each input file:

Transformed matrix → ln1p/
List of transformed columns → *.ln1p_cols.txt
Step 3 — Data cleaning
Columns with all NA are removed
Rows containing any NA are removed

This ensures compatibility with PCA and MCOA.

Step 4 — Standardization

All matrices are standardized using:

vegan::decostand(method = "standardize")
Step 5 — Multitable co-inertia analysis (MCOA)

For each dataset:

PCA is performed using ade4::dudi.pca
Tables are combined into a multitable structure
MCOA is performed using ade4::mcoa
Step 6 — Strategy space construction
Species coordinates are extracted (SynVar)
Axis directions are flipped for consistency:
SynVar1 = -SynVar1
SynVar2 = -SynVar2
Output:
samples.tsv

Contains:

Species ID
MCOA coordinates (SynVar1, SynVar2, …)
Step 7 — Variable contribution (Tco)
Variable contributions to axes are extracted from MCOA
Axis signs are adjusted to match SynVar
Output:
variables_contribution_Tco.tsv
Step 8 — Variable–axis association analysis

For each variable:

Linear models are fitted:
feature ~ SynVar1
feature ~ SynVar2

Metrics computed:

Regression coefficient
p-value
R²
Output:
variables_vs_axes_lm.tsv
Step 9 — Feature selection

For each axis:

Selection criteria:

p-value < 0.001
R² > 0.2
Ranked by |SV| (MCOA contribution)

Top 50 variables are retained.

Outputs:
axis1.tsv
axis2.tsv
📊 Output Structure

All outputs are stored in:

check/
├── ln1p/                         # transformed matrices + logs
│   ├── *.tsv
│   ├── *.ln1p_cols.txt
│
├── samples.tsv                  # species coordinates
├── variables_contribution_Tco.tsv
├── variables_vs_axes_lm.tsv
├── axis1.tsv
├── axis2.tsv
📦 Software Requirements

Required R packages:

data.table
ade4
vegan
Installation example:
micromamba activate r_ade4
micromamba install r-data.table r-ade4 r-vegan
▶️ How to Run
Rscript all_in_one_mcoa_pipeline.R
⚠️ Important Notes
The pipeline is deterministic (no random sampling is performed)
set.seed(123) is included but not required for current execution
Adaptive transformation introduces data-dependent preprocessing, which should be documented when used in publications
Sample size may decrease after NA filtering
🧠 Key Features of This Pipeline
Fully automated end-to-end workflow
Adaptive preprocessing based on data distribution
Built-in quality control
Reproducible and transparent transformation logging
Direct linkage between MCOA axes and biological features
🧭 Recommended Usage
Suitable for large-scale comparative genomics analyses
Recommended for exploratory analysis and robust feature selection
For publication, clearly describe the transformation criteria used