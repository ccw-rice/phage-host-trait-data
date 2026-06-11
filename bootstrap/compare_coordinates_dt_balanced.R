#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(vegan)
})

# ==================================================
# paths
# ==================================================

base_dir <- "/root/mcoa_project/mcoa_bootstrap_dt_balanced"

ref_file <- "/root/mcoa_project/check/samples.tsv"

out_dir <- "/root/mcoa_project/bootstrap_sample_comparison_dt_balanced"

dir.create(
  out_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ==================================================
# load reference
# ==================================================

ref <- fread(ref_file)

if(!("id" %in% colnames(ref))){
  stop("Reference samples.tsv has no id column")
}

setkey(ref, id)

runs <- list.dirs(
  base_dir,
  recursive = FALSE,
  full.names = TRUE
)

cat(
  "Found",
  length(runs),
  "bootstrap runs\n"
)

# ==================================================
# compare
# ==================================================

res <- list()

shift_list <- list()

for(r in runs){

  run_name <- basename(r)

  samp_file <- file.path(
    r,
    "samples.tsv"
  )

  if(!file.exists(samp_file)){
    cat("Missing:", samp_file, "\n")
    next
  }

  samp <- fread(samp_file)

  if(!("id" %in% colnames(samp))){
    cat("No id column:", samp_file, "\n")
    next
  }

  setkey(samp, id)

  common <- intersect(
    ref$id,
    samp$id
  )

  if(length(common) < 10){
    cat(
      run_name,
      "too few shared genomes\n"
    )
    next
  }

  x_ref <- ref[common]
  x_run <- samp[common]

  # ==========================================
  # align directions
  # ==========================================

  if(cor(
      x_ref$SynVar1,
      x_run$SynVar1
    ) < 0){

    x_run$SynVar1 <- -x_run$SynVar1
  }

  if(cor(
      x_ref$SynVar2,
      x_run$SynVar2
    ) < 0){

    x_run$SynVar2 <- -x_run$SynVar2
  }

  # ==========================================
  # correlations
  # ==========================================

  cor_x <- cor(
    x_ref$SynVar1,
    x_run$SynVar1,
    method = "spearman"
  )

  cor_y <- cor(
    x_ref$SynVar2,
    x_run$SynVar2,
    method = "spearman"
  )

  # ==========================================
  # Procrustes
  # ==========================================

  proc <- protest(
    as.matrix(
      x_ref[,.(SynVar1, SynVar2)]
    ),
    as.matrix(
      x_run[,.(SynVar1, SynVar2)]
    )
  )

  # ==========================================
  # shifts
  # ==========================================

  shift <- sqrt(

    (x_ref$SynVar1 -
       x_run$SynVar1)^2 +

    (x_ref$SynVar2 -
       x_run$SynVar2)^2

  )

  shift_list[[run_name]] <- data.table(

    run = run_name,

    id = common,

    shift = shift

  )

  # ==========================================
  # summary
  # ==========================================

  res[[run_name]] <- data.table(

    run = run_name,

    cor_x = cor_x,

    cor_y = cor_y,

    procrustes_r = proc$t0,

    mean_shift = mean(shift),

    median_shift = median(shift),

    max_shift = max(shift)

  )

}

# ==================================================
# combine
# ==================================================

summary_dt <- rbindlist(res)

shift_dt <- rbindlist(shift_list)

# ==================================================
# save
# ==================================================

fwrite(
  summary_dt,
  file.path(
    out_dir,
    "sample_stability_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  shift_dt,
  file.path(
    out_dir,
    "sample_shift_all.tsv"
  ),
  sep = "\t"
)

# ==================================================
# Table 3
# ==================================================

table3 <- summary_dt[, .(

  Bootstrap_run = run,

  MCOA1_correlation =
    round(cor_x, 2),

  MCOA2_correlation =
    round(cor_y, 2),

  Mean_displacement =
    round(mean_shift, 2)

)]

fwrite(
  table3,
  file.path(
    out_dir,
    "Table3_genome_position_stability.tsv"
  ),
  sep = "\t"
)

# ==================================================
# overall stats
# ==================================================

overall <- data.table(

  mean_cor_x =
    mean(summary_dt$cor_x),

  mean_cor_y =
    mean(summary_dt$cor_y),

  mean_shift =
    mean(summary_dt$mean_shift),

  median_shift =
    median(summary_dt$mean_shift)

)

fwrite(
  overall,
  file.path(
    out_dir,
    "overall_summary.tsv"
  ),
  sep = "\t"
)

cat(
  "\nDONE\n",
  "Output:\n",
  out_dir,
  "\n"
)
