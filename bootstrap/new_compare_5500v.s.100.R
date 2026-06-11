#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# =========================
# paths
# =========================

boot_dir <- "/root/mcoa_project/mcoa_bootstrap_dt_balanced"

full_dir <- "/root/mcoa_project/check"

out_dir <- "/root/mcoa_project/bootstrap_summary_dt_balanced"

dir.create(
  out_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# =========================
# bootstrap runs
# =========================

runs <- list.dirs(
  boot_dir,
  recursive = FALSE,
  full.names = TRUE
)

cat("Found", length(runs), "runs\n")

# =========================
# collect all axis results
# =========================

all_list <- list()

for (r in runs) {

  run_name <- basename(r)

  for (ax in 1:2) {

    f <- file.path(
      r,
      paste0("axis", ax, ".tsv")
    )

    if (!file.exists(f)) {
      cat("Missing:", f, "\n")
      next
    }

    dt <- fread(f)

    if (nrow(dt) == 0)
      next

    if (!("var" %in% colnames(dt))) {
      cat("No var column:", f, "\n")
      next
    }

    dt[, run := run_name]
    dt[, axis := ax]
    dt[, rank := seq_len(.N)]

    all_list[[length(all_list) + 1]] <-
      dt[, .(var, axis, run, rank)]
  }
}

if (length(all_list) == 0) {
  stop("No axis tables found.")
}

boot_all <- rbindlist(all_list)

cat(
  "Loaded",
  nrow(boot_all),
  "records\n"
)

# =========================
# frequency + rank
# =========================

freq_table <- boot_all[
  ,
  .(
    freq = .N / length(runs),
    mean_rank = mean(rank)
  ),
  by = .(var, axis)
]

# =========================
# original full dataset
# =========================

get_full <- function(ax){

  f <- file.path(
    full_dir,
    paste0("axis", ax, ".tsv")
  )

  if (!file.exists(f)) {
    stop(
      paste(
        "Cannot find:",
        f
      )
    )
  }

  dt <- fread(f)

  if (!("var" %in% colnames(dt))) {
    stop(
      paste(
        "No var column in:",
        f
      )
    )
  }

  unique(dt$var)
}

full_axis1 <- get_full(1)
full_axis2 <- get_full(2)

# =========================
# in original dataset?
# =========================

freq_table[
  ,
  in_full := FALSE
]

freq_table[
  axis == 1 &
    var %in% full_axis1,
  in_full := TRUE
]

freq_table[
  axis == 2 &
    var %in% full_axis2,
  in_full := TRUE
]

# =========================
# Table 2
# =========================

table2 <- copy(freq_table)

setnames(
  table2,
  old = c(
    "var",
    "axis",
    "freq",
    "mean_rank",
    "in_full"
  ),
  new = c(
    "Variable",
    "Axis",
    "Frequency_across_50_runs",
    "Mean_rank_among_top_contributors",
    "Whether_in_original_dataset"
  )
)

setorder(
  table2,
  Axis,
  -Frequency_across_50_runs,
  Mean_rank_among_top_contributors
)

fwrite(
  table2,
  file.path(
    out_dir,
    "Table2_variable_stability.tsv"
  ),
  sep = "\t"
)

# =========================
# robust traits
# =========================

robust <- freq_table[
  freq >= 0.6
]

fwrite(
  robust,
  file.path(
    out_dir,
    "robust_traits.tsv"
  ),
  sep = "\t"
)

# =========================
# classification
# =========================

freq_table[
  ,
  category := "noise"
]

freq_table[
  freq >= 0.6 &
    in_full,
  category := "core"
]

freq_table[
  freq >= 0.6 &
    !in_full,
  category := "masked"
]

freq_table[
  freq < 0.6 &
    in_full,
  category := "full_only"
]

fwrite(
  freq_table,
  file.path(
    out_dir,
    "trait_classification.tsv"
  ),
  sep = "\t"
)

# =========================
# overlap summary
# =========================

summary_list <- list()

for (ax in 1:2) {

  full_vars <-
    if (ax == 1)
      full_axis1
    else
      full_axis2

  boot_vars <-
    unique(
      robust[
        axis == ax
      ]$var
    )

  inter <- intersect(
    full_vars,
    boot_vars
  )

  uni <- union(
    full_vars,
    boot_vars
  )

  summary_list[[ax]] <- data.table(

    axis = ax,

    n_full = length(full_vars),

    n_boot = length(boot_vars),

    overlap = length(inter),

    jaccard =
      length(inter) /
      length(uni)

  )
}

overlap_summary <-
  rbindlist(summary_list)

fwrite(
  overlap_summary,
  file.path(
    out_dir,
    "overlap_summary.tsv"
  ),
  sep = "\t"
)

cat(
  "\nFinished.\n",
  "Output directory:\n",
  out_dir,
  "\n"
)
