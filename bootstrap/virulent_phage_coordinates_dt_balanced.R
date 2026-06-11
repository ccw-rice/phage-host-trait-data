#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# =========================
# Paths (only changed here)
# =========================
boot_dir  <- "/root/mcoa_project/mcoa_bootstrap_dt_balanced"
vir_file  <- "/root/mcoa_project/data_raw/Virulent.merged_mean.tsv"
full_file <- "/root/mcoa_project/check/samples.tsv"
out_dir   <- "/root/mcoa_project/virulent_trend_aligned_dt"

dir.create(out_dir, showWarnings=FALSE)

# =========================
# Utility functions
# =========================
clean_id <- function(x){
  x <- trimws(x)
  x <- gsub("\\s+", " ", x)
  return(x)
}

# =========================
# Full dataset (reference coordinate system)
# =========================
full <- fread(full_file)
full <- full[, .(id, x_full = SynVar1, y_full = SynVar2)]
full[, id := clean_id(id)]

# =========================
# Virulent phage data (unchanged)
# =========================
vir <- fread(vir_file)

rn <- vir[[1]]
vir[[1]] <- NULL

vir[] <- lapply(vir, function(x)
  suppressWarnings(as.numeric(as.character(x)))
)

# Core step: rowSums (already correct)
vir$count <- rowSums(vir, na.rm=TRUE)

vir <- data.table(
  id = clean_id(rn),
  virulent = vir$count
)

# Retain only genomes with virulent > 0
vir_pos <- vir[virulent > 0]

cat("Total virulent genomes:", nrow(vir_pos), "\n")

# =========================
# Bootstrap analysis
# =========================
runs <- list.dirs(boot_dir, recursive=FALSE)

res_list <- list()

for (r in runs){

  run_name <- basename(r)
  f <- file.path(r, "samples.tsv")

  if (!file.exists(f)) next

  cat("Processing:", run_name, "\n")

  dt <- fread(f)
  dt <- dt[, .(id, x = SynVar1, y = SynVar2)]
  dt[, id := clean_id(id)]

  # =========================
  # Align to full dataset
  # =========================
  df_align <- merge(dt, full, by="id")

  if (nrow(df_align) < 50){
    warning("Too few matched genomes: ", run_name)
    next
  }

  if (cor(df_align$x, df_align$x_full) < 0){
    dt$x <- -dt$x
  }

  if (cor(df_align$y, df_align$y_full) < 0){
    dt$y <- -dt$y
  }

  # =========================
  # Virulent > 0 genomes
  # =========================
  df <- merge(dt, vir_pos, by="id")

  if (nrow(df) < 10){
    warning("Too few virulent genomes in ", run_name)
    next
  }

  # =========================
  # Centroid
  # =========================
  res_list[[run_name]] <- data.table(
    run = run_name,
    x_mean = mean(df$x),
    y_mean = mean(df$y),
    n = nrow(df)
  )
}

res_all <- rbindlist(res_list)

# =========================
# Save results
# =========================
fwrite(res_all,
       file.path(out_dir, "virulent_centroids.tsv"),
       sep="\t")

# =========================
# Plot
# =========================
p <- ggplot(res_all, aes(x=x_mean, y=y_mean)) +

  geom_point(alpha=0.6, color="gray") +

  coord_fixed(xlim = c(-3, 3), ylim = c(-3, 3)) +

  geom_hline(yintercept = 0, linetype="dashed", color="gray") +
  geom_vline(xintercept = 0, linetype="dashed", color="gray") +

  labs(
    x="MCOA1",
    y="MCOA2",
    title="Virulent genomes (dt-balanced)"
  ) +

  theme_bw(base_size=14)

ggsave(
  file.path(out_dir, "virulent_centroid.png"),
  plot = p,
  width = 6,
  height = 5,
  dpi = 300
)

cat("\nDONE: virulent bootstrap cloud (dt-balanced)\n")

