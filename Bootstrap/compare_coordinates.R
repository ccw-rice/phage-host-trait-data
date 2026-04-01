#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(vegan)
})

# =========================
# 路径
# =========================
boot_dir <- "/root/mcoa_project/mcoa_bootstrap"
full_file <- "/root/mcoa_project/check/samples.tsv"
out_dir <- "/root/mcoa_project/bootstrap_sample_comparison"

dir.create(out_dir, showWarnings=FALSE)

# =========================
# 读取 full data
# =========================
full <- fread(full_file)

full <- full[, .(id, x = SynVar1, y = SynVar2)]

results <- list()
shift_list <- list()

runs <- list.dirs(boot_dir, recursive=FALSE)

for (r in runs){

  run_name <- basename(r)
  f <- file.path(r, "samples.tsv")

  if (!file.exists(f)) next

  boot <- fread(f)
  boot <- boot[, .(id, x = SynVar1, y = SynVar2)]

  # =========================
  # 对齐样本
  # =========================
  df <- merge(full, boot, by="id", suffixes=c("_full","_boot"))

  if (nrow(df) < 50) next  # 防止异常

  # =========================
  # ⭐方向对齐（关键）
  # =========================
  if (cor(df$x_full, df$x_boot) < 0){
    df$x_boot <- -df$x_boot
  }

  if (cor(df$y_full, df$y_boot) < 0){
    df$y_boot <- -df$y_boot
  }

  # =========================
  # 相关性
  # =========================
  cor_x <- cor(df$x_full, df$x_boot)
  cor_y <- cor(df$y_full, df$y_boot)

  # =========================
  # Procrustes
  # =========================
  proc <- protest(
    as.matrix(df[,.(x_full,y_full)]),
    as.matrix(df[,.(x_boot,y_boot)]),
    permutations = 0
  )

  # =========================
  # 位移
  # =========================
  df[, shift := sqrt((x_full - x_boot)^2 + (y_full - y_boot)^2)]
  df[, run := run_name]

  shift_list[[run_name]] <- df[,.(id, run, shift)]

  results[[run_name]] <- data.table(
    run = run_name,
    cor_x = cor_x,
    cor_y = cor_y,
    procrustes_R = proc$t0,
    mean_shift = mean(df$shift),
    median_shift = median(df$shift)
  )
}

# =========================
# 汇总
# =========================
summary_table <- rbindlist(results)
shift_table   <- rbindlist(shift_list)

fwrite(summary_table,
       file.path(out_dir, "sample_stability_summary.tsv"),
       sep="\t")

fwrite(shift_table,
       file.path(out_dir, "sample_shift_all.tsv"),
       sep="\t")

cat("\n🔥 DONE: sample stability analysis complete\n")
