#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# =========================
# 路径（只改这里）
# =========================
boot_dir  <- "/root/mcoa_project/mcoa_bootstrap_dt_balanced"
prop_file <- "/root/mcoa_project/data_raw/Prophage.merged_mean.tsv"
full_file <- "/root/mcoa_project/check/samples.tsv"
out_dir   <- "/root/mcoa_project/prophage_trend_aligned_dt"

dir.create(out_dir, showWarnings=FALSE)

# =========================
# 工具函数
# =========================
clean_id <- function(x){
  x <- trimws(x)
  x <- gsub("\\s+", " ", x)
  return(x)
}

# =========================
# full data（参考坐标系）
# =========================
full <- fread(full_file)
full <- full[, .(id, x_full=SynVar1, y_full=SynVar2)]
full[, id := clean_id(id)]

# =========================
# 🔥 prophage（完全按旧逻辑）
# =========================
prop <- fread(prop_file)

rn <- prop[[1]]
prop[[1]] <- NULL

# 全部转 numeric
prop[] <- lapply(prop, function(x)
  suppressWarnings(as.numeric(as.character(x)))
)

# ✔ 核心：rowSums
prop$count <- rowSums(prop, na.rm=TRUE)

prop <- data.table(
  id = clean_id(rn),
  prophage = round(prop$count)
)

# ✔ 分组（完全复制旧逻辑）
prop[, group := ifelse(prophage >= 6, 6, prophage)]

# =========================
# bootstrap
# =========================
runs <- list.dirs(boot_dir, recursive=FALSE)

all_centroids <- list()

for (r in runs){

  run_name <- basename(r)
  f <- file.path(r, "samples.tsv")
  if (!file.exists(f)) next

  dt <- fread(f)
  dt <- dt[, .(id, x=SynVar1, y=SynVar2)]
  dt[, id := clean_id(id)]

  # =========================
  # 对齐到 full（关键）
  # =========================
  df_align <- merge(dt, full, by="id")

  if (cor(df_align$x, df_align$x_full) < 0){
    dt$x <- -dt$x
  }

  if (cor(df_align$y, df_align$y_full) < 0){
    dt$y <- -dt$y
  }

  # =========================
  # 合并 prophage
  # =========================
  df <- merge(dt, prop, by="id")

  # =========================
  # 计算 centroid
  # =========================
  cent <- df[, .(
    x_mean = mean(x),
    y_mean = mean(y)
  ), by=group]

  cent[, run := run_name]

  all_centroids[[run_name]] <- cent
}

centroids_all <- rbindlist(all_centroids)

# =========================
# 平均轨迹
# =========================
centroid_mean <- centroids_all[, .(
  x_mean = mean(x_mean),
  y_mean = mean(y_mean),
  x_sd   = sd(x_mean),
  y_sd   = sd(y_mean)
), by=group]

# =========================
# 保存
# =========================
fwrite(centroids_all,
       file.path(out_dir, "centroids_all.tsv"), sep="\t")

fwrite(centroid_mean,
       file.path(out_dir, "centroid_mean.tsv"), sep="\t")

# =========================
# 画图（完全保留）
# =========================
p <- ggplot(centroid_mean, aes(x=x_mean, y=y_mean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=y_mean-y_sd, ymax=y_mean+y_sd), width=0.05) +
  geom_errorbarh(aes(xmin=x_mean-x_sd, xmax=x_mean+x_sd), height=0.05) +
  geom_text(aes(label=group), vjust=-1) +
  labs(x="MCOA1", y="MCOA2",
       title="Prophage trajectory (dt-balanced)") +
  theme_bw(base_size=14)

ggsave(file.path(out_dir, "trajectory_mean.png"),
       p, width=6, height=5, dpi=300)

cat("\n🔥 DONE: prophage (dt-balanced, aligned)\n")
