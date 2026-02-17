#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(dplyr) })

# --- 配置 ---
in_lm   <- "variables_vs_axes_lm.tsv"          # 宽表：var, coef1,p1,R2_1, coef2,p2,R2_2, coef3,p3,R2_3
in_tco  <- "variables_contribution_Tco.tsv"    # 宽表：SV1,SV2,SV3,..., var, group
out_dir <- "best_vars_filtered"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

p_cut  <- 1e-3
r2_cut <- 0.2
topN   <- 50

# --- 读入 ---
lm_wide <- fread(in_lm)
tco_wide <- fread(in_tco)

stopifnot("var" %in% names(lm_wide))
stopifnot("var" %in% names(tco_wide))

# --- 1) 自动检测回归结果中的轴（支持 coef1 / coef_cor_MCOA1 两种命名）---
axes_map <- list()
nms <- names(lm_wide)

coef_cols1 <- grep("^coef[0-9]+$", nms, value = TRUE)
if (length(coef_cols1) > 0){
  for (cc in coef_cols1){
    ax <- as.integer(sub("^coef", "", cc))
    pc <- paste0("p", ax)
    rc <- paste0("R2_", ax)
    if (pc %in% nms && rc %in% nms){
      axes_map[[as.character(ax)]] <- list(coef=cc, p=pc, r2=rc)
    }
  }
}

coef_cols2 <- grep("^coef_cor_MCOA[0-9]+$", nms, value = TRUE)
if (length(coef_cols2) > 0){
  for (cc in coef_cols2){
    ax <- as.integer(sub("^coef_cor_MCOA", "", cc))
    pc <- paste0("pval_cor_MCOA", ax)
    rc <- paste0("R2_cor_MCOA", ax)
    if (pc %in% nms && rc %in% nms){
      axes_map[[as.character(ax)]] <- list(coef=cc, p=pc, r2=rc)
    }
  }
}

if (length(axes_map) == 0) {
  stop("Cannot detect axis columns in ", in_lm)
}
axes <- sort(as.integer(names(axes_map)))
message("[INFO] Detected axes in LM: ", paste(axes, collapse=", "))

# --- 2) 宽 -> 长：回归结果（得到 var/axis/coef/pval/r2）---
lm_long <- rbindlist(lapply(axes, function(ax){
  cols <- axes_map[[as.character(ax)]]
  data.frame(
    var  = lm_wide$var,
    axis = ax,
    coef = lm_wide[[cols$coef]],
    pval = lm_wide[[cols$p]],
    r2   = lm_wide[[cols$r2]],
    stringsAsFactors = FALSE
  )
}), fill = TRUE)

# --- 3) 宽 -> 长：Tco（把 SV1,SV2,... 融合成 var/axis/SV）---
sv_cols <- grep("^SV[0-9]+$", names(tco_wide), value = TRUE)
if (length(sv_cols) == 0) stop("No SV1/SV2/... columns in ", in_tco)

# data.table::melt：把每个 SV 列融成长表
tco_long <- melt(
  tco_wide[, c("var", sv_cols), with = FALSE],
  id.vars = "var",
  variable.name = "sv_name",
  value.name = "SV"
)
# sv_name: "SV1" -> 1
tco_long$axis <- as.integer(sub("^SV", "", tco_long$sv_name))
tco_long$sv_name <- NULL

message("[INFO] Detected axes in Tco: ", paste(sort(unique(tco_long$axis)), collapse=", "))

# --- 4) 合并：按 var & axis 对齐，得到 df(var, axis, coef, pval, r2, SV) ---
df <- lm_long %>% left_join(tco_long, by = c("var", "axis"))

# --- 5) 筛选函数（AND：p < p_cut & r2 > r2_cut），再按 |SV| 取前 topN ---
filter_top <- function(axis_id, topN=50, r2_cut=0.2, p_cut=1e-3){
  df %>%
    filter(axis == axis_id) %>%
    filter(is.finite(pval), is.finite(r2), is.finite(SV)) %>%
    filter(pval < p_cut, r2 > r2_cut) %>%
    arrange(desc(abs(SV))) %>%
    slice_head(n = topN)
}

# --- 6) 遍历输出 ---
summary_list <- list()
for (ax in axes){
  res <- filter_top(ax, topN, r2_cut, p_cut)
  out_file <- file.path(out_dir, sprintf("best_vars_axis%d_filtered.tsv", ax))
  fwrite(res, out_file, sep = "\t")
  message("[OK] Wrote: ", out_file, " (", nrow(res), " rows)")
  summary_list[[as.character(ax)]] <- data.frame(axis=ax,
                                                 kept=nrow(res),
                                                 total=sum(df$axis==ax))
}
summary_df <- bind_rows(summary_list)
fwrite(summary_df, file.path(out_dir, "filter_summary.tsv"), sep="\t")
message("[DONE] Summary -> ", file.path(out_dir, "filter_summary.tsv"))
