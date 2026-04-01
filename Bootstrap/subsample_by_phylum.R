#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(123)

# =========================
# 路径
# =========================
meta_file <- "/root/mcoa_project/data_raw/genome_and_phylum.csv"
in_dir    <- "/root/mcoa_project/data_raw"
out_dir   <- "/root/mcoa_project/data_raw_subset"

dir.create(out_dir, showWarnings=FALSE)

# =========================
# 读取 meta（⭐固定列名）
# =========================
meta <- fread(meta_file)

# ⭐直接指定列（关键修复）
meta <- meta[, .(
  species = `GTDB species name`,
  phylum  = phylum
)]

# 去掉NA
meta <- meta[!is.na(species) & !is.na(phylum)]

# =========================
# 每个 phylum 抽样
# =========================
selected <- meta[, {
  if (.N <= 100) {
    .SD
  } else {
    .SD[sample(.N, 100)]
  }
}, by = phylum]

selected_species <- unique(selected$species)

cat("[INFO] Selected species:", length(selected_species), "\n")

# =========================
# 筛选所有数据表（按 species）
# =========================
files <- list.files(in_dir, pattern="\\.tsv$", full.names=TRUE)

# 排除meta文件
files <- files[!grepl("genome_and_phylum", files)]

for (f in files){

  dt <- fread(f)

  rn <- dt[[1]]
  dt[[1]] <- NULL

  # ⭐按 species 筛选
  keep_idx <- rn %in% selected_species

  dt_sub <- dt[keep_idx,]

  out_f <- file.path(out_dir, basename(f))

  fwrite(cbind(ID = rn[keep_idx], dt_sub),
         out_f, sep="\t")

  cat("[OK]", basename(f), "→", nrow(dt_sub), "rows\n")
}

# =========================
# 保存抽样ID
# =========================
writeLines(selected_species,
           file.path(out_dir, "selected_species.txt"))

cat("[DONE] Subsampling complete\n")
