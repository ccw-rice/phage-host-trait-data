#!/usr/bin/env Rscript

options(error=function(){traceback(2); quit(status=1)})

suppressPackageStartupMessages({
  library(data.table)
  library(ade4)
  library(vegan)
})

set.seed(123)

# =========================
# 🔥 用你新文件（关键变化）
# =========================
meta_file <- "/root/mcoa_project/data_raw/genome_phylum_dt.csv"

raw_dir   <- "/root/mcoa_project/data_raw"
base_out  <- "/root/mcoa_project/mcoa_bootstrap_dt_balanced"

dir.create(base_out, showWarnings=FALSE)

# =========================
# 工具函数
# =========================
get_rec_cols <- function(X){

  rec_cols <- character(0)

  for (j in seq_along(X)){
    v <- X[[j]]
    v <- v[is.finite(v)]

    if (length(v) == 0) next
    if (!is.numeric(v)) next

    if (all(v >= 0) &&
        (max(v) > 20 ||
         (max(v) > 1 &&
          quantile(v, 0.9) / max(v) < 0.1))){
      rec_cols <- c(rec_cols, colnames(X)[j])
    }
  }

  return(rec_cols)
}

clean_matrix <- function(X){
  X <- X[, colSums(is.na(X)) < nrow(X), drop=FALSE]
  X <- X[rowSums(is.na(X)) == 0, , drop=FALSE]
  return(X)
}

# =========================
# 抽样函数（无ID bug）
# =========================
sample_group_advanced <- function(sub_dt, target_n){

  sub_dt <- copy(sub_dt)

  counts <- sub_dt[, .N, by=phylum]
  P <- nrow(counts)
  quota <- floor(target_n / P)

  sampled_list <- list()

  for (i in 1:P){
    ph <- counts$phylum[i]
    sub <- sub_dt[phylum == ph]

    n_take <- min(nrow(sub), quota)
    sampled_list[[i]] <- sub[sample(.N, n_take)]
  }

  sampled <- rbindlist(sampled_list)

  remaining <- sub_dt[!species %in% sampled$species]

  need <- target_n - nrow(sampled)

  if (need > 0){

    pools <- split(remaining, by="phylum")
    phyla_names <- names(pools)

    idx <- 1
    add_list <- list()

    while (length(add_list) < need){

      ph <- phyla_names[idx]

      if (nrow(pools[[ph]]) > 0){

        pick_idx <- sample(nrow(pools[[ph]]), 1)

        add_list[[length(add_list)+1]] <- pools[[ph]][pick_idx]

        pools[[ph]] <- pools[[ph]][-pick_idx]
      }

      idx <- idx + 1
      if (idx > length(phyla_names)) idx <- 1

      if (all(sapply(pools, nrow) == 0)) break
    }

    if (length(add_list) > 0){
      sampled <- rbind(sampled, rbindlist(add_list))
    }
  }

  if (nrow(sampled) > target_n){
    sampled <- sampled[sample(.N, target_n)]
  }

  return(sampled)
}

# =========================
# 🔥 读取新 meta（关键）
# =========================
meta <- fread(meta_file)

meta <- meta[, .(
  species = `GTDB species name`,
  phylum  = phylum,
  dt      = `Doubling time (hour)`
)]

meta <- meta[is.finite(dt)]

meta[, group := ifelse(dt < 5, "<5h", ">5h")]

# =========================
# 主循环
# =========================
for (i in 1:50){

  cat("\n========== RUN", i, "==========\n")

  set.seed(1000 + i)

  run_dir <- file.path(base_out, sprintf("run_%02d", i))
  dir.create(run_dir, showWarnings=FALSE, recursive=TRUE)

  ln_dir <- file.path(run_dir, "ln1p")
  dir.create(ln_dir, showWarnings=FALSE)

  # =========================
  # 抽样
  # =========================
  fast <- sample_group_advanced(meta[group=="<5h"], 500)
  slow <- sample_group_advanced(meta[group==">5h"], 500)

  selected <- rbindlist(list(fast, slow))
  selected_species <- selected$species

  cat("Selected:", length(selected_species), "\n")

  files <- list.files(raw_dir, pattern="\\.tsv$", full.names=TRUE)

  ln_data <- list()

  # =========================
  # 后续完全不变
  # =========================
  for (f in files){

    dt <- fread(file = f)

    rn <- dt[[1]]
    dt[[1]] <- NULL

    keep <- rn %in% selected_species
cat("\nChecking matching...\n")

cat("Selected species:", length(selected_species), "\n")
cat("Matched species :", sum(keep), "\n")

missing_species <- setdiff(selected_species, rn)

cat("Missing species :", length(missing_species), "\n")

if(length(missing_species) > 0){
    print(head(missing_species, 100))
}
    X <- as.data.frame(dt[keep,])
    rownames(X) <- rn[keep]

    X[] <- lapply(X, function(x)
      suppressWarnings(as.numeric(as.character(x))))

    rec_cols <- get_rec_cols(X)

    for (col in rec_cols){
      X[[col]] <- log1p(X[[col]])
    }

    base_name <- basename(f)

    fwrite(
      cbind(ID=rownames(X), X),
      file.path(ln_dir, base_name),
      sep="\t"
    )

    ln_data[[base_name]] <- X
  }

std <- function(X) vegan::decostand(X, method="standardize")

# ======================================
# DEBUG SECTION
# ======================================

X_KO_raw     <- ln_data[["Ko.merged_mean.tsv"]]
X_COG_raw    <- ln_data[["Cog.merged_mean.tsv"]]
X_CAZy_raw   <- ln_data[["Cazy.merged_mean.tsv"]]
X_TRAITS_raw <- ln_data[["Traits.merged_mean.tsv"]]

cat("\n==============================\n")
cat("DEBUG REPORT\n")
cat("==============================\n")

cat("\nSelected species:\n")
cat(length(selected_species), "\n")

cat("\nRows BEFORE clean_matrix:\n")
cat("KO     :", nrow(X_KO_raw), "\n")
cat("COG    :", nrow(X_COG_raw), "\n")
cat("CAZY   :", nrow(X_CAZy_raw), "\n")
cat("TRAITS :", nrow(X_TRAITS_raw), "\n")

# --------------------------------------
# NA statistics
# --------------------------------------

cat("\nRows containing NA:\n")

cat(
  "KO     :",
  sum(rowSums(is.na(X_KO_raw)) > 0),
  "\n"
)

cat(
  "COG    :",
  sum(rowSums(is.na(X_COG_raw)) > 0),
  "\n"
)

cat(
  "CAZY   :",
  sum(rowSums(is.na(X_CAZy_raw)) > 0),
  "\n"
)

cat(
  "TRAITS :",
  sum(rowSums(is.na(X_TRAITS_raw)) > 0),
  "\n"
)

# --------------------------------------
# clean matrix
# --------------------------------------

X_KO     <- std(clean_matrix(X_KO_raw))
X_COG    <- std(clean_matrix(X_COG_raw))
X_CAZy   <- std(clean_matrix(X_CAZy_raw))
X_TRAITS <- std(clean_matrix(X_TRAITS_raw))

cat("\nRows AFTER clean_matrix:\n")
cat("KO     :", nrow(X_KO), "\n")
cat("COG    :", nrow(X_COG), "\n")
cat("CAZY   :", nrow(X_CAZy), "\n")
cat("TRAITS :", nrow(X_TRAITS), "\n")

# --------------------------------------
# lost species
# --------------------------------------

lost_ko <- setdiff(selected_species, rownames(X_KO))
lost_cog <- setdiff(selected_species, rownames(X_COG))
lost_cazy <- setdiff(selected_species, rownames(X_CAZy))
lost_traits <- setdiff(selected_species, rownames(X_TRAITS))

cat("\nSpecies removed:\n")
cat("KO     :", length(lost_ko), "\n")
cat("COG    :", length(lost_cog), "\n")
cat("CAZY   :", length(lost_cazy), "\n")
cat("TRAITS :", length(lost_traits), "\n")

# --------------------------------------
# intersection
# --------------------------------------

ids <- Reduce(intersect, list(
  rownames(X_KO),
  rownames(X_COG),
  rownames(X_CAZy),
  rownames(X_TRAITS)
))

cat("\nCommon species after intersection:\n")
cat(length(ids), "\n")

cat("\nSpecies lost by intersection:\n")
cat(length(selected_species) - length(ids), "\n")

# --------------------------------------
# subset
# --------------------------------------

X_KO     <- X_KO[ids,]
X_COG    <- X_COG[ids,]
X_CAZy   <- X_CAZy[ids,]
X_TRAITS <- X_TRAITS[ids,]

cat("\nFinal matrix dimensions:\n")

cat(
  "KO     :",
  nrow(X_KO),
  "x",
  ncol(X_KO),
  "\n"
)

cat(
  "COG    :",
  nrow(X_COG),
  "x",
  ncol(X_COG),
  "\n"
)

cat(
  "CAZY   :",
  nrow(X_CAZy),
  "x",
  ncol(X_CAZy),
  "\n"
)

cat(
  "TRAITS :",
  nrow(X_TRAITS),
  "x",
  ncol(X_TRAITS),
  "\n"
)

cat("\n==============================\n")
cat("END DEBUG REPORT\n")
cat("==============================\n\n")
acp_KO     <- dudi.pca(X_KO, scannf=FALSE, nf=3)
acp_COG    <- dudi.pca(X_COG, scannf=FALSE, nf=3)
acp_CAZy   <- dudi.pca(X_CAZy, scannf=FALSE, nf=3)
acp_TRAITS <- dudi.pca(X_TRAITS, scannf=FALSE, nf=3)

acom <- mcoa(ktab.list.dudi(list(
  acp_KO, acp_COG, acp_CAZy, acp_TRAITS
)), scannf=FALSE)

# =========================
# 保存 Tco
# =========================

Tco <- as.data.table(acom$Tco)

setnames(
  Tco,
  paste0(
    "SV",
    seq_len(ncol(Tco))
  )
)

Tco[, var := rownames(acom$Tco)]

# 与 samples.tsv 保持方向一致
if("SV1" %in% colnames(Tco))
  Tco[, SV1 := -SV1]

if("SV2" %in% colnames(Tco))
  Tco[, SV2 := -SV2]

fwrite(
  Tco,
  file.path(
    run_dir,
    "variables_contribution_Tco.tsv"
  ),
  sep="\t"
)

# =========================
# 保存 samples
# =========================

S <- as.data.frame(acom$SynVar)

colnames(S) <- paste0(
  "SynVar",
  seq_len(ncol(S))
)

S$id <- rownames(S)

S$SynVar1 <- -S$SynVar1
S$SynVar2 <- -S$SynVar2

fwrite(
  S,
  file.path(
    run_dir,
    "samples.tsv"
  ),
  sep="\t"
)
  cat("Run", i, "DONE\n")
}

cat("\n🔥 ALL RUNS COMPLETE\n")
