#!/usr/bin/env Rscript

options(error=function(){traceback(2); quit(status=1)})

suppressPackageStartupMessages({
  library(data.table)
  library(ade4)
  library(vegan)
})

set.seed(123)

meta_file <- "/root/mcoa_project/data_raw/genome_and_phylum.csv"
raw_dir   <- "/root/mcoa_project/data_raw"
base_out  <- "/root/mcoa_project/mcoa_bootstrap"

dir.create(base_out, showWarnings=FALSE)

# =========================
# clean_names
# =========================
clean_names <- function(x){
  x <- tolower(x)
  x <- sub("^x(?=[0-9])", "", x, perl=TRUE)
  x <- gsub("\\(", "", x)
  x <- gsub("\\)", "", x)
  x <- gsub("[^a-z0-9]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  return(x)
}

# =========================
# 动态判断 ln1p
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

# =========================
# NA清理
# =========================
clean_matrix <- function(X){
  X <- X[, colSums(is.na(X)) < nrow(X), drop=FALSE]
  X <- X[rowSums(is.na(X)) == 0, , drop=FALSE]
  return(X)
}

# =========================
# 读取 meta
# =========================
meta <- fread(file = meta_file)

meta <- meta[, .(
  species = `GTDB species name`,
  phylum  = phylum
)]

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
  # 1️⃣ subsample
  # =========================
  selected <- meta[, {
    if (.N <= 100) .SD else .SD[sample(.N, 100)]
  }, by=phylum]

  selected_species <- selected$species

  files <- list.files(raw_dir, pattern="\\.tsv$", full.names=TRUE)

  ln_data <- list()

  # =========================
  # 2️⃣ dynamic ln1p + 保存
  # =========================
  for (f in files){

    dt <- fread(file = f)

    rn <- dt[[1]]
    dt[[1]] <- NULL

    keep <- rn %in% selected_species

    X <- as.data.frame(dt[keep,])
    rownames(X) <- rn[keep]

    # ⭐ numeric 转换
    X[] <- lapply(X, function(x)
      suppressWarnings(as.numeric(as.character(x))))

    # ⭐ dynamic 判断
    rec_cols <- get_rec_cols(X)

    # ⭐ 应用 ln1p
    for (col in rec_cols){
      X[[col]] <- log1p(X[[col]])
    }

    # =========================
    # ⭐保存 ln1p 数据
    # =========================
    base_name <- basename(f)

    fwrite(
      cbind(ID=rownames(X), X),
      file.path(ln_dir, base_name),
      sep="\t"
    )

    # ⭐记录哪些列被转换
    writeLines(
      rec_cols,
      file.path(ln_dir, paste0(base_name, ".ln1p_cols.txt"))
    )

    ln_data[[base_name]] <- X
  }

  # =========================
  # 3️⃣ MCOA
  # =========================
  std <- function(X) vegan::decostand(X, method="standardize")

  X_KO     <- std(clean_matrix(ln_data[["Ko.merged_mean.tsv"]]))
  X_COG    <- std(clean_matrix(ln_data[["Cog.merged_mean.tsv"]]))
  X_CAZy   <- std(clean_matrix(ln_data[["Cazy.merged_mean.tsv"]]))
  X_TRAITS <- std(clean_matrix(ln_data[["Traits.merged_mean.tsv"]]))

  ids <- Reduce(intersect, list(
    rownames(X_KO), rownames(X_COG),
    rownames(X_CAZy), rownames(X_TRAITS)
  ))

  X_KO     <- X_KO[ids,]
  X_COG    <- X_COG[ids,]
  X_CAZy   <- X_CAZy[ids,]
  X_TRAITS <- X_TRAITS[ids,]

  acp_KO     <- dudi.pca(X_KO, scannf=FALSE, nf=3)
  acp_COG    <- dudi.pca(X_COG, scannf=FALSE, nf=3)
  acp_CAZy   <- dudi.pca(X_CAZy, scannf=FALSE, nf=3)
  acp_TRAITS <- dudi.pca(X_TRAITS, scannf=FALSE, nf=3)

  acom <- mcoa(ktab.list.dudi(list(
    acp_KO, acp_COG, acp_CAZy, acp_TRAITS
  )), scannf=FALSE)

  # =========================
  # 4️⃣ samples
  # =========================
  S <- as.data.frame(acom$SynVar)
  colnames(S) <- paste0("SynVar", seq_len(ncol(S)))
  S$id <- rownames(S)

  S$SynVar1 <- -S$SynVar1
  S$SynVar2 <- -S$SynVar2

  fwrite(S, file.path(run_dir,"samples.tsv"), sep="\t")

  # =========================
  # 5️⃣ Tco
  # =========================
  Tco <- as.data.table(acom$Tco)
  setnames(Tco, paste0("SV", seq_len(ncol(Tco))))
  Tco[, var := rownames(acom$Tco)]

  Tco[, SV1 := -SV1]
  Tco[, SV2 := -SV2]

  fwrite(Tco, file.path(run_dir,"variables_contribution_Tco.tsv"), sep="\t")

  # =========================
  # 6️⃣ LM
  # =========================
  DF <- data.frame(cbind(X_KO, X_COG, X_CAZy, X_TRAITS))

  vars <- colnames(DF)

  out <- data.frame(
    var=vars,
    coef1=NA, p1=NA, R2_1=NA,
    coef2=NA, p2=NA, R2_2=NA
  )

  for (j in seq_along(vars)){
    y <- DF[[vars[j]]]

    sm1 <- summary(lm(y ~ S$SynVar1))
    out$coef1[j] <- sm1$coefficients[2,1]
    out$p1[j]    <- sm1$coefficients[2,4]
    out$R2_1[j]  <- sm1$r.squared

    sm2 <- summary(lm(y ~ S$SynVar2))
    out$coef2[j] <- sm2$coefficients[2,1]
    out$p2[j]    <- sm2$coefficients[2,4]
    out$R2_2[j]  <- sm2$r.squared
  }

  fwrite(out, file.path(run_dir,"variables_vs_axes_lm.tsv"), sep="\t")

  # =========================
  # 7️⃣ 筛选
  # =========================
  out$var <- clean_names(out$var)
  Tco$var <- clean_names(Tco$var)

  lm_long <- rbindlist(list(
    data.table(var=out$var, axis=1, coef=out$coef1, pval=out$p1, r2=out$R2_1),
    data.table(var=out$var, axis=2, coef=out$coef2, pval=out$p2, r2=out$R2_2)
  ))

  tco_long <- melt(Tco[,.(var,SV1,SV2)],
                   id.vars="var",
                   variable.name="axis",
                   value.name="SV")

  tco_long[,axis:=as.integer(sub("SV","",axis))]

  df <- merge(lm_long,tco_long,by=c("var","axis"))

  for (ax in 1:2){
    sub <- df[axis==ax & pval<1e-3 & r2>0.2]
    sub <- sub[order(-abs(SV))]

    fwrite(head(sub,50),
           file.path(run_dir,paste0("axis",ax,".tsv")),
           sep="\t")
  }

  cat("Run", i, "DONE\n")
}

cat("\n🔥 ALL 50 RUNS COMPLETE\n")
