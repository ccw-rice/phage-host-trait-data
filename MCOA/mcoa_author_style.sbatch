#!/bin/bash
#SBATCH -J mcoa_author_style
#SBATCH -p commons
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH -o /projects/alvarez/scripts/logs/mcoa_author_style.%A.out
#SBATCH -e /projects/alvarez/scripts/logs/mcoa_author_style.%A.err

set -euo pipefail
set -x

echo "[INFO] start at: $(date)"
hostname
whoami

module purge
export MAMBA_ROOT_PREFIX=/scratch/cw106/software
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook --shell bash)"
micromamba activate r_ade4

IN_DIR="/scratch/cw106/MCOA_origin_preprocessed2/MCOA_origin_data"
ls -lh "$IN_DIR" || true

R --vanilla <<'RS'
options(error=function(){traceback(2); quit(status=1)})

suppressPackageStartupMessages({
  pkgs <- c("ade4","vegan","tibble","ggplot2","ggpubr","viridis","data.table")
  for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) stop("Missing package: ", p)
  lapply(pkgs, library, character.only=TRUE)
})

in_dir  <- "/scratch/cw106/MCOA_origin_preprocessed2/MCOA_origin_data"
out_dir <- file.path(in_dir, "mcoa_author_style")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# 若 ln1p 文件不存在则回退到原始文件
fallback <- function(f_ln1p, f_raw){
  p1 <- file.path(in_dir, f_ln1p)
  if (file.exists(p1)) return(p1)
  message("[WARN] Not found ln1p: ", p1, " -> use raw: ", f_raw)
  file.path(in_dir, f_raw)
}

f_KO     <- fallback("Ko.merged_mean.ln1p.tsv",    "Ko.merged_mean.tsv")
f_COG    <- fallback("Cog.merged_mean.ln1p.tsv",   "Cog.merged_mean.tsv")
f_CAZy   <- fallback("Cazy.merged_mean.ln1p.tsv",  "Cazy.merged_mean.tsv")
f_TRAITS <- fallback("Traits.merged_mean.ln1p.tsv","Traits.merged_mean.tsv")

print(data.frame(KO=f_KO, COG=f_COG, CAZy=f_CAZy, TRAITS=f_TRAITS), row.names=FALSE)

read_tsv_matrix <- function(path){
  dt <- data.table::fread(path)
  rn <- dt[[1]]; dt[[1]] <- NULL
  X <- as.data.frame(dt); rownames(X) <- rn
  X[] <- lapply(X, function(x) suppressWarnings(as.numeric(as.character(x))))
  if (any(!is.finite(as.matrix(X)))) stop("Non-finite in ", basename(path))
  X
}
std_vegan <- function(X) vegan::decostand(X, method="standardize")

# 读、标准化
X_KO     <- std_vegan(read_tsv_matrix(f_KO))
X_COG    <- std_vegan(read_tsv_matrix(f_COG))
X_CAZy   <- std_vegan(read_tsv_matrix(f_CAZy))
X_TRAITS <- std_vegan(read_tsv_matrix(f_TRAITS))

# 行交集
ids <- Reduce(intersect, list(rownames(X_KO), rownames(X_COG), rownames(X_CAZy), rownames(X_TRAITS)))
if (length(ids)==0) stop("No common IDs across tables.")
X_KO <- X_KO[ids,,drop=FALSE]; X_COG <- X_COG[ids,,drop=FALSE]
X_CAZy <- X_CAZy[ids,,drop=FALSE]; X_TRAITS <- X_TRAITS[ids,,drop=FALSE]
message("[INFO] intersection size: ", length(ids))

# 每表 PCA (nf=3)
acp_KO     <- ade4::dudi.pca(as.data.frame(X_KO),     scannf=FALSE, nf=3)
acp_COG    <- ade4::dudi.pca(as.data.frame(X_COG),    scannf=FALSE, nf=3)
acp_CAZy   <- ade4::dudi.pca(as.data.frame(X_CAZy),   scannf=FALSE, nf=3)
acp_TRAITS <- ade4::dudi.pca(as.data.frame(X_TRAITS), scannf=FALSE, nf=3)

# 组装并 MCOA
ktab <- ade4::ktab.list.dudi(list(acp_KO, acp_COG, acp_CAZy, acp_TRAITS),
                             tabnames=c("KO","COG","CAZy","Traits"))
acom <- ade4::mcoa(ktab, option="inertia", scannf=FALSE)

# 方差份额
var_share <- acom$pseudoeig/sum(acom$pseudoeig)
data.table::fwrite(data.frame(axis=seq_along(var_share), var_share),
                   file.path(out_dir,"axis_variation_share.tsv"), sep="\t")

# 变量贡献（Tco）
Tco <- as.data.frame(acom$Tco); colnames(Tco) <- paste0("SV", seq_len(ncol(Tco)))
Tco$var <- rownames(Tco)
groups <- acom$TC[,1]; names(groups) <- rownames(acom$TC)
Tco$group <- groups[match(Tco$var, names(groups))]
data.table::fwrite(Tco, file.path(out_dir,"variables_contribution_Tco.tsv"), sep="\t")

# 样本坐标
S <- as.data.frame(acom$SynVar); colnames(S) <- paste0("SynVar", seq_len(ncol(S))); S$id <- rownames(S)
data.table::fwrite(S[,c("id","SynVar1","SynVar2","SynVar3")], file.path(out_dir,"samples_SynVar123.tsv"), sep="\t")

# 变量~轴 显著性 (lm)
DF <- data.frame(cbind(X_KO, X_COG, X_CAZy, X_TRAITS))
stopifnot(identical(rownames(DF), rownames(S)))
Syn1 <- acom$SynVar[,"SynVar1"]; Syn2 <- acom$SynVar[,"SynVar2"]; Syn3 <- acom$SynVar[,"SynVar3"]

vars <- colnames(DF)
cor_tab <- data.frame(var=vars, coef1=NA_real_, p1=NA_real_, R2_1=NA_real_,
                      coef2=NA_real_, p2=NA_real_, R2_2=NA_real_,
                      coef3=NA_real_, p3=NA_real_, R2_3=NA_real_)
for (i in seq_along(vars)){
  y <- as.numeric(DF[[vars[i]]])
  sm1 <- summary(lm(y ~ Syn1)); cor_tab$coef1[i] <- sm1$coefficients["Syn1","Estimate"]; cor_tab$p1[i] <- sm1$coefficients["Syn1","Pr(>|t|)"]; cor_tab$R2_1[i] <- sm1$r.squared
  sm2 <- summary(lm(y ~ Syn2)); cor_tab$coef2[i] <- sm2$coefficients["Syn2","Estimate"]; cor_tab$p2[i] <- sm2$coefficients["Syn2","Pr(>|t|)"]; cor_tab$R2_2[i] <- sm2$r.squared
  sm3 <- summary(lm(y ~ Syn3)); cor_tab$coef3[i] <- sm3$coefficients["Syn3","Estimate"]; cor_tab$p3[i] <- sm3$coefficients["Syn3","Pr(>|t|)"]; cor_tab$R2_3[i] <- sm3$r.squared
}
data.table::fwrite(cor_tab, file.path(out_dir,"variables_vs_axes_lm.tsv"), sep="\t")

best1 <- cor_tab$var[cor_tab$p1 < 1e-3]
best2 <- cor_tab$var[cor_tab$p2 < 1e-3]
best3 <- cor_tab$var[cor_tab$p3 < 1e-3]
writeLines(best1, file.path(out_dir,"best_vars_axis1_Plt0.001.txt"))
writeLines(best2, file.path(out_dir,"best_vars_axis2_Plt0.001.txt"))
writeLines(best3, file.path(out_dir,"best_vars_axis3_Plt0.001.txt"))

# 简单散点，确认走完整个流程
pdf(file.path(out_dir,"MCOA_samples_SynVar1_vs_2.pdf"), width=7.5, height=6)
plot(acom$SynVar[,"SynVar1"], acom$SynVar[,"SynVar2"], pch=16, cex=0.6,
     xlab="MCOA1", ylab="MCOA2", main="MCOA samples (SynVar)")
dev.off()

message("[DONE] Outputs in: ", out_dir)
RS

echo "[INFO] end at: $(date)"
