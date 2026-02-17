#!/bin/bash
#SBATCH -J mt_batch_v3i
#SBATCH -p commons
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH -o /scratch/%u/logs/mt.batch.%A_%a.out
#SBATCH -e /scratch/%u/logs/mt.batch.%A_%a.err

set -Eeuo pipefail

# === 可通过 --export 覆盖的参数 ===
ENV=${ENV:-/scratch/$USER/envs/microtrait_env}
HMM_DIR=${HMM_DIR:-/scratch/$USER/microtrait-db}
MANIFEST=${MANIFEST:-/scratch/$USER/microtrait_manifest.gcf.tsv}
OUTROOT=${OUTROOT:-/scratch/$USER/microtrait_batch_gcf}

# 阈值
COV_GENE=${COV_GENE:-50}
COV_HMM=${COV_HMM:-55}
IE_OR=${IE_OR:-1e-5}
PROT_G=${PROT_G:-45}
PROT_H=${PROT_H:-50}

HMMSCAN="$ENV/bin/hmmscan"
RSCRIPT="$ENV/bin/Rscript"

# === 读取 manifest：两列 <sample>\t<protein.faa> ===
idx=${SLURM_ARRAY_TASK_ID:-1}
read sample faa <<<"$(awk -v i="$idx" 'BEGIN{FS="\t"} $0!~/^#/ && NF && NR==i {print $1, $2}' "$MANIFEST")" || true

if [[ -z "${sample:-}" || -z "${faa:-}" ]]; then
  echo "[FATAL] manifest line $idx invalid." >&2; exit 1
fi
if [[ ! -s "$faa" ]]; then
  echo "[FATAL] FAA not found or empty: $faa" >&2; exit 1
fi

OUT="$OUTROOT/$sample"
mkdir -p "$OUT/sanitized" "$OUT/traits"

echo "[INFO] sample=$sample"
echo "[INFO] faa=$faa"
echo "[INFO] out=$OUT"
echo "[INFO] env hmmscan: $HMMSCAN"

# === microtrait HMM DB ===
MICRO_DB="$HMM_DIR/hmm/microtrait-hmmdb/microtrait.hmmdb"
[[ -s "$MICRO_DB" ]] || { echo "[FATAL] microtrait DB not found: $MICRO_DB" >&2; exit 1; }

CUT_MICRO=""
if grep -q -m1 '^TC' "$MICRO_DB"; then
  CUT_MICRO="--cut_tc"
  echo "[INFO] microtrait DB has TC -> --cut_tc"
else
  echo "[WARN] microtrait DB has no TC -> no --cut_tc"
fi

# === dbCAN DB（可选）===
DBCAN_DB=""
if [[ -s "$HMM_DIR/hmm/dbcan/dbcan.noGA.hmmdb" ]]; then
  DBCAN_DB="$HMM_DIR/hmm/dbcan/dbcan.noGA.hmmdb"
elif [[ -s "$HMM_DIR/hmm/dbcan/dbcan.select.v8.hmmdb" ]]; then
  DBCAN_DB="$HMM_DIR/hmm/dbcan/dbcan.select.v8.hmmdb"
else
  echo "[WARN] dbCAN DB not found; skip dbCAN."
fi

# === 1) hmmscan: microtrait ===
"$HMMSCAN" $CUT_MICRO --cpu "${SLURM_CPUS_PER_TASK:-8}" --noali \
  --domtblout "$OUT/microtrait.domtblout" \
  "$MICRO_DB" "$faa" > "$OUT/microtrait.hmmscan.log" 2>&1

# === 2) hmmscan: dbCAN（可选）===
if [[ -n "$DBCAN_DB" ]]; then
  "$HMMSCAN" --cpu "${SLURM_CPUS_PER_TASK:-8}" --noali \
    --domtblout "$OUT/dbcan.domtblout" \
    "$DBCAN_DB" "$faa" > "$OUT/dbcan.hmmscan.log" 2>&1 || true
fi

# （可选）检查 dbCAN 日志
if [[ -e "$OUT/dbcan.hmmscan.log" ]] && ! grep -q "CPU time" "$OUT/dbcan.hmmscan.log"; then
  echo "[WARN] hmmscan may not have finished cleanly for dbCAN: $sample"
fi

# === 3) 规范化 domtblout 为 23 列 TSV ===
awk 'BEGIN{OFS="\t"} !/^#/ {d=$23; for(i=24;i<=NF;i++) d=d" "$i;
  print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,d
}' "$OUT/microtrait.domtblout" > "$OUT/sanitized/microtrait.domtblout.tsv"

if [[ -s "$OUT/dbcan.domtblout" ]]; then
  awk 'BEGIN{OFS="\t"} !/^#/ {d=$23; for(i=24;i<=NF;i++) d=d" "$i;
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,d
  }' "$OUT/dbcan.domtblout" > "$OUT/sanitized/dbcan.domtblout.tsv"
else
  : > "$OUT/sanitized/dbcan.domtblout.tsv"
fi

# === 4) R：并集覆盖 + 阈值；补全 hmm_qlen；空 dbCAN 用空表兜底 ===
"$RSCRIPT" --vanilla - "$OUT" "$COV_GENE" "$COV_HMM" "$IE_OR" "$PROT_G" "$PROT_H" <<'RS'
suppressPackageStartupMessages({
  library(data.table)
  library(microtrait)
})

args <- commandArgs(trailingOnly=TRUE)
out   <- args[1]
cov_g <- as.numeric(args[2]); cov_h <- as.numeric(args[3])
ie_or <- as.numeric(args[4]); prot_g <- as.numeric(args[5]); prot_h <- as.numeric(args[6])

san <- file.path(out, "sanitized")

cn <- c("target_name","target_accession","tlen",
        "query_name","query_accession","qlen",
        "full_evalue","full_score","full_bias",
        "dom_num","dom_of","c_evalue","i_evalue","dom_score","dom_bias",
        "hmm_from","hmm_to","ali_from","ali_to","env_from","env_to",
        "acc","target_description")

# 读取 microtrait
mt <- fread(file.path(san,"microtrait.domtblout.tsv"), sep="\t", header=FALSE, col.names=cn, quote="")

# 读取 dbCAN：存在且 size>0 才读
dbfile <- file.path(san,"dbcan.domtblout.tsv")
db <- if (file.exists(dbfile) && file.info(dbfile)$size > 0) {
  fread(dbfile, sep="\t", header=FALSE, col.names=cn, quote="")
} else {
  NULL
}

numi <- function(v) suppressWarnings(as.integer(v))
numd <- function(v) suppressWarnings(as.numeric(v))

norm <- function(x){
  setDT(x)
  x[, `:=`(
    hmm_name   = target_name,
    gene_name  = query_name,
    ali_from   = numi(ali_from), ali_to = numi(ali_to),
    hmm_from   = numi(hmm_from), hmm_to = numi(hmm_to),
    gene_len   = numi(qlen),     hmm_len= numi(tlen),
    i_evalue   = numd(i_evalue), hmm_score = numd(dom_score),
    gene_score = numd(dom_score)
  )]
  x <- x[is.finite(ali_from) & is.finite(ali_to) & is.finite(hmm_from) & is.finite(hmm_to)]
  x <- x[gene_len > 0 & hmm_len > 0]
  x[, `:=`(
    g1 = pmin(ali_from, ali_to), g2 = pmax(ali_from, ali_to),
    h1 = pmin(hmm_from, hmm_to), h2 = pmax(hmm_from, hmm_to)
  )]
  x[]
}

mt <- norm(mt)
if (!is.null(db) && nrow(db) > 0) db <- norm(db) else db <- NULL

union_len <- function(starts, ends){
  ok <- is.finite(starts) & is.finite(ends)
  if (!any(ok)) return(0L)
  s <- as.integer(starts[ok]); e <- as.integer(ends[ok])
  w <- which(e < s); if (length(w)) { tmp <- s[w]; s[w] <- e[w]; e[w] <- tmp[w] }
  if (!length(s)) return(0L)
  o <- order(s, e); s <- s[o]; e <- e[o]
  tot <- 0L; cs <- s[1]; ce <- e[1]
  for (i in seq_along(s)[-1]) {
    if (s[i] <= ce + 1L) { if (e[i] > ce) ce <- e[i] } else { tot <- tot + (ce - cs + 1L); cs <- s[i]; ce <- e[i] }
  }
  tot + (ce - cs + 1L)
}

summ <- function(x){
  x[, .(
    gene_len   = max(gene_len, na.rm=TRUE),
    hmm_len    = max(hmm_len,  na.rm=TRUE),
    gene_cov_u = 100 * union_len(g1, g2) / pmax(1, max(gene_len, na.rm=TRUE)),
    hmm_cov_u  = 100 * union_len(h1, h2) / pmax(1, max(hmm_len,  na.rm=TRUE)),
    i_evalue   = suppressWarnings(min(i_evalue, na.rm=TRUE)),
    hmm_score  = suppressWarnings(max(hmm_score, na.rm=TRUE)),
    gene_score = suppressWarnings(max(gene_score, na.rm=TRUE)),
    gene_from  = suppressWarnings(min(g1, na.rm=TRUE)),
    gene_to    = suppressWarnings(max(g2, na.rm=TRUE)),
    hmm_from   = suppressWarnings(min(h1, na.rm=TRUE)),
    hmm_to     = suppressWarnings(max(h2, na.rm=TRUE))
  ), by=.(gene_name, hmm_name)]
}

mts <- summ(mt)
dbs <- if (is.null(db)) NULL else summ(db)

# 关键修复：为两侧都补上 hmm_qlen（microtrait 侧也需要）
mts[, hmm_qlen := hmm_len]
if (!is.null(dbs)) dbs[, hmm_qlen := hmm_len]

# 阈值过滤
is_prot <- grepl("^(pep|lon|clp|hsl|ftsH|degP|prt|pepT)", mts$hmm_name, ignore.case=TRUE)
keep_mt <- ( (mts$gene_cov_u >= cov_g & mts$hmm_cov_u >= cov_h) | (mts$i_evalue <= ie_or) )
keep_mt[is_prot] <- ( (mts$gene_cov_u[is_prot] >= prot_g & mts$hmm_cov_u[is_prot] >= prot_h) |
                      (mts$i_evalue[is_prot] <= ie_or) )
mts_f <- mts[keep_mt]

dbs_f <- NULL
if (!is.null(dbs) && nrow(dbs) > 0) {
  dbs_f <- dbs[ (gene_cov_u >= pmin(45, cov_g) & hmm_cov_u >= pmin(45, cov_h)) | (i_evalue <= ie_or) ]
}

# 兜底空表，并确保列齐全
ensure_cols <- function(dt){
  if (is.null(dt)) return(NULL)
  if (!("hmm_qlen"   %in% names(dt))) dt[, hmm_qlen := hmm_len]
  if (!("gene_score" %in% names(dt))) dt[, gene_score := hmm_score]
  dt[]
}

mts_f <- ensure_cols(mts_f)

empty_dbcan <- data.table(
  gene_name=character(), hmm_name=character(),
  gene_len=integer(), hmm_len=integer(), hmm_qlen=integer(),
  gene_from=integer(), gene_to=integer(),
  hmm_from=integer(), hmm_to=integer(),
  i_evalue=double(), hmm_score=double(), gene_score=double()
)
dbs_f <- ensure_cols(dbs_f)
if (is.null(dbs_f) || nrow(dbs_f) == 0L) dbs_f <- empty_dbcan

# 映射到 trait
res <- microtrait:::map.traits.fromdomtblout(mts_f, dbs_f)

dir.create(file.path(out,"traits"), showWarnings=FALSE)
write.table(res$trait_counts_atgranularity1, file.path(out,"traits","trait_counts_g1.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(res$trait_counts_atgranularity2, file.path(out,"traits","trait_counts_g2.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
if (!is.null(res$genes_detected)) {
  write.table(data.frame(gene=res$genes_detected), file.path(out,"traits","genes_detected.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

# 简单自检打印
g1 <- file.path(out,"traits","trait_counts_g1.tsv")
g2 <- file.path(out,"traits","trait_counts_g2.tsv")
cat("[g1 rows>0] ", if (file.exists(g1)) max(0, nrow(fread(g1)) - 1) else 0, "\n", sep="")
cat("[g2 rows>0] ", if (file.exists(g2)) max(0, nrow(fread(g2)) - 1) else 0, "\n", sep="")
RS

echo "[OK] $sample -> $OUT"
