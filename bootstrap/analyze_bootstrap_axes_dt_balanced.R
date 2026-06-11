#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(vegan)
})

# =========================
# paths
# =========================

boot_dir <- "/root/mcoa_project/mcoa_bootstrap_dt_balanced"

raw_dir  <- "/root/mcoa_project/data_raw"

# =========================
# tools
# =========================

get_rec_cols <- function(X){

  rec_cols <- character(0)

  for(j in seq_along(X)){

    v <- X[[j]]
    v <- v[is.finite(v)]

    if(length(v)==0) next
    if(!is.numeric(v)) next

    if(all(v >= 0) &&
       (max(v) > 20 ||
        (max(v) > 1 &&
         quantile(v,0.9)/max(v) < 0.1))){

      rec_cols <- c(
        rec_cols,
        colnames(X)[j]
      )
    }
  }

  rec_cols
}

clean_matrix <- function(X){

  X <- X[
    ,
    colSums(is.na(X)) < nrow(X),
    drop=FALSE
  ]

  X <- X[
    rowSums(is.na(X)) == 0,
    ,
    drop=FALSE
  ]

  X
}

clean_names <- function(x){

  x <- tolower(x)

  x <- sub(
    "^x(?=[0-9])",
    "",
    x,
    perl=TRUE
  )

  x <- gsub("\\(","",x)
  x <- gsub("\\)","",x)

  x <- gsub(
    "[^a-z0-9]",
    "_",
    x
  )

  x <- gsub("_+","_",x)

  x <- gsub("^_|_$","",x)

  x
}

# =========================
# all bootstrap
# =========================

runs <- list.dirs(
  boot_dir,
  recursive=FALSE,
  full.names=TRUE
)

for(r in runs){

  run_name <- basename(r)

  cat(
    "\n==========",
    run_name,
    "==========\n"
  )

  sample_file <- file.path(
    r,
    "samples.tsv"
  )

  if(!file.exists(sample_file))
    next

  # =========================
  # read coordinates
  # =========================

  S <- fread(sample_file)

  ids <- S$id

  # =========================
  # read matrix
  # =========================

  files <- c(
    "Ko.merged_mean.tsv",
    "Cog.merged_mean.tsv",
    "Cazy.merged_mean.tsv",
    "Traits.merged_mean.tsv"
  )

  data_list <- list()

  for(f in files){

    dt <- fread(
      file.path(
        raw_dir,
        f
      )
    )

    rn <- dt[[1]]

    dt[[1]] <- NULL

    keep <- rn %in% ids

    X <- as.data.frame(
      dt[keep,]
    )

    rownames(X) <- rn[keep]

    X[] <- lapply(
      X,
      function(x)
        suppressWarnings(
          as.numeric(
            as.character(x)
          )
        )
    )

    rec_cols <- get_rec_cols(X)

    for(col in rec_cols){
      X[[col]] <- log1p(
        X[[col]]
      )
    }

    data_list[[f]] <- clean_matrix(X)
  }

  # =========================
  # integrate matrix
  # =========================

  DF <- data.frame(

    cbind(

      data_list[["Ko.merged_mean.tsv"]],

      data_list[["Cog.merged_mean.tsv"]],

      data_list[["Cazy.merged_mean.tsv"]],

      data_list[["Traits.merged_mean.tsv"]]

    )

  )

  common_ids <- intersect(
    rownames(DF),
    S$id
  )

  DF <- DF[
    common_ids,
    ,
    drop=FALSE
  ]

  S <- S[
    match(
      common_ids,
      S$id
    )
  ]

  vars <- colnames(DF)

  # =========================
  # LM
  # =========================

  out <- data.frame(

    var=vars,

    coef1=NA,
    p1=NA,
    R2_1=NA,

    coef2=NA,
    p2=NA,
    R2_2=NA

  )

  for(j in seq_along(vars)){

    y <- DF[[vars[j]]]

    sm1 <- summary(
      lm(
        y ~ S$SynVar1
      )
    )

    out$coef1[j] <- sm1$coefficients[2,1]
    out$p1[j]    <- sm1$coefficients[2,4]
    out$R2_1[j]  <- sm1$r.squared

    sm2 <- summary(
      lm(
        y ~ S$SynVar2
      )
    )

    out$coef2[j] <- sm2$coefficients[2,1]
    out$p2[j]    <- sm2$coefficients[2,4]
    out$R2_2[j]  <- sm2$r.squared
  }

  fwrite(
    out,
    file.path(
      r,
      "variables_vs_axes_lm.tsv"
    ),
    sep="\t"
  )

  # =========================
  # if Tco exits
  # =========================

  tco_file <- file.path(
    r,
    "variables_contribution_Tco.tsv"
  )

  if(file.exists(tco_file)){

    Tco <- fread(
      tco_file
    )

    out$var <- clean_names(
      out$var
    )

    Tco$var <- clean_names(
      Tco$var
    )

    lm_long <- rbindlist(

      list(

        data.table(
          var=out$var,
          axis=1,
          coef=out$coef1,
          pval=out$p1,
          r2=out$R2_1
        ),

        data.table(
          var=out$var,
          axis=2,
          coef=out$coef2,
          pval=out$p2,
          r2=out$R2_2
        )

      )

    )

    tco_long <- melt(

      Tco[,.(var,SV1,SV2)],

      id.vars="var",

      variable.name="axis",

      value.name="SV"

    )

    tco_long[
      ,
      axis := as.integer(
        sub(
          "SV",
          "",
          axis
        )
      )
    ]

    df_merge <- merge(

      lm_long,

      tco_long,

      by=c(
        "var",
        "axis"
      )

    )

    for(ax in 1:2){

      sub <- df_merge[

        axis==ax &
        pval < 1e-3 &
        r2 > 0.2

      ]

      sub <- sub[
        order(
          -abs(SV)
        )
      ]

      fwrite(

        head(
          sub,
          50
        ),

        file.path(
          r,
          paste0(
            "axis",
            ax,
            ".tsv"
          )
        ),

        sep="\t"

      )
    }
  }

  cat(
    "DONE:",
    run_name,
    "\n"
  )
}

cat(
  "\nALL RUNS COMPLETE\n"
)
