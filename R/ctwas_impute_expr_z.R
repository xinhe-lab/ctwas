#' Impute expression
#' @param zdf a data frame, with columns "id" and "z". Z scores for every SNP.
#' @param weight a string, pointing to the fusion/twas format of weights.
#'   Note the effect size are obtained with standardized genotype and phenotype.
#' @param method a string,  blup/bslmm/lasso/top1/enet/best
#'   "best" means the method giving the best cross validation R^2
#' @param checksnps T/F, if need to check SNP consistency between weights
#'    file and genotype. fusion format of weights gives snp information from plink .bim file. checking include chr, pos, A1, A2.
#' @importFrom logging addHandler loginfo
#'
#' @export
impute_expr_z <- function(zdf, ld_pgenf,
                        weight,
                        method = "lasso",
                        outputdir = getwd(),
                        outname = NULL,
                        logfile = NULL,
                        compress = T,
                        checksnps = F){

  if (!is.null(logfile)){
    addHandler(writeToFile, file= logfile, level='DEBUG')
  }

  outname <- file.path(outputdir, outname)

  ld_pvarf <- prep_pvar(ld_pgenf, outputdir = outputdir)
  ld_snpinfo <- read_pvar(ld_pvarf)

  ld_pgen <- prep_pgen(ld_pgenf, ld_pvarf)

  b <- unique(ld_snpinfo$chrom)
  if (length(b) !=1){
    stop("Input LD reference genotype not splited by chromosome")
  }

  exprf <- paste0(outname, "_chr", b, ".expr")
  exprvarf <- paste0(outname, "_chr", b, ".exprvar")
  exprqcf <-  paste0(outname, "_chr", b, ".exprqc.Rd")

  loginfo('expression z score inmputation started for chr %s.', b)

  exprlist <- list()
  qclist <- list()
  wgtdir <- dirname(weight)
  wgtposfile <- file.path(wgtdir, paste0(basename(weight), ".pos"))

  wgtpos <- read.table(wgtposfile, header = T, stringsAsFactors = F)
  wgtpos <- transform(wgtpos,
                      "ID" = ifelse(duplicated(ID) | duplicated(ID, fromLast=TRUE),
                          paste(ID, ave(ID, ID, FUN=seq_along), sep='_ID'), ID))

  loginfo("number of genes with weights provided: %s", nrow(wgtpos))

  for (i in 1: nrow(wgtpos)){
    wf <- file.path(wgtdir, wgtpos[i, "WGT"])
    load(wf)
    gname <- wgtpos[i, "ID"]

    g.method = method
    if (g.method == "best"){
      g.method = names(which.max(cv.performance["rsq", ]))
    }
    if (!(g.method %in% names(cv.performance[1,]))) next

    wgt.matrix <- wgt.matrix[abs(wgt.matrix[, g.method]) > 0, , drop = F]
    wgt.matrix <- wgt.matrix[complete.cases(wgt.matrix), ,drop = F]

    if (nrow(wgt.matrix) == 0) next

    snpnames <- Reduce(intersect,
                       list(rownames(wgt.matrix), ld_snpinfo$id, zdf$id))

    if (length(snpnames) == 0) next

    wgt.idx <- match(snpnames, rownames(wgt.matrix))
    ld.idx <-  match(snpnames, ld_snpinfo$id)
    zdf.idx <- match(snpnames, zdf$id)

    if (checksnps){
      # `snps` from FUSION follows .bim format
      snps$V3 <- NULL
      colnames(snps) <- c("chrom", "id", "pos", "alt", "ref")
      ldsnps <- ld_snpinfo[ld.idx, ]

      if (!identical(ldsnps, snps)){
        stop("LD reference SNP and eQTL info inconsistent. STOP.")
      }
    }

    wgt <-  wgt.matrix[wgt.idx, g.method, drop = F]
    X.g <- read_pgen(ld_pgen, variantidx = ld.idx)

    # genotypes are standardized, as weights are standardized
    X.g <- scale(X.g)

    gexpr <- X.g %*% wgt
    if (abs(max(gexpr) - min(gexpr)) < 1e-8) next

    z.s <- as.matrix(zdf[zdf.idx, "z"])
    var.s <- sqrt(apply(X.g, 2, var))
    Gamma.g <- cov(X.g)

    z.g <-  (t(wgt) * var.s) %*% z.s/ sqrt(t(wgt) %*% Gamma.g %*% wgt)

    exprlist[[gname]][["expr"]] <- gexpr
    exprlist[[gname]][["z.g"]] <- z.g
    exprlist[[gname]][["chrom"]] <- wgtpos[wgtpos$ID == gname, "CHR"]
    exprlist[[gname]][["p0"]] <- wgtpos[wgtpos$ID == gname, "P0"]
    exprlist[[gname]][["p1"]] <- wgtpos[wgtpos$ID == gname, "P1"]
    exprlist[[gname]][["wgt"]] <- wgt

    qclist[[gname]][["n"]] <- nrow(wgt.matrix)
    qclist[[gname]][["nmiss"]] <-  nrow(wgt.matrix) - length(snpnames)
    qclist[[gname]][["missrate"]] <-
      qclist[[gname]][["nmiss"]]/qclist[[gname]][["n"]]

  }

  z.g <- unlist(lapply(exprlist,'[[', "z.g"))
  expr <- do.call(cbind, lapply(exprlist, '[[', "expr"))
  gnames <- names(exprlist)
  chrom <- unlist(lapply(exprlist,'[[', "chrom"))
  p0 <- unlist(lapply(exprlist,'[[', "p0"))
  p1 <- unlist(lapply(exprlist,'[[', "p1"))
  wgtlist <- lapply(exprlist, '[[', "wgt")

  if (length(exprlist) == 0){
    expr <- data.table::data.table(NULL)
    geneinfo <- data.table::data.table(NULL)
  }

  loginfo("Number of genes with imputed expression: %s for chr %s",
          ncol(expr), b)

  data.table::fwrite(expr, file = exprf,
                     row.names = F, col.names = F,
                     sep = "\t", quote = F)

  if (isTRUE(compress)){
    system(paste0("gzip -f ", exprf))
    exprf <- paste0(exprf, '.gz')
  }

  geneinfo <- data.frame("chrom" = chrom,
                         "id" = gnames,
                         "p0" = p0,
                         "p1" = p1)
  data.table::fwrite(geneinfo, file = exprvarf, sep = "\t", quote = F)

  save(wgtlist, qclist, file = exprqcf)

  loginfo('expression inmputation done for chr %s.', b)

  zdf.g <- data.frame("id" = gnames,
                       "z" = z.g)

  return(list("zdf" = zdf.g, "ld_exprf"= exprf))
}

