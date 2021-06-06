#' Impute expression
#' @param weight a string, pointing to the fusion/twas format of weights.or predictdb format.
#'   Note the effect size are obtained with standardized genotype and phenotype.
#' @param method a string,  blup/bslmm/lasso/top1/enet/best
#'   "best" means the method giving the best cross validation R^2, only used for fusion weight format.
#' @param harmonize T/F, if need to harmonize SNP data, if T, will harmonize GWAS and eQTL genotype alleles.
#' @importFrom logging addHandler loginfo
#'
#' @export
impute_expr <- function(pgenf,
                        weight,
                        method = "lasso",
                        outputdir = getwd(),
                        outname = NULL,
                        logfile = NULL,
                        compress = T,
                        harmonize = T){

  if (!is.null(logfile)){
    addHandler(writeToFile, file= logfile, level='DEBUG')
  }

  if (isTRUE(harmonize)){
    loginfo("Harmonize set to True: will flip weights to match LD reference")
  }

  outname <- file.path(outputdir, outname)

  pvarf <- prep_pvar(pgenf, outputdir = outputdir)
  snpinfo <- read_pvar(pvarf)

  pgen <- prep_pgen(pgenf, pvarf)

  b <- unique(snpinfo$chrom)
  if (length(b) !=1){
    stop("Input genotype not splited by chromosome")
  }

  exprf <- paste0(outname, "_chr", b, ".expr")
  exprvarf <- paste0(outname, "_chr", b, ".exprvar")
  exprqcf <-  paste0(outname, "_chr", b, ".exprqc.Rd")

  loginfo('expression inmputation started for chr %s.', b)

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

    if (isTRUE(harmonize)){
      w <- harmonize_wgt_ld(wgt.matrix, snps, snpinfo)
      wgt.matrix <- w[["wgt"]]
      snps <- w[["snps"]]
    }

    g.method = method
    if (g.method == "best"){
      g.method = names(which.max(cv.performance["rsq", ]))
    }
    if (!(g.method %in% names(cv.performance[1,]))) next

    wgt.matrix <- wgt.matrix[abs(wgt.matrix[, g.method]) > 0, , drop = F]
    wgt.matrix <- wgt.matrix[complete.cases(wgt.matrix), ,drop = F]

    if (nrow(wgt.matrix) == 0) next

    snpnames <- intersect(rownames(wgt.matrix), snpinfo$id)
    if (length(snpnames) == 0) next

    wgt.idx <- match(snpnames, rownames(wgt.matrix))
    gwas.idx <-  match(snpnames, snpinfo$id)

    wgt <-  wgt.matrix[wgt.idx, g.method, drop = F]

    g <- read_pgen(pgen, variantidx = gwas.idx)

    # genotypes are standardized
    g <- scale(g)

    gexpr <- g %*% wgt

    if (abs(max(gexpr) - min(gexpr)) < 1e-8) next

    exprlist[[gname]][["expr"]] <- gexpr
    exprlist[[gname]][["chrom"]] <- wgtpos[wgtpos$ID == gname, "CHR"]
    exprlist[[gname]][["p0"]] <-  min(snps[snps[, "id"] %in% snpnames, "pos"])  # eQTL positions
    exprlist[[gname]][["p1"]] <-  max(snps[snps[, "id"] %in% snpnames, "pos"])  # eQTL positions
    exprlist[[gname]][["wgt"]] <- wgt

    qclist[[gname]][["n"]] <- nrow(wgt.matrix)
    qclist[[gname]][["nmiss"]] <-  nrow(wgt.matrix) - length(snpnames)
    qclist[[gname]][["missrate"]] <-
        qclist[[gname]][["nmiss"]]/qclist[[gname]][["n"]]

  }

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

  return(exprf)
}

