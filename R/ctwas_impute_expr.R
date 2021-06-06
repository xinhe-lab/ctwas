#' Impute expression
#' @param weight a string, pointing to the fusion/twas format of weights.or predictdb format.
#'   Note the effect size are obtained with standardized genotype and phenotype.
#' @param method a string,  blup/bslmm/lasso/top1/enet/best
#'   "best" means the method giving the best cross validation R^2, only used for fusion weight format.
#' @param harmonize T/F, if need to harmonize SNP data, if T, will harmonize GWAS and eQTL genotype alleles.
#' @importFrom logging addHandler loginfo
#'
#' @export
impute_expr <- function(pgenfs,
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

  exprfs <- vector()

  for (b in 1:22){
    pgenf <- pgenfs[b]
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

    loginfo("Reading weights for chromosome ", b)
    if (dir.exists(weight)){
      weightall <- read_weight_fusion(weight, b, snpinfo, z_snp = NULL, method = method, harmonize = T)
    } else if (file_ext(weight)=='db'){
      weightall <- read_weight_predictdb(weight, b, snpinfo, z_snp = NULL, harmonize = T)
    } else{
      stop("Unrecognized weight format, need to use either FUSION format or predict.db format")
    }

    exprlist <- weightall[["exprlist"]]
    qclist <- weightall[["qclist"]]

    if (length(exprlist) > 0) {
      loginfo("Start gene z score imputation ...")
      for (gname in names(exprlist)){
        wgt <- exprlist[[gname]][["wgt"]]
        snpnames <- rownames(wgt)
        gwas.idx <-  match(snpnames, snpinfo$id)
        g <- read_pgen(pgen, variantidx = gwas.idx)
        g <- scale(g)  # genotypes are standardized
        gexpr <- g %*% wgt
        exprlist[[gname]][["expr"]] <- gexpr
      }
    }

    loginfo ("Imputation done, writing results to output...")
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

    exprfs[b] <- exprf
  }

  return(exprfs)
}

