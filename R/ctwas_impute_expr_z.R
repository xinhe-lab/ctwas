#' Impute expression
#' @param z_snp A data frame with two columns: "id", "A1", "A2", "z". giving the z scores for
#' snps. "A1" is effect allele. "A2" is the other allele. If `harmonize= False`, A1 and A2 are not required.
#' @param weight a string, pointing to the fusion/twas format of weights.
#'   Note the effect size are obtained with standardized genotype and phenotype.
#' @param method a string,  blup/bslmm/lasso/top1/enet/best
#'   "best" means the method giving the best cross validation R^2
#' @param harmonize T/F, if need to harmonize SNP data, if T, will harmonize gwas, LD and eQTL genotype alleles.
#' @importFrom logging addHandler loginfo
#'
#' @export
impute_expr_z <- function(z_snp,
                        ld_pgenf,
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

  outname <- file.path(outputdir, outname)

  ld_pvarf <- prep_pvar(ld_pgenf, outputdir = outputdir)
  ld_snpinfo <- read_pvar(ld_pvarf)

  ld_pgen <- prep_pgen(ld_pgenf, ld_pvarf)

  b <- unique(ld_snpinfo$chrom)
  if (length(b) !=1){
    stop("Input LD reference genotype not splited by chromosome")
  }

  if (isTRUE(harmonize)){
    loginfo("flipping z scores to match LD reference")
    z_snp <- harmonize_z_ld(z_snp, ld_snpinfo)
    loginfo("will also flip weights to match LD reference for each gene")
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

    if (isTRUE(harmonize)){
      w <- harmonize_wgt_ld(wgt.matrix, snps, ld_snpinfo)
      wgt.matrix <- w[["wgt"]]
      snps <- w[["snps"]]
    }

    snpnames <- Reduce(intersect,
                       list(rownames(wgt.matrix), ld_snpinfo$id, zdf$id))

    if (length(snpnames) == 0) next

    wgt.idx <- match(snpnames, rownames(wgt.matrix))
    ld.idx <-  match(snpnames, ld_snpinfo$id)
    z.idx <- match(snpnames, z_snp$id)

    wgt <-  wgt.matrix[wgt.idx, g.method, drop = F]
    X.g <- read_pgen(ld_pgen, variantidx = ld.idx)

    # genotypes are standardized, as weights are standardized
    X.g <- scale(X.g)

    gexpr <- X.g %*% wgt
    if (abs(max(gexpr) - min(gexpr)) < 1e-8) next

    z.s <- as.matrix(z_snp[z.idx, "z"])
    var.s <- sqrt(apply(X.g, 2, var))
    Gamma.g <- cov(X.g)

    z.g <-  (t(wgt) * var.s) %*% z.s/ sqrt(t(wgt) %*% Gamma.g %*% wgt)

    exprlist[[gname]][["expr"]] <- gexpr
    exprlist[[gname]][["z.g"]] <- z.g
    exprlist[[gname]][["chrom"]] <- wgtpos[wgtpos$ID == gname, "CHR"]
    exprlist[[gname]][["p0"]] <-  min(snps[snps[, "id"] %in% snpnames, "pos"])  # eQTL positions
    exprlist[[gname]][["p1"]] <-  max(snps[snps[, "id"] %in% snpnames, "pos"])  # eQTL positions
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

  z_gene <- data.frame("id" = gnames,
                       "z" = z.g)

  return(list("z_gene" = z_gene, "ld_exprf"= exprf))
}

