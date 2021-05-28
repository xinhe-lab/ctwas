#' Impute expression
#' @param z_snp A data frame with two columns: "id", "A1", "A2", "z". giving the z scores for
#' snps. "A1" is effect allele. "A2" is the other allele. If `harmonize= False`, A1 and A2 are not required.
# @param ld_Rf a file listing all LD matrix R files in one chromosome. It gives the file path for the .RDS file which contains LD matrix for a region/block, one file per line. All .RDS file from the target chromosome should be included here. For each RDS file, a file with same base name but ended with .Rvar needs to be present in the same folder. the .Rvar file has 5 columns: "chrom", "id", "pos", "alt", "ref". The order of rows needs to match the order of rows in .RDS file.
#' @param weight a string, pointing to the fusion/twas format of weights.
#'   Note the effect size are obtained with standardized genotype and phenotype.
#' @param method a string,  blup/bslmm/lasso/top1/enet/best
#'   "best" means the method giving the best cross validation R^2
#' @param harmonize T/F, if need to harmonize SNP data, if T, will harmonize gwas, LD and eQTL genotype alleles.
#' @importFrom logging addHandler loginfo
#'
#' @export

impute_expr_z <- function (z_snp,
                           ld_pgenf = NULL,
                           ld_Rf = NULL,
                           weight,
                           method = "lasso",
                           outputdir = getwd(),
                           outname = NULL,
                           logfile = NULL,
                           compress = T,
                           harmonize = T){


  dir.create(outputdir, showWarnings=F)

  if (!is.null(logfile)) {
    addHandler(writeToFile, file = logfile, level = "DEBUG")
  }
  outname <- file.path(outputdir, outname)

  if (is.null(ld_pgenfs) & is.null(ld_Rf)){
    stop("Stopped: missing LD information.
         LD information needs to be provided either in genotype form
         (see parameter description for ld_pgenf) or R matrix form
         (see parameter description for ld_Rf) ")
  }

  if (!is.null(ld_pgenf)){
    ld_pvarf <- prep_pvar(ld_pgenf, outputdir = outputdir)
    ld_snpinfo <- read_pvar(ld_pvarf)
  } else {
    ld_Rinfo <- read_ld_Rinfo(ld_Rf)
    ld_snpinfo <- read_ld_Rvar(ld_Rf)
  }

  b <- unique(ld_snpinfo$chrom)
  if (length(b) != 1) {
    stop("Input LD reference not split by chromosome")
  }

  if (isTRUE(harmonize)) {
    logging::loginfo("flipping z scores to match LD reference")
    z_snp <- harmonize_z_ld(z_snp, ld_snpinfo)
    logging::loginfo("will also flip weights to match LD reference for each gene")
  }

  exprlist <- list()
  qclist <- list()
  wgtdir <- dirname(weight)
  wgtposfile <- file.path(wgtdir, paste0(basename(weight), ".pos"))

  wgtpos <- read.table(wgtposfile, header = T, stringsAsFactors = F)
  wgtpos <- transform(wgtpos,
                      ID = ifelse(duplicated(ID) | duplicated(ID, fromLast = TRUE),
                                  paste(ID, ave(ID, ID, FUN = seq_along), sep = "_ID"), ID))
  loginfo("number of genes with weights provided: %s", nrow(wgtpos))

  wgtpos <- wgtpos[wgtpos$CHR==b,]
  loginfo("number of genes on chromosome %s: %s", b, nrow(wgtpos))

  loginfo("collecting gene weight information ...")
  if (nrow(wgtpos) > 0){
    for (i in 1:nrow(wgtpos)) {
      wf <- file.path(wgtdir, wgtpos[i, "WGT"])
      load(wf)
      gname <- wgtpos[i, "ID"]
      if (isTRUE(harmonize)) {
        w <- harmonize_wgt_ld(wgt.matrix, snps, ld_snpinfo)
        wgt.matrix <- w[["wgt"]]
        snps <- w[["snps"]]
      }
      g.method = method
      if (g.method == "best") {
        g.method = names(which.max(cv.performance["rsq",]))
      }
      if (!(g.method %in% names(cv.performance[1, ])))
        next
      wgt.matrix <- wgt.matrix[abs(wgt.matrix[, g.method]) > 0, , drop = F]
      wgt.matrix <- wgt.matrix[complete.cases(wgt.matrix), , drop = F]
      if (nrow(wgt.matrix) == 0)
        next
      snpnames <- Reduce(intersect, list(rownames(wgt.matrix), ld_snpinfo$id, z_snp$id))
      if (length(snpnames) == 0)
        next
      wgt.idx <- match(snpnames, rownames(wgt.matrix))
      wgt <- wgt.matrix[wgt.idx, g.method, drop = F]

      p0 <-  min(snps[snps[, "id"] %in% snpnames, "pos"])
      p1 <- max(snps[snps[, "id"] %in% snpnames, "pos"])

      exprlist[[gname]] <- list("chrom" = b,
                                "p0" = p0,
                                "p1" = p1,
                                "wgt" = wgt)

      if (is.null(ld_pgenf)){
        ifreg <- ifelse(p1 >= ld_Rinfo[, "start"] & p0 < ld_Rinfo[, "stop"], T, F)
        exprlist[[gname]][["reg"]] <- paste(sort(ld_Rinfo[ifreg, "region_name"]),
                                            collapse = ";")
      }

      nwgt <- nrow(wgt.matrix)
      nmiss <- nrow(wgt.matrix) - length(snpnames)
      qclist[[gname]] <- list("n" = nwgt,
                              "nmiss" = nmiss,
                              "missrate" = nwgt/nmiss)

    }

    loginfo("Start gene z score imputation ...")
    if (!is.null(ld_pgenf)){
      loginfo("ld genotype is given, using genotypes to impute gene z score.")
      ld_pgen <- prep_pgen(ld_pgenf, ld_pvarf)
      gnames <- names(exprlist)
      for (i in 1:length(gnames)){
        gname <- gnames[i]
        wgt <- exprlist[[gname]][["wgt"]]
        snpnames <- rownames(wgt)
        ld.idx <- match(snpnames, ld_snpinfo$id)
        z.idx <- match(snpnames, z_snp$id)
        X.g <- read_pgen(ld_pgen, variantidx = ld.idx)
        X.g <- scale(X.g)
        gexpr <- X.g %*% wgt
        if (abs(max(gexpr) - min(gexpr)) < 1e-08) {
          exprlist[[gname]] <- NULL
        } else {
          z.s <- as.matrix(z_snp[z.idx, "z"])
          var.s <- sqrt(apply(X.g, 2, var))
          Gamma.g <- cov(X.g)
          z.g <- (t(wgt) * var.s) %*% z.s/sqrt(t(wgt) %*% Gamma.g %*% wgt)
          exprlist[[gname]][["expr"]] <- gexpr
          exprlist[[gname]][["z.g"]] <- z.g
        }
      }
    } else {
      loginfo("Using given LD matrices to impute gene z score.")
      regs <- data.frame("gid" = names(exprlist),
                         "reg" = unlist(lapply(exprlist, "[[", "reg")))
      batches <- unique(regs$reg)
      for (batch in batches){
        gnames <- regs[regs$reg == batch, "gid"]
        regnames <- strsplit(batch, ";")[[1]]
        regRDS <- ld_Rinfo[match(regnames, ld_Rinfo$region_name), "RDS_file"]
        R_snp <- lapply(regRDS, readRDS)
        R_snp <- as.matrix(Matrix::bdiag(R_snp))
        R_snp_anno <- do.call(rbind, lapply(regRDS, read_ld_Rvar_RDS))
        for (i in 1:length(gnames)){
          wgt <- exprlist[[gname]][["wgt"]]
          snpnames <- rownames(wgt)
          ld.idx <- match(snpnames, R_snp_anno$id)
          zdf.idx <- match(snpnames, z_snp$id)
          R.s <- R_snp[ld.idx,ld.idx]
          z.s <-  as.matrix(z_snp[zdf.idx, "z"])
          z.g <- crossprod(wgt,z.s)/sqrt(crossprod(wgt,R.s)%*%wgt)
          exprlist[[gname]][["z.g"]] <- z.g
        }
      }
    }
  }


  loginfo ("Imputation done, writing results to output...")
  z.g <- unlist(lapply(exprlist,'[[', "z.g"))
  gnames <- names(exprlist)
  chrom <- unlist(lapply(exprlist,'[[', "chrom"))
  p0 <- unlist(lapply(exprlist,'[[', "p0"))
  p1 <- unlist(lapply(exprlist,'[[', "p1"))
  wgtlist <- lapply(exprlist, '[[', "wgt")

  exprvarf <- paste0(outname, "_chr", b, ".exprvar")
  if (length(exprlist) == 0){
    geneinfo <- data.table::data.table(NULL)
  } else {
    geneinfo <- data.frame("chrom" = chrom,
                           "id" = gnames,
                           "p0" = p0,
                           "p1" = p1)
  }
  data.table::fwrite(geneinfo, file = exprvarf, sep = "\t", quote = F)

  exprqcf <- paste0(outname, "_chr", b, ".exprqc.Rd")
  save(wgtlist, qclist, file = exprqcf)

  exprf <- paste0(outname, "_chr", b, ".expr")
  if (!is.null(ld_pgenf)){
    if (length(exprlist) == 0){
      expr <- data.table::data.table(NULL)
    } else {
      expr <- do.call(cbind, lapply(exprlist, '[[', "expr"))
    }
  } else {
    # an empty .expr file when using LD R matrices
    expr <- data.table::data.table(NA)
  }

  data.table::fwrite(expr, file = exprf,
                     row.names = F, col.names = F,
                     sep = "\t", quote = F)
  if (isTRUE(compress)){
    system(paste0("gzip -f ", exprf))
    exprf <- paste0(exprf, '.gz')
  }

  loginfo("Imputation done: number of genes with imputed expression: %s for chr %s",
          length(gnames), b)

  z_gene <- data.frame("id" = gnames,
                       "z" = z.g)

  return(list("z_gene" = z_gene, "ld_exprf"= exprf))
}

