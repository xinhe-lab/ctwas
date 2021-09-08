#' Impute expression
#' @param z_snp A data frame with two columns: "id", "A1", "A2", "z". giving the z scores for
#' snps. "A1" is effect allele. "A2" is the other allele. If `harmonize= False`, A1 and A2 are not required.
# @param ld_Rf a file listing all LD matrix R files in one chromosome. It gives the file path for the .RDS file which contains LD matrix for a region/block, one file per line. All .RDS file from the target chromosome should be included here. For each RDS file, a file with same base name but ended with .Rvar needs to be present in the same folder. the .Rvar file has 5 columns: "chrom", "id", "pos", "alt", "ref". The order of rows needs to match the order of rows in .RDS file.
#' @param weight a string, pointing to the fusion/twas format of weights, or predictdb format.
#'   Note the effect size are obtained with standardized genotype and phenotype.
#'
#' @param method a string,  blup/bslmm/lasso/top1/enet/best
#'   "best" means the method giving the best cross validation R^2, this is only used for fusion weights.
#' @param harmonize T/F, if need to harmonize SNP data, if T, will harmonize gwas, LD and eQTL genotype alleles.
#' @importFrom logging addHandler loginfo
#' @importFrom tools file_ext
#'
#' @export
impute_expr_z <- function (z_snp, weight, ld_pgenfs = NULL, ld_R_dir = NULL, method = "lasso", 
                           outputdir = getwd(), outname = NULL, logfile = NULL, compress = T, 
                           harmonize_z = T, harmonize_wgt = T, 
                           strand_ambig_action_z = c("drop", "none", "recover"), recover_strand_ambig_wgt = T){
  dir.create(outputdir, showWarnings = F)
  if (!is.null(logfile)) {
    addHandler(writeToFile, file = logfile, level = "DEBUG")
  }
  if (is.null(ld_pgenfs) & is.null(ld_R_dir)) {
    stop("Stopped: missing LD information.\n         LD information needs to be provided either in genotype form\n         (see parameter description for ld_pgenfs) or R matrix form\n         (see parameter description for ld_R_dir) ")
  } else if (is.null(ld_pgenfs)) {
    ld_Rfs <- write_ld_Rf(ld_R_dir, outname = outname, outputdir = outputdir)
  }
  outname <- file.path(outputdir, outname)
  ld_exprfs <- vector()
  z_genelist <- list()
  ld_snplist <- c() #list to store names of snps in ld reference
  for (b in 1:22) {
    if (!is.null(ld_pgenfs)) {
      ld_pgenf <- ld_pgenfs[b]
      ld_pvarf <- prep_pvar(ld_pgenf, outputdir = outputdir)
      ld_snpinfo <- read_pvar(ld_pvarf)
    } else {
      ld_Rf <- ld_Rfs[b]
      ld_Rinfo <- data.table::fread(ld_Rf, header = T)
      ld_snpinfo <- read_ld_Rvar(ld_Rf)
    }
    chrom <- unique(ld_snpinfo$chrom)
    if (length(chrom) != 1) {
      stop("Input LD reference not split by chromosome")
    }
    ld_snplist <- c(ld_snplist, ld_snpinfo$id) #store names of snps in ld reference
    if (isTRUE(harmonize_z)) {
      loginfo("Flipping z scores to match LD reference")
      z_snp <- harmonize_z_ld(z_snp, ld_snpinfo,
                              strand_ambig_action = strand_ambig_action_z, 
                              ld_pgenfs = ld_pgenfs, 
                              ld_Rinfo = ld_Rinfo)
    }
    loginfo("Reading weights for chromosome %s", b)
    if (dir.exists(weight)) {
      weightall <- read_weight_fusion(weight, b, ld_snpinfo, z_snp, method = method, 
                                      harmonize_wgt=harmonize_wgt)
    } else if (file_ext(weight) == "db") {
      weightall <- read_weight_predictdb(weight, b, ld_snpinfo, z_snp, 
                                         harmonize_wgt=harmonize_wgt, ld_Rinfo=ld_Rinfo,
                                         recover_strand_ambig=recover_strand_ambig_wgt)
    } else {
      stop("Unrecognized weight format, need to use either FUSION format or predict.db format")
    }
    exprlist <- weightall[["exprlist"]]
    qclist <- weightall[["qclist"]]
    if (length(exprlist) > 0) {
      loginfo("Start gene z score imputation ...")
      if (!is.null(ld_pgenfs)) {
        loginfo("ld genotype is given, using genotypes to impute gene z score.")
        ld_pgen <- prep_pgen(ld_pgenf, ld_pvarf)
        gnames <- names(exprlist)
        for (i in 1:length(gnames)) {
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
            z.g <- (t(wgt) * var.s) %*% z.s/sqrt(t(wgt) %*% 
                                                   Gamma.g %*% wgt)
            exprlist[[gname]][["expr"]] <- gexpr
            exprlist[[gname]][["z.g"]] <- z.g
          }
        }
      } else {
        loginfo("Using given LD matrices to impute gene z score.")
        for (gname in names(exprlist)) {
          p0 <- exprlist[[gname]][["p0"]]
          p1 <- exprlist[[gname]][["p1"]]
          ifreg <- ifelse(p1 >= ld_Rinfo[, "start"] & 
                            p0 < ld_Rinfo[, "stop"], T, F)
          exprlist[[gname]][["reg"]] <- paste(sort(ld_Rinfo[ifreg, 
                                                            "region_name"]), collapse = ";")
        }
        regs <- data.frame(gid = names(exprlist), reg = unlist(lapply(exprlist, 
                                                                      "[[", "reg")), stringsAsFactors = F)
        batches <- unique(regs$reg)
        for (batch in batches) {
          gnames <- regs[regs$reg == batch, "gid"]
          regnames <- strsplit(batch, ";")[[1]]
          regRDS <- ld_Rinfo[match(regnames, ld_Rinfo$region_name), 
                             "RDS_file"]
          R_snp <- lapply(regRDS, readRDS)
          R_snp <- as.matrix(Matrix::bdiag(R_snp))
          R_snp_anno <- do.call(rbind, lapply(regRDS, 
                                              read_ld_Rvar_RDS))
          for (i in 1:length(gnames)) {
            gname <- gnames[i]
            wgt <- exprlist[[gname]][["wgt"]]
            snpnames <- rownames(wgt)
            ld.idx <- match(snpnames, R_snp_anno$id)
            zdf.idx <- match(snpnames, z_snp$id)
            R.s <- R_snp[ld.idx, ld.idx]
            z.s <- as.matrix(z_snp[zdf.idx, "z"])
            z.g <- crossprod(wgt, z.s)/sqrt(crossprod(wgt, 
                                                      R.s) %*% wgt)
            exprlist[[gname]][["z.g"]] <- z.g
          }
        }
      }
    }
    loginfo("Imputation done, writing results to output...")
    z.g <- unlist(lapply(exprlist, "[[", "z.g"))
    gnames <- names(exprlist)
    chrom <- unlist(lapply(exprlist, "[[", "chrom"))
    p0 <- unlist(lapply(exprlist, "[[", "p0"))
    p1 <- unlist(lapply(exprlist, "[[", "p1"))
    wgtlist <- lapply(exprlist, "[[", "wgt")
    exprvarf <- paste0(outname, "_chr", b, ".exprvar")
    if (length(exprlist) == 0) {
      geneinfo <- data.table::data.table(NULL)
    } else {
      geneinfo <- data.frame(chrom = chrom, id = gnames, 
                             p0 = p0, p1 = p1)
    }
    data.table::fwrite(geneinfo, file = exprvarf, sep = "\t", 
                       quote = F)
    z_gene_chr <- data.frame(id = gnames, z = z.g)
    exprqcf <- paste0(outname, "_chr", b, ".exprqc.Rd")
    save(wgtlist, qclist, z_gene_chr, file = exprqcf)
    exprf <- paste0(outname, "_chr", b, ".expr")
    if (!is.null(ld_pgenfs)) {
      if (length(exprlist) == 0) {
        expr <- data.table::data.table(NULL)
      } else {
        expr <- do.call(cbind, lapply(exprlist, "[[", 
                                      "expr"))
      }
    } else {
      expr <- data.table::data.table(NA)
    }
    data.table::fwrite(expr, file = exprf, row.names = F, 
                       col.names = F, sep = "\t", quote = F)
    if (isTRUE(compress)) {
      system(paste0("gzip -f ", exprf))
      exprf <- paste0(exprf, ".gz")
    }
    loginfo("Imputation done: number of genes with imputed expression: %s for chr %s", 
            length(gnames), b)
    z_genelist[[b]] <- z_gene_chr
    ld_exprfs[b] <- exprf
  }
  z_gene <- do.call(rbind, z_genelist)
  z_snp <- z_snp[z_snp$id %in% ld_snplist,] #subset z_snp to snps in ld reference
  return(list(z_gene = z_gene, ld_exprfs = ld_exprfs, z_snp = z_snp))
}
