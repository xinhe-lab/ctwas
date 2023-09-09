#' Impute expression z scores
#' 
#' @param z_snp A data frame with two columns: "id", "A1", "A2", "z". giving the z scores for
#' snps. "A1" is effect allele. "A2" is the other allele. If `harmonize= False`, A1 and A2 are not required.
#' 
#' @param weight a string, pointing to a directory with the fusion/twas format of weights, or a .db file in predictdb format.
#' A vector of multiple sets of weights in PredictDB format can also be specified; genes will have their filename appended
#' to their gene name to ensure IDs are unique.
#' 
#' @param ld_pgenfs a character vector of .pgen or .bed files. One file for one
#' chromosome, in the order of 1 to 22. Therefore, the length of this vector
#' needs to be 22. If .pgen files are given, then .pvar and .psam are assumed
#' to present in the same directory. If .bed files are given, then .bim and
#' .fam files are assumed to present in the same directory.
#'  
#' @param LD_R_dir a string, pointing to a directory containing all LD matrix files and variant information. Expects .RDS files which contain LD correlation matrices for a region/block.
#' For each RDS file, a file with same base name but ended with .Rvar needs to be present in the same folder. the .Rvar file has 5 required columns: "chrom", "id", "pos", "alt", "ref". 
#' If using PredictDB format weights and \code{scale_by_ld_variance=T}, a 6th column is also required: "variance", which is the variance of the each SNP.
#' The order of rows needs to match the order of rows in .RDS file.
#' 
#' @param method a string, blup/bslmm/lasso/top1/enet/best. This option is only used for fusion weights. 
#' "best" means the method giving the best cross #' validation R^2. Note that top1 uses only the weight 
#' with largest effect. 
#'   
#' @param outputdir a string, the directory to store output
#' 
#' @param outname a string, the output name
#' 
#' @param logfile the log file, if NULL will print log info on screen
#' 
#' @param compress TRUE/FALSE. If TRUE, the imputed expression files are compressed
#' 
#' @param harmonize_z TRUE/FALSE. If TRUE, GWAS and eQTL genotype alleles are harmonized
#' 
#' @param harmonize_wgt TRUE/FALSE. If TRUE, GWAS and eQTL genotype alleles are harmonized
#' 
#' @param strand_ambig_action_z the action to take to harmonize strand ambiguous variants (A/T, G/C) between 
#' the z scores and LD reference. "drop" removes the ambiguous variant from the z scores. "none" treats the variant 
#' as unambiguous, flipping the z score to match the LD reference and then taking no additional action. "recover" 
#' imputes the sign of ambiguous z scores using unambiguous z scores and the LD reference and flips the z scores 
#' if there is a mismatch between the imputed sign and the observed sign of the z score. This option is computationally intensive
#' 
#' @param recover_strand_ambig_wgt TRUE/FALSE. If TRUE, a procedure is used to recover strand ambiguous variants. If FALSE, 
#' these variants are dropped from the prediction model. This procedure compares correlations between variants in the the 
#' LD reference and prediction models, and it can only be used with PredictDB format prediction models, which include this
#' information.
#' 
#' @param ncore The number of cores used to parallelize imputation over weights
#' 
#' @param chrom a numeric vector of chromosomes to perform z score imputation over. Useful for large jobs requiring batches
#' 
#' @param scale_by_ld_variance TRUE/FALSE. If TRUE, PredictDB weights are scaled by genotype variance, which is the default 
#' behavior for PredictDB
#' 
#' @importFrom logging addHandler loginfo
#' @importFrom tools file_ext
#'
#' @export
impute_expr_z <- function (z_snp, weight, ld_pgenfs = NULL, ld_R_dir = NULL, method = "lasso", 
                           outputdir = getwd(), outname = NULL, logfile = NULL, compress = T, 
                           harmonize_z = T, harmonize_wgt = T, 
                           strand_ambig_action_z = c("drop", "none", "recover"), recover_strand_ambig_wgt = T,
                           ncore=1, chrom=1:22,
                           scale_by_ld_variance=T){
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
  for (b in chrom) {
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
    if (length(chrom) > 1) {
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
    if (isTRUE(dir.exists(weight))) {
      if (recover_strand_ambig_wgt){
        loginfo("Ignoring recover_strand_ambig_wgt=T; not implemented for FUSION format")
      }
      weightall <- read_weight_fusion(weight, b, ld_snpinfo, z_snp, method = method, 
                                      harmonize_wgt=harmonize_wgt)
    } else if (all(file_ext(weight) == "db")) {
      weightall <- read_weight_predictdb(weight, b, ld_snpinfo, z_snp, 
                                         harmonize_wgt=harmonize_wgt, ld_Rinfo=ld_Rinfo,
                                         recover_strand_ambig=recover_strand_ambig_wgt,
                                         ncore=ncore,
                                         scale_by_ld_variance=scale_by_ld_variance)
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
          ifreg <- ifelse(p1 >= ld_Rinfo[, "start"] & p0 < ld_Rinfo[, "stop"], T, F)
          exprlist[[gname]][["reg"]] <- paste(sort(ld_Rinfo$region_name[ifreg]), collapse = ";")
        }
        
        regs <- data.frame(gid = names(exprlist), reg = unlist(lapply(exprlist, "[[", "reg")), stringsAsFactors = F)
        batches <- names(sort(-table(regs$reg)))
        
        corelist <- lapply(1:ncore, function(core){batches_core <- batches[0:ceiling(length(batches)/ncore-1)*ncore+core]; batches_core[!is.na(batches_core)]})
        names(corelist) <- 1:ncore
        
        cl <- parallel::makeCluster(ncore, outfile = "")
        doParallel::registerDoParallel(cl)
        
        outlist <- foreach(core = 1:ncore, .combine = "c", .packages = c("ctwas", "tools")) %dopar% {
          
          batches <- corelist[[core]]
          
          outlist_core <- list()
          
          for (batch in batches) {
            gnames <- regs[regs$reg == batch, "gid"]
            regnames <- strsplit(batch, ";")[[1]]
            regRDS <- ld_Rinfo$RDS_file[match(regnames, ld_Rinfo$region_name)]
            R_snp <- lapply(regRDS, readRDS)
            R_snp <- suppressWarnings({Matrix::bdiag(R_snp)})
            R_snp_anno <- do.call(rbind, lapply(regRDS, read_ld_Rvar_RDS))
            for (i in 1:length(gnames)) {
              gname <- gnames[i]
              wgt <- exprlist[[gname]][["wgt"]]
              snpnames <- rownames(wgt)
              ld.idx <- match(snpnames, R_snp_anno$id)
              zdf.idx <- match(snpnames, z_snp$id)
              R.s <- R_snp[ld.idx, ld.idx]
              z.s <- as.matrix(z_snp$z[zdf.idx])
              z.g <- as.matrix(crossprod(wgt, z.s)/sqrt(t(wgt)%*%R.s%*% wgt))
              dimnames(z.g) <- NULL
              outlist_core[[gname]][["z.g"]] <- z.g
            }
          }
          outlist_core
        }
        
        parallel::stopCluster(cl)
        
        for (gname in names(outlist)){
          exprlist[[gname]][["z.g"]] <- outlist[[gname]][["z.g"]]
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
    
    gene_name <- lapply(exprlist, "[[", "gname")
    weight_name <- lapply(exprlist, "[[", "weight_name")
    
    if (length(exprlist) == 0) {
      geneinfo <- data.table::data.table(NULL)
    } else {
      geneinfo <- data.frame(chrom = chrom, id = gnames, p0 = p0, p1 = p1)
      geneinfo$gene_name <- gene_name
      geneinfo$weight_name <- weight_name
    }
    data.table::fwrite(geneinfo, file = exprvarf, sep = "\t", quote = F)
    z_gene_chr <- data.frame(id = gnames, z = z.g)
    exprqcf <- paste0(outname, "_chr", b, ".exprqc.Rd")
    save(wgtlist, qclist, z_gene_chr, file = exprqcf)
    exprf <- paste0(outname, "_chr", b, ".expr")
    if (!is.null(ld_pgenfs)) {
      if (length(exprlist) == 0) {
        expr <- data.table::data.table(NULL)
      } else {
        expr <- do.call(cbind, lapply(exprlist, "[[", "expr"))
      }
    } else {
      expr <- data.table::data.table(NA)
    }
    data.table::fwrite(expr, file = exprf, row.names = F, col.names = F, sep = "\t", quote = F)
    if (isTRUE(compress)) {
      system(paste0("gzip -f ", exprf))
      exprf <- paste0(exprf, ".gz")
    }
    loginfo("Imputation done: number of genes with imputed expression: %s for chr %s", length(gnames), b)
    z_genelist[[b]] <- z_gene_chr
    ld_exprfs[b] <- exprf
  }
  z_gene <- do.call(rbind, z_genelist)
  z_snp <- z_snp[z_snp$id %in% ld_snplist,] #subset z_snp to snps in ld reference
  return(list(z_gene = z_gene, ld_exprfs = ld_exprfs, z_snp = z_snp))
}
