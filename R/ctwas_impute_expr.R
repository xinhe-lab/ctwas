#' Impute expression
#' 
#' @param pgenfs A character vector of .pgen or .bed files. One file for one
#' chromosome, in the order of 1 to 22. Therefore, the length of this vector
#' needs to be 22. If .pgen files are given, then .pvar and .psam are assumed
#' to present in the same directory. If .bed files are given, then .bim and
#' .fam files are assumed to present in the same directory.
#' 
#' @param weight a string, pointing to a directory with the fusion/twas format of weights, 
#' or a vector of one or more .db files in PredictDB format.
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
#' @param harmonize_wgt TRUE/FALSE. If TRUE, GWAS and eQTL genotype alleles are harmonized
#' 
#' @param strand_ambig_action_wgt the action to take to harmonize strand ambiguous variants (A/T, G/C) between 
#' the weights and LD reference. "drop" removes the ambiguous variant from the prediction models. "none" treats the variant 
#' as unambiguous, flipping the weights to match the LD reference and then taking no additional action. "recover" uses a procedure
#' to recover strand ambiguous variants. This procedure compares correlations between variants in the 
#' LD reference and prediction models, and it can only be used with PredictDB format prediction models, which include this
#' information.
#' 
#' @param ncore The number of cores used to parallelize imputation over weights
#' 
#' @param scale_by_ld_variance TRUE/FALSE. If TRUE, PredictDB weights are scaled by genotype variance, which is the default 
#' behavior for PredictDB
#' 
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
                        harmonize_wgt = T,
                        ncore=1,
                        scale_by_ld_variance=T,
                        strand_ambig_action_wgt = c("drop", "none", "recover")){

  if (!is.null(logfile)){
    addHandler(writeToFile, file= logfile, level='DEBUG')
  }

  if (isTRUE(harmonize_wgt)){
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

    loginfo("Reading weights for chromosome %s", b)
    
    if (isTRUE(dir.exists(weight))) {
      weight_type <- "fusion"
      weightall <- read_weight_fusion(weight, 
                                      b, 
                                      snpinfo, 
                                      z_snp=NULL, 
                                      method = method, 
                                      harmonize_wgt=harmonize_wgt)
    } else if (all(file_ext(weight) == "db")) {
      weight_type <- "predictdb"
      weightall <- read_weight_predictdb(weight, 
                                         b, 
                                         snpinfo, 
                                         z_snp=NULL, 
                                         harmonize_wgt=harmonize_wgt, 
                                         strand_ambig_action=strand_ambig_action_wgt,
                                         ncore=ncore,
                                         scale_by_ld_variance=scale_by_ld_variance)
    } else {
      stop("Unrecognized weight format, need to use either FUSION format or predict.db format")
    }

    exprlist <- weightall[["exprlist"]]
    qclist <- weightall[["qclist"]]
    
    if (length(exprlist) > 0) {
      loginfo("Start gene experssion imputation ...")
      for (gname in names(exprlist)){
        wgt <- exprlist[[gname]][["wgt"]]
        snpnames <- rownames(wgt)
        gwas.idx <-  match(snpnames, snpinfo$id)
        g <- read_pgen(pgen, variantidx = gwas.idx)
        
        if (weight_type=="fusion"){
          g <- scale(g)  # genotypes are standardized for FUSION
        }
        
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
