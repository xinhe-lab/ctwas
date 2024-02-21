#' Compute gene z-scores
#'
#' @param z_snp A data frame with columns: "id", "A1", "A2", "z". giving the z scores for
#' snps. "A1" is effect allele. "A2" is the other allele. For harmonized data, A1 and A2 are not required.
#'
#' @param weight a string, pointing to a directory with the FUSION/TWAS format of weights, or a .db file in predictdb format.
#' A vector of multiple sets of weights in PredictDB format can also be specified; genes will have their filename appended
#' to their gene name to ensure IDs are unique.
#'
#' @param region_info a data frame of region definition and associated file names
#'
#' @param weight_format a string, the format of weight, PredictDB or FUSION
#'
#' @param method a string, blup/bslmm/lasso/top1/enet/best. This option is only used for FUSION weights.
#' "best" means the method giving the best cross #' validation R^2. Note that top1 uses only the weight
#' with largest effect.
#'
#' @param ncore The number of cores used to parallelize imputation over weights
#'
#' @param chr a numeric vector of chromosomes to perform z score imputation over. Useful for large jobs requiring batches
#'
#' @param scale_by_ld_variance TRUE/FALSE. If TRUE, PredictDB weights are scaled by genotype variance, which is the default
#' behavior for PredictDB
#'
#' @param logfile the log file, if NULL will print log info on screen
#'
#' @return a list of gene z-scores and SNP z-scores
#'
#' @importFrom logging addHandler loginfo
#'
#' @export
compute_gene_z <- function (z_snp,
                            weight,
                            region_info,
                            weight_format = c("PredictDB", "FUSION"),
                            method = c("lasso", "blup", "bslmm", "top1", "enet", "best"),
                            ncore=1,
                            chr=1:22,
                            scale_by_ld_variance=TRUE,
                            logfile = NULL){

  if (!is.null(logfile)) {
    addHandler(writeToFile, file = logfile, level = "DEBUG")
  }

  weight_format <- match.arg(weight_format)
  method <- match.arg(method)

  # read and check Rvar from the region_info table
  ld_Rinfo_list <- read_region_ld_Rinfo(region_info)

  z_genelist <- list()
  ld_snplist <- c() # list to store names of snps in ld reference

  for (b in chr) {

    # region info in the chromosome
    ld_Rinfo <- ld_Rinfo_list[[b]]

    if(!is.null(ld_Rinfo)){
      loginfo("Impute gene z scores for chromosome %s", b)

      # read snp info in all the regions in the chromosome
      ld_snpinfo <- read_ld_Rvar_snp_info(ld_Rinfo$R_snp_info)

      chrom <- unique(ld_snpinfo$chrom)
      if (length(chrom) > 1) {
        stop("Input LD reference not split by chromosome")
      }
      ld_snplist <- c(ld_snplist, ld_snpinfo$id) #store names of snps in ld reference

      loginfo("Reading weights with %s format for chromosome %s", weight_format, b)
      if (weight_format == "FUSION") {
        weightall <- read_weight_FUSION(weight,
                                        b,
                                        ld_snpinfo,
                                        z_snp,
                                        method = method,
                                        harmonize_wgt=F)
      } else if (weight_format == "PredictDB") {
        weightall <- read_weight_predictdb(weight,
                                           b,
                                           ld_snpinfo,
                                           z_snp,
                                           harmonize_wgt=F,
                                           ld_Rinfo=ld_Rinfo,
                                           ncore=ncore,
                                           scale_by_ld_variance=scale_by_ld_variance)
      } else {
        stop("Unrecognized weight format, need to use either FUSION format or predict.db format")
      }

      exprlist <- weightall[["exprlist"]]
      qclist <- weightall[["qclist"]]
      if (length(exprlist) > 0) {
        loginfo("Start gene z score imputation ...")
        loginfo("Using given LD matrices to impute gene z score.")

        for (gname in names(exprlist)) {
          p0 <- exprlist[[gname]][["p0"]]
          p1 <- exprlist[[gname]][["p1"]]
          ifreg <- ifelse(p1 >= ld_Rinfo[, "start"] & p0 < ld_Rinfo[, "stop"], T, F)
          exprlist[[gname]][["reg"]] <- paste(sort(ld_Rinfo$region_name[ifreg]), collapse = ";")
        }

        regs <- data.frame(gid = names(exprlist), reg = unlist(lapply(exprlist, "[[", "reg")), stringsAsFactors = F)
        batches <- names(sort(-table(regs$reg)))

        corelist <- lapply(1:ncore, function(core){
          batches_core <- batches[0:ceiling(length(batches)/ncore-1)*ncore+core];
          batches_core[!is.na(batches_core)]})
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
            R_snp <- suppressWarnings({as.matrix(Matrix::bdiag(R_snp))})
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
      loginfo("Gene z score imputation done.")
      gnames <- names(exprlist)
      z.g <- unlist(lapply(exprlist, "[[", "z.g"))
      chrom <- unlist(lapply(exprlist, "[[", "chrom"))
      p0 <- unlist(lapply(exprlist, "[[", "p0"))
      p1 <- unlist(lapply(exprlist, "[[", "p1"))
      wgtlist <- lapply(exprlist, "[[", "wgt")
      gene_name <- lapply(exprlist, "[[", "gname")
      weight_name <- lapply(exprlist, "[[", "weight_name")

      if (length(exprlist) == 0) {
        gene_info <- data.table::data.table(NULL)
      } else {
        gene_info <- data.frame(chrom = chrom, id = gnames, p0 = p0, p1 = p1)
        gene_info$gene_name <- gene_name
        gene_info$weight_name <- weight_name
      }

      z_gene_chr <- data.frame(id = gnames, z = z.g)
      loginfo("Number of genes with imputed expression: %s for chr %s", length(gnames), b)
      z_genelist[[b]] <- z_gene_chr
    }

  }
  z_gene <- do.call(rbind, z_genelist)

  # filter z_snp to snps in LD reference
  z_snp <- z_snp[z_snp$id %in% ld_snplist,]

  # combine z-scores of genes and SNPs
  zdf <- combine_z(z_gene, z_snp)

  return(list(z_gene = z_gene,
              z_snp = z_snp,
              zdf = zdf,
              gene_info = gene_info))
}

# combine z-scores of genes and SNPs
combine_z <- function(z_gene, z_snp){

  z_snp$type <- "SNP"
  z_snp$QTLtype <- "SNP"
  if (is.null(z_gene$type)){
    z_gene$type <- "gene"
  }
  if (is.null(z_gene$QTLtype)){
    z_gene$QTLtype <- "gene"
  }
  zdf <- rbind(z_snp[, c("id", "z", "type", "QTLtype")],
               z_gene[, c("id", "z", "type", "QTLtype")])

  return(zdf)
}

