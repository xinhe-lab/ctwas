#' Compute gene z-scores
#'
#' @param z_snp A data frame with columns: "id", "A1", "A2", "z". giving the z scores for
#' snps. "A1" is effect allele. "A2" is the other allele. For harmonized data, A1 and A2 are not required.
#'
#' @param weight a string, pointing to a directory with the weights (.db file) in predictdb format.
#' A vector of multiple sets of weights in PredictDB format can also be specified;
#' genes will have their filename appended to their gene name to ensure IDs are unique.
#'
#' @param region_info a data frame of region definition and associated file names.
#'
#' @param ncore The number of cores used to parallelize imputation over weights
#'
#' @param chr a numeric vector of chromosomes to perform z score imputation over. Useful for large jobs requiring batches
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
                            ncore=1,
                            chr=1:22,
                            logfile = NULL){

  if (!is.null(logfile)) {
    addHandler(writeToFile, file = logfile, level = "DEBUG")
  }

  # read and check Rvar files from the region_info table,
  # return a list of updated region_info table (with region names) for each chromosome
  ld_Rinfo_list <- read_region_ld_Rinfo(region_info)

  z_gene_list <- vector("list", length = 22)
  gene_info_list <- vector("list", length = 22)
  ld_ref_snps <- c() # store all SNPs in LD reference

  for (b in chr) {

    # updated region info table in the chromosome
    ld_Rinfo <- ld_Rinfo_list[[b]]

    if(!is.null(ld_Rinfo)){
      loginfo("Impute gene z scores for chromosome %s", b)

      # read SNP info from Rvar file of all the regions in the chromosome
      # ld_snpinfo <- read_ld_Rvar_snp_info(ld_Rinfo$Rvar_file)
      ld_snpinfo <- do.call(rbind, lapply(ld_Rinfo$Rvar_file, data.table::fread))

      chrom <- unique(ld_snpinfo$chrom)
      if (length(chrom) > 1) {
        stop("Input LD reference not split by chromosome")
      }
      ld_ref_snps <- c(ld_ref_snps, ld_snpinfo$id) # store names of SNPs in LD reference

      loginfo("Reading weights for chromosome %s", b)
      weightall <- read_weight_predictdb(weight,
                                         chrom = b,
                                         ld_snpinfo = ld_snpinfo,
                                         z_snp = z_snp,
                                         harmonize_wgt=FALSE,
                                         ld_Rinfo = ld_Rinfo,
                                         ncore = ncore)

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

        outlist <- foreach(core = 1:ncore, .combine = "c", .packages = c("ctwas")) %dopar% {

          batches <- corelist[[core]]
          outlist_core <- list()

          for (batch in batches) {
            gnames <- regs[regs$reg == batch, "gid"]
            regnames <- strsplit(batch, ";")[[1]]
            reg_idx <- match(regnames, ld_Rinfo$region_name)
            reg_RDS_files <- ld_Rinfo$RDS_file[reg_idx]
            reg_Rvar_files <- ld_Rinfo$Rvar_file[reg_idx]
            R_snp <- lapply(reg_RDS_files, readRDS)
            R_snp <- suppressWarnings({as.matrix(Matrix::bdiag(R_snp))})
            # R_snp_info <- do.call(rbind, lapply(reg_Rvar_files, read_ld_Rvar_file))
            R_snp_info <- do.call(rbind, lapply(reg_Rvar_files, data.table::fread))

            for (i in 1:length(gnames)) {
              gname <- gnames[i]
              wgt <- exprlist[[gname]][["wgt"]]
              snpnames <- rownames(wgt)
              ld.idx <- match(snpnames, R_snp_info$id)
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
      loginfo("Number of genes with imputed expression: %s for chr %s", length(gnames), b)

      if (length(gnames) > 0) {
        gene_info_chr <- data.frame(chrom = chrom,
                                    id = gnames,
                                    p0 = p0,
                                    p1 = p1,
                                    gene_name = gene_name,
                                    weight_name = weight_name)
      } else {
        gene_info_chr <- data.table::data.table(NULL)
      }

      gene_info_list[[b]] <- gene_info_chr

      # data frame with gene ids, and imputed gene z-scores
      z_gene_list[[b]] <- data.frame(id = gnames, z = z.g)
    }

  }

  # gene z-score data frame with gene ids, and imputed gene z-scores
  z_gene <- do.call(rbind, z_gene_list)

  # gene info data frame with gene ids, gene names, gene coordinates and weight names
  gene_info <- do.call(rbind, gene_info_list)

  # SNP z-score data frame, only include SNPs in LD reference
  z_snp <- z_snp[z_snp$id %in% ld_ref_snps,]

  return(list(z_gene = z_gene,
              z_snp = z_snp,
              gene_info = gene_info))
}

