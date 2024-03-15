#' Compute gene z-scores
#'
#' @param z_snp A data frame with columns: "id", "A1", "A2", "z". giving the z scores for
#' snps. "A1" is effect allele. "A2" is the other allele. For harmonized data, A1 and A2 are not required.
#'
#' @param region_info a data frame of region definition and associated file names.
#'
#' @param weight_list a list of weights
#'
#' @param weight_info a data frame of weight information
#'
#' @param ncore The number of cores used to parallelize imputation over weights
#'
#' @param logfile the log file, if NULL will print log info on screen
#'
#' @return a list of gene z-scores and gene info table
#'
#' @importFrom logging addHandler loginfo writeToFile
#'
#' @export
compute_gene_z <- function (z_snp,
                            region_info,
                            weight_list,
                            weight_info,
                            ncore=1,
                            logfile = NULL){

  if (!is.null(logfile)) {
    addHandler(writeToFile, file = logfile, level = "DEBUG")
  }

  z_gene_list <- list()

  # check to make sure all SNPs are in LD reference
  ld_snpinfo <- read_LD_SNP_files(region_info$SNP_info)
  if (!all(z_snp$id %in% ld_snpinfo$id)){
    stop("Not all SNPs are in LD reference. Have you done harmonization?")
  }
  rm(ld_snpinfo)

  if (is.null(region_info$region_tag)){
    region_info$region_tag <- paste0(region_info$chrom, ":", region_info$start, "-", region_info$stop)
  }

  for (b in unique(region_info$chrom)) {

    loginfo("Impute gene z scores for chromosome %s", b)

    # region info in the chromosome
    regioninfo <- region_info[region_info$chrom == b, ]
    regioninfo <- regioninfo[order(regioninfo$start), ]

    # genes in the chromosome
    weightinfo <- weight_info[weight_info$chrom == b, ]
    wgtlist <- weight_list[weightinfo$id]

    if (length(wgtlist) > 0) {
      loginfo("Start gene z score imputation ...")

      # fine the regions for each gene
      for (gname in weightinfo$id) {
        p0 <- weightinfo[gname, "p0"]
        p1 <- weightinfo[gname, "p1"]
        ifreg <- ifelse(p1 >= regioninfo[, "start"] & p0 < regioninfo[, "stop"], T, F)
        weightinfo[gname, "region_tag"] <- paste(sort(regioninfo$region_tag[ifreg]), collapse = ";")
      }

      regs <- data.frame(gname = weightinfo$id,
                         region_tag = weightinfo$region_tag,
                         stringsAsFactors = F)

      batches <- names(sort(-table(regs$region_tag)))

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
          gnames <- regs[regs$region_tag == batch, "gname"]
          region_tags <- strsplit(batch, ";")[[1]]
          reg_idx <- match(region_tags, regioninfo$region_tag)
          R_snp <- lapply(regioninfo$LD_matrix[reg_idx], load_LD)
          R_snp <- suppressWarnings({as.matrix(Matrix::bdiag(R_snp))})
          ld_snpinfo <- read_LD_SNP_files(regioninfo$SNP_info[reg_idx])

          for (gname in gnames) {
            wgt <- wgtlist[[gname]]
            snpnames <- rownames(wgt.matrix)
            # keep SNPs in weight, z_snp and LD reference
            snpnames <- Reduce(intersect, list(snpnames, z_snp$id, ld_snpinfo$id))
            ld.idx <- match(snpnames, snpinfo$id)
            z.idx <- match(snpnames, z_snp$id)
            R.s <- R_snp[ld.idx, ld.idx]
            z.s <- as.matrix(z_snp$z[z.idx])
            z.g <- as.matrix(crossprod(wgt, z.s)/sqrt(t(wgt)%*%R.s%*% wgt))
            dimnames(z.g) <- NULL
            outlist_core[[gname]][["z.g"]] <- z.g
          }
        }
        outlist_core
      }
      parallel::stopCluster(cl)
    }

    loginfo("Gene z-score imputation done.")
    z.g <- unlist(lapply(outlist, "[[", "z.g"))
    gnames <- names(outlist)
    loginfo("Number of genes with imputed expression: %d for chr%s", length(gnames), b)

    # data frame with gene|weight ids, and imputed gene z-scores
    z_gene_list[[b]] <- data.frame(id = gnames, z = z.g)

    # if (length(gnames) > 0) {
    #   gene_info_chr <- weight_info[gnames, c("chrom", "id", "p0", "p1", "gene_name", "weight_name")]
    # } else {
    #   gene_info_chr <- data.table::data.table(NULL)
    # }
    #
    # gene_info_list[[b]] <- gene_info_chr

  }

  # gene z-score data frame with gene ids, and imputed gene z-scores
  z_gene <- do.call(rbind, z_gene_list)

  # gene info data frame with gene ids, gene names, gene coordinates and weight names
  # gene_info <- do.call(rbind, gene_info_list)
  gene_info <- weight_info[match(z_gene$id, weight_info$id), c("chrom", "id", "p0", "p1", "gene_name", "weight_name")]

  return(list(z_gene = z_gene, gene_info = gene_info))
}

