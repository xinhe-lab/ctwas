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

  if (is.null(region_info$region_tag)){
    region_info$region_tag <- paste0(region_info$chrom, ":", region_info$start, "-", region_info$stop)
  }

  # find the regions overlapping with each gene
  weight_info <- weight_info[weight_info$chrom %in% region_info$chrom, ]
  for (i in 1:nrow(weight_info)) {
    chrom <- weight_info[i, "chrom"]
    p0 <- weight_info[i, "p0"]
    p1 <- weight_info[i, "p1"]
    idx <- which(region_info$chrom == chrom & region_info$start <= p1 & region_info$stop > p0)
    weight_info[i, "region_tag"] <- paste(sort(region_info[idx, "region_tag"]), collapse = ";")
    weight_info[i, "cross_boundary"] <- ifelse(length(idx) > 1, 1, 0)
  }

  z_snp <- z_snp[,c("id", "z")]

  # impute gene z-scores for each chromosome
  cl <- parallel::makeCluster(ncore, outfile = "")
  doParallel::registerDoParallel(cl)

  for (b in unique(region_info$chrom)) {
    loginfo("Impute gene z scores for chr%s", b)
    weightinfo <- weight_info[weight_info$chrom == b, ]

    if (nrow(weightinfo) > 0) {
      batches <- names(sort(-table(weightinfo$region_tag)))

      corelist <- lapply(1:ncore, function(core){
        batches_core <- batches[0:ceiling(length(batches)/ncore-1)*ncore+core];
        batches_core[!is.na(batches_core)]})
      names(corelist) <- 1:ncore

      outlist <- foreach(core = 1:ncore, .combine = "c", .packages = c("ctwas")) %dopar% {

        batches <- corelist[[core]]
        outlist_core <- list()

        for (batch in batches) {
          # load the R_snp and SNP info for the region
          region_tags <- strsplit(batch, ";")[[1]]
          reg_idx <- match(region_tags, region_info$region_tag)
          if (length(reg_idx) > 1){
            R_snp <- lapply(region_info$LD_matrix[reg_idx], load_LD)
            R_snp <- suppressWarnings({as.matrix(Matrix::bdiag(R_snp))})
          }else{
            R_snp <- load_LD(region_info$LD_matrix[reg_idx])
          }
          ld_snpinfo <- read_LD_SNP_files(region_info$SNP_info[reg_idx])

          # compute gene z-scores
          ids <- weightinfo[weightinfo$region_tag == batch, "id"]
          for (id in ids) {
            wgt <- weight_list[[id]]
            snpnames <- rownames(wgt)
            ld.idx <- match(snpnames, ld_snpinfo$id)
            z.idx <- match(snpnames, z_snp$id)
            R.s <- R_snp[ld.idx, ld.idx]
            z.s <- as.matrix(z_snp$z[z.idx])
            z.g <- as.matrix(crossprod(wgt, z.s)/sqrt(t(wgt)%*%R.s%*% wgt))
            dimnames(z.g) <- NULL
            outlist_core[[id]] <- data.frame(id = id, z = z.g)
          }
        }
        outlist_core
      }
    }
    loginfo("Number of genes with imputed expression: %d for chr%s", length(outlist), b)
    z_gene_list[[b]] <- do.call(rbind, outlist)
  }

  parallel::stopCluster(cl)

  # gene z-score data frame with gene ids, and imputed gene z-scores
  z_gene <- do.call(rbind, z_gene_list)

  # gene info data frame with gene ids, gene names, gene coordinates and weight names
  gene_info <- weight_info[match(z_gene$id, weight_info$id),
                           c("chrom", "id", "p0", "p1", "gene_name", "weight_name", "region_tag", "cross_boundary")]
  rownames(gene_info) <- NULL

  return(list(z_gene = z_gene, gene_info = gene_info))
}

