#' Compute gene z-scores
#'
#' @param z_snp A data frame with columns: "id", "z", giving the z-scores for SNPs.
#'
#' @param weights a list of weights
#'
#' @param logfile the log file, if NULL will print log info on screen
#'
#' @return a data frame of gene z-scores
#'
#' @importFrom logging addHandler loginfo writeToFile
#'
#' @export
compute_gene_z <- function (z_snp, weights, logfile = NULL){

  if (!is.null(logfile)) {
    addHandler(writeToFile, file = logfile, level = "DEBUG")
  }

  z_snp <- z_snp[,c("id", "z")]

  # compute gene z-scores
  # return a data frame with gene ids, and imputed gene z-scores
  z_gene <- data.frame()
  weight_ids <- names(weights)
  for (id in weight_ids) {
    wgt <- weights[[id]][["wgt"]]
    snpnames <- rownames(wgt)
    R.s <- weights[[id]][["LD_matrix"]]
    z.idx <- match(snpnames, z_snp$id)
    z.s <- as.matrix(z_snp$z[z.idx])
    z.g <- as.matrix(crossprod(wgt, z.s)/sqrt(t(wgt)%*%R.s%*% wgt))
    dimnames(z.g) <- NULL
    z_gene <- rbind(z_gene, data.frame(id = id, z = z.g))
  }
  rownames(z_gene) <- NULL

  return(z_gene)
}


#' add gene info to z_gene
add_gene_info <- function(z_gene, weights, region_info){
  # # find the regions overlapping with each gene
  # weight_info <- weight_info[weight_info$chrom %in% region_info$chrom, ]
  # for (i in 1:nrow(weight_info)) {
  #   chrom <- weight_info[i, "chrom"]
  #   p0 <- weight_info[i, "p0"]
  #   p1 <- weight_info[i, "p1"]
  #   idx <- which(region_info$chrom == chrom & region_info$start <= p1 & region_info$stop > p0)
  #   weight_info[i, "region_tag"] <- paste(sort(region_info[idx, "region_tag"]), collapse = ";")
  #   weight_info[i, "cross_boundary"] <- ifelse(length(idx) > 1, 1, 0)
  # }
  # gene_info <- weight_info[match(z_gene$id, weight_info$id), c("chrom", "id", "p0", "p1")]
  return(gene_info)
}

#' Combine SNP and gene z-scores
#'
#' @param z_snp a data frame with four columns: "id", "A1", "A2", "z".
#' giving the z scores for SNPs. "A1" is effect allele. "A2" is the other allele.
#'
#' @param z_gene a data frame with two columns: "id", "z". giving the z scores for genes.
#' Optionally, a "type" column can also be supplied; this is for using multiple sets of weights
#'
#' @export
#'
combine_z <- function(z_snp, z_gene){
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
  rownames(zdf) <- NULL
  return(zdf)
}
