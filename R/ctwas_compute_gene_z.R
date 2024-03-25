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
  for (id in names(weights)) {
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


#' get gene info
get_gene_info <- function(weights, region_info=NULL){

  gene_info <- lapply(names(weights), function(x){
    as.data.frame(weights[[x]][c("chrom", "p0","p1", "gene_name", "weight_name")])})
  gene_info <- do.call(rbind, gene_info)
  gene_info$id <- names(weights)
  gene_info <- gene_info[, c("chrom", "id", "p0", "p1", "gene_name", "weight_name")]
  rownames(gene_info) <- NULL

  # find the regions overlapping with each gene
  if (!is.null(region_info)) {
    for (i in 1:nrow(gene_info)) {
      chrom <- gene_info[i, "chrom"]
      p0 <- gene_info[i, "p0"]
      p1 <- gene_info[i, "p1"]
      region_idx <- which(region_info$chrom == chrom & region_info$start <= p1 & region_info$stop > p0)
      gene_info[i, "region_tag"] <- paste(sort(region_info[region_idx, "region_tag"]), collapse = ";")
      gene_info[i, "n_regions"] <- length(region_idx)
    }
  }

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
