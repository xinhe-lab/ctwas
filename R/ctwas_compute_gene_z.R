#' @title Computes gene z-scores.
#'
#' @param z_snp A data frame with columns: "id", "z", giving the z-scores for SNPs.
#'
#' @param weights a list of weights
#'
#' @param ncore The number of cores used to parallelize computation over weights
#'
#' @param logfile the log file, if NULL will print log info on screen
#'
#' @return a data frame of gene z-scores,
#'
#' @importFrom logging addHandler loginfo writeToFile
#' @importFrom parallel mclapply
#'
#' @export
compute_gene_z <- function (z_snp,
                            weights,
                            ncore = 1,
                            logfile = NULL){

  if (!is.null(logfile)) {
    addHandler(writeToFile, file = logfile, level = "DEBUG")
  }

  # check input data
  if (!inherits(weights,"list")){
    stop("'weights' should be a list.")
  }

  if (any(sapply(weights, is.null))) {
    warning("weights contain NULL, remove empty weights.")
    loginfo("Remove empty weights")
    weights[which(sapply(weights, is.null))] <- NULL
  }

  loginfo("Computing gene z-scores ...")

  weight_snpnames <- unique(unlist(lapply(weights, function(x){rownames(x[["wgt"]])})))
  z_snp <- z_snp[z_snp$id %in% weight_snpnames, c("id", "z")]

  z_gene <- mclapply_check(names(weights), function(id) {
    wgt <- weights[[id]][["wgt"]]
    snpnames <- rownames(wgt)
    R.s <- weights[[id]][["R_wgt"]]
    z.idx <- match(snpnames, z_snp$id)
    z.s <- as.matrix(z_snp$z[z.idx])
    z.g <- as.matrix(crossprod(wgt, z.s)/sqrt(t(wgt)%*%R.s%*% wgt))
    dimnames(z.g) <- NULL
    data.frame(id = id, z = z.g)
  }, mc.cores = ncore)

  z_gene <- do.call("rbind", z_gene)
  rownames(z_gene) <- NULL

  if (anyNA(z_gene)){
    stop("z_gene contains missing values!\n")
  }

  return(z_gene)
}

# get gene info from weights
get_gene_info <- function(weights){
  # check input data
  if (!inherits(weights,"list")){
    stop("'weights' should be a list.")
  }
  gene_info <- lapply(names(weights), function(x){
    as.data.frame(weights[[x]][c("chrom", "p0","p1", "gene_name", "weight_name")])})
  gene_info <- do.call(rbind, gene_info)
  gene_info$id <- names(weights)
  gene_info <- gene_info[, c("chrom", "id", "p0", "p1", "gene_name", "weight_name")]
  gene_info[, c("chrom","p0", "p1")] <- sapply(gene_info[, c("chrom","p0", "p1")], as.integer)
  # gene_info <- gene_info[with(gene_info, order(chrom, p0)), ]
  rownames(gene_info) <- NULL

  return(gene_info)
}

# get regions for each gene
get_gene_regions <- function(gene_info, region_info){

  gene_info <- gene_info[gene_info$chrom %in% unique(region_info$chrom), , drop=FALSE]
  for (i in 1:nrow(gene_info)) {
    chrom <- gene_info[i, "chrom"]
    p0 <- gene_info[i, "p0"]
    p1 <- gene_info[i, "p1"]
    region_idx <- which(region_info$chrom == chrom & region_info$start <= p1 & region_info$stop > p0)
    gene_info[i, "region_start"] <- min(region_info[region_idx,"start"])
    gene_info[i, "region_stop"] <- max(region_info[region_idx,"stop"])
    gene_info[i, "region_id"] <- paste(sort(region_info[region_idx, "region_id"]), collapse = ";")
    gene_info[i, "n_regions"] <- length(region_idx)
  }
  return(gene_info)
}

# Combines SNP and gene z-scores
#
# @param z_snp a data frame of SNP z-scores, with columns: "id", "z".
#
# @param z_gene a data frame of gene z-scores, with columns: "id", "z".
#
combine_z <- function(z_snp, z_gene){
  z_snp$group <- "SNP"
  z_gene$group <- "gene"
  zdf <- rbind(z_snp[, c("id", "z", "group")],
               z_gene[, c("id", "z", "group")])
  rownames(zdf) <- NULL
  return(zdf)
}
