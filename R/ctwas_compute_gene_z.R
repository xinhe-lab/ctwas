#' @title Computes z-scores of molecular traits
#'
#' @param z_snp A data frame with columns: "id", "z", giving the z-scores for SNPs.
#'
#' @param weights a list of preprocessed weights.
#'
#' @param ncore The number of cores used to parallelize computation over weights.
#'
#' @param logfile The log filename. If NULL, print log info on screen.
#'
#' @return a data frame of z-scores of molecular traits
#'
#' @importFrom logging addHandler loginfo logwarn writeToFile
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
    stop("weights contain NULL, remove empty weights!")
  }

  if (anyNA(z_snp)){
    stop("z_snp contains missing values!")
  }

  loginfo("Computing gene z-scores ...")
  weight_snp_ids <- unique(unlist(lapply(weights, function(x){rownames(x[["wgt"]])})))

  if (any(!weight_snp_ids %in% z_snp$id)){
    stop(paste("Some SNPs in weights not found in z_snp!\n",
         "Please run preprocess_z_snp() before running preprocess_weights() with 'gwas_snp_ids = z_snp$id'!"))
  }

  z_snp <- z_snp[z_snp$id %in% weight_snp_ids, c("id", "z")]

  z_gene <- mclapply_check(names(weights), function(id) {
    wgt <- weights[[id]][["wgt"]]
    wgt_snp_ids <- rownames(wgt)
    R.s <- weights[[id]][["R_wgt"]]
    type <- weights[[id]][["type"]]
    context <- weights[[id]][["context"]]
    group <- paste0(context,"|",type)
    z.idx <- match(wgt_snp_ids, z_snp$id)
    z.s <- as.matrix(z_snp$z[z.idx])
    z.g <- as.matrix(crossprod(wgt, z.s)/sqrt(t(wgt)%*%R.s%*% wgt))
    dimnames(z.g) <- NULL
    data.frame(id = id, z = z.g, type = type, context = context, group = group)
  }, mc.cores = ncore)

  z_gene <- do.call("rbind", z_gene)
  rownames(z_gene) <- NULL

  return(z_gene)
}

#' @title Get gene info from weights
#'
#' @param weights a list of preprocessed weights
#'
#' @return a data frame of gene info
#'
#' @export
get_gene_info <- function(weights){

  if (!inherits(weights,"list")){
    stop("'weights' should be a list.")
  }
  gene_info <- lapply(names(weights), function(x){
    as.data.frame(weights[[x]][c("chrom", "p0","p1", "molecular_id", "weight_name")])
  })
  gene_info <- do.call(rbind, gene_info)
  gene_info$id <- names(weights)
  gene_info <- gene_info[, c("chrom", "id", "p0", "p1", "molecular_id", "weight_name")]
  gene_info[, c("chrom", "p0", "p1")] <- sapply(gene_info[, c("chrom", "p0", "p1")], as.integer)
  rownames(gene_info) <- NULL

  return(gene_info)
}

#' @title Get regions for each gene
#'
#' @param gene_info a data frame of gene info
#'
#' @param region_info a data frame of region definitions
#'
#' @return a data frame of gene info with regions overlapping each gene
#'
#' @export
get_gene_regions <- function(gene_info, region_info){

  gene_info <- gene_info[gene_info$chrom %in% unique(region_info$chrom), , drop=FALSE]
  for (i in 1:nrow(gene_info)) {
    chrom <- gene_info[i, "chrom"]
    p0 <- gene_info[i, "p0"]
    p1 <- gene_info[i, "p1"]
    region.idx <- which(region_info$chrom == chrom & region_info$start <= p1 & region_info$stop > p0)
    if (length(region.idx) > 0) {
      gene_info[i, "region_start"] <- min(region_info[region.idx,"start"])
      gene_info[i, "region_stop"] <- max(region_info[region.idx,"stop"])
      gene_info[i, "region_id"] <- paste(sort(region_info[region.idx, "region_id"]), collapse = ",")
      gene_info[i, "n_regions"] <- length(region.idx)
    } else {
      warning(paste("No regions overlapping with", gene_info[i, "id"]))
      gene_info[i, "region_start"] <- NA
      gene_info[i, "region_stop"] <- NA
      gene_info[i, "region_id"] <- NA
      gene_info[i, "n_regions"] <- length(region.idx)
    }
  }
  return(gene_info)
}

# Combines SNP and gene z-scores
#
# @param z_snp a data frame of SNP z-scores, with columns: "id", "A1", "A2", "z".
# ("A1" and "A2" are optional)
#
# @param z_gene a data frame of gene z-scores, with columns: "id", "z", "type",
# "context", "group".
#
combine_z <- function(z_snp, z_gene){

  if (is.null(z_gene$type)){
    z_gene$type <- "gene"
  }
  if (is.null(z_gene$context)){
    z_gene$context <- "gene"
  }
  if (is.null(z_gene$group)){
    z_gene$group <- "gene"
  }

  z_snp$type <- "SNP"
  z_snp$context <- "SNP"
  z_snp$group <- "SNP"

  z_df <- rbind(z_gene[, c("id", "z", "type", "context", "group")],
                z_snp[, c("id", "z", "type", "context", "group")])
  rownames(z_df) <- NULL
  return(z_df)
}


#' @title Filter z_gene by group size
#'
#' @param z_gene a data frame of gene z-scores, with columns: "id", "z", "type",
#' "context", "group".
#'
#' @param min_group_size Minimum number of variables in a group.
#'
#' @return a data frame of gene z-scores.
#'
#' @export
filter_z_gene_by_group_size <- function(z_gene, min_group_size){
  gene_group_size <- table(z_gene$group)
  if (any(gene_group_size < min_group_size)){
    groups_dropped <- names(gene_group_size)[gene_group_size < min_group_size]
    loginfo("Filter groups with group size < %d from z_gene", min_group_size)
    z_gene <- z_gene[!z_gene$group %in% groups_dropped, , drop=FALSE]
    if (nrow(z_gene) == 0){
      stop("No genes left after group size filtering!")
    }
    gene_group_size <- table(z_gene$group)
    loginfo("Group sizes in z_gene after filtering \n {%s}: {%s}",
            names(gene_group_size), gene_group_size)
  }
  return(z_gene)
}
