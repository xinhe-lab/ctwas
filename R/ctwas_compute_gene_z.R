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

  loginfo("Computing gene z-scores...")
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
    z.g <- as.matrix(crossprod(wgt, z.s)/sqrt(t(wgt)%*%R.s%*%wgt))
    dimnames(z.g) <- NULL
    data.frame(id = id, z = z.g, type = type, context = context, group = group)
  }, mc.cores = ncore)

  z_gene <- do.call("rbind", z_gene)
  rownames(z_gene) <- NULL

  if (any(abs(z_gene$z) == Inf)) {
    logwarn("z_gene has Inf z-scores!")
    drop_idx <- which(abs(z_gene$z) == Inf)
    loginfo("Remove %d genes with Inf z-scores.", length(drop_idx))
    z_gene <- z_gene[-drop_idx, ]
  }

  return(z_gene)
}

#' @title Get gene info from weights.
#'
#' @param weights a list of preprocessed weights.
#'
#' @param ncore The number of cores used to parallelize over genes.
#'
#' @return a data frame of gene info.
#'
#' @export
get_gene_info <- function(weights, ncore = 1){

  if (!inherits(weights,"list")){
    stop("'weights' should be a list.")
  }
  gene_info <- mclapply_check(names(weights), function(x){
    as.data.frame(weights[[x]][c("chrom", "p0", "p1", "molecular_id", "weight_name")])
  }, mc.cores = ncore)
  gene_info <- do.call(rbind, gene_info)
  gene_info$id <- names(weights)
  gene_info <- gene_info[, c("chrom", "p0", "p1", "id", "molecular_id", "weight_name")]
  gene_info[, c("chrom", "p0", "p1")] <- sapply(gene_info[, c("chrom", "p0", "p1")], as.integer)
  rownames(gene_info) <- NULL

  return(gene_info)
}

#' @title Get regions for each gene.
#'
#' @param gene_info a data frame of gene info.
#'
#' @param region_info a data frame of region definitions.
#'
#' @param ncore The number of cores used to parallelize over genes.
#'
#' @return a data frame of gene info with regions overlapping each gene
#'
#' @importFrom logging loginfo logwarn
#'
#' @export
get_gene_regions <- function(gene_info,
                             region_info,
                             ncore = 1){

  gene_info <- gene_info[gene_info$chrom %in% unique(region_info$chrom), , drop=FALSE]

  if (!all(c("chrom", "p0", "p1") %in% colnames(gene_info))){
    stop("Columns chrom, p0, p1 are required in gene_info!")
  }

  gene_region_info_list <- mclapply_check(1:nrow(gene_info), function(i){
    tmp_gene_region_info <- gene_info[i, , drop=FALSE]
    chrom <- tmp_gene_region_info$chrom
    p0 <- tmp_gene_region_info$p0
    p1 <- tmp_gene_region_info$p1
    # find all the regions that overlap with the QTL range of the gene
    gene_regions.idx <- which(region_info$chrom == chrom & region_info$start <= p1 & region_info$stop > p0)
    if (length(gene_regions.idx) > 0) {
      tmp_gene_region_info$region_start <- min(region_info[gene_regions.idx,"start"])
      tmp_gene_region_info$region_stop <- max(region_info[gene_regions.idx,"stop"])
      tmp_gene_region_info$region_id <- paste(region_info[gene_regions.idx, "region_id"], collapse = ",")
      tmp_gene_region_info$n_regions <- length(gene_regions.idx)
    } else {
      logwarn("No regions overlapping with %s!", tmp_gene_region_info$id)
      tmp_gene_region_info$region_start <- NA
      tmp_gene_region_info$region_stop <- NA
      tmp_gene_region_info$region_id <- NA
      tmp_gene_region_info$n_regions <- length(gene_regions.idx)
    }
    tmp_gene_region_info
  }, mc.cores = ncore)

  gene_region_info <- do.call(rbind, gene_region_info_list)

  return(gene_region_info)
}

#' @title Get cross-boundary genes.
#'
#' @param region_info a data frame of region definitions.
#'
#' @param weights a list of preprocessed weights.
#'
#' @param gene_ids a vector of gene IDs (z_gene$id). If available, limits to these genes.
#'
#' @param mapping_table a data frame of mapping between molecular traits and genes,
#' with required columns: "molecular_id", "gene_name".
#'
#' @param map_by column name to be mapped by (default: "molecular_id").
#'
#' @param ncore The number of cores used to parallelize over genes.
#'
#' @return a data frame of cross-boundary genes.
#'
#' @importFrom logging loginfo
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join select summarise group_by ungroup rename
#' @importFrom rlang .data
#'
#' @export
get_boundary_genes <- function(region_info,
                               weights,
                               gene_ids = NULL,
                               mapping_table = NULL,
                               map_by = "molecular_id",
                               ncore = 1){

  loginfo("Get boundary genes...")

  # get gene info from weights
  gene_info <- get_gene_info(weights, ncore=ncore)

  # limit genes to gene_ids
  if (!is.null(gene_ids)) {
    gene_info <- gene_info[gene_info$id %in% gene_ids, , drop=FALSE]
  }

  gene_info_columns <- c("chrom", "p0", "p1", "id", "molecular_id", "weight_name")

  gene_info <- gene_info[, gene_info_columns]

  if (!is.null(mapping_table)) {
    # get combined gene_info
    # map molecular traits to genes, allow for many-to-many matching, return all matches
    # a gene could have molecular traits, and a molecular trait could be linked to multiple genes
    loginfo("Map molecular traits to genes")

    mapped_gene_info <- gene_info %>%
      left_join(mapping_table[,c(map_by, "gene_name")],
                by = map_by, multiple = "all", relationship = "many-to-many")

    gene_info_columns <- c(gene_info_columns, "gene_name")

    # get molecular trait level positions
    molecular_trait_level_info <- mapped_gene_info %>%
      group_by(.data$id) %>%
      summarise(chrom = .data$chrom[1],
                p0 = min(.data$p0),
                p1 = max(.data$p1),
                gene_name = paste(unique(.data$gene_name), collapse = ","),
                molecular_id = paste(unique(.data$molecular_id), collapse = ","),
                weight_name = paste(unique(.data$weight_name), collapse = ",")) %>%
      ungroup() %>%
      select(all_of(gene_info_columns)) %>%
      as.data.frame()

    # get gene level positions
    gene_level_info <- mapped_gene_info %>%
      group_by(.data$gene_name) %>%
      summarise(chrom = .data$chrom[1],
                id = paste(unique(.data$id), collapse = ","),
                p0 = min(.data$p0),
                p1 = max(.data$p1),
                molecular_id = paste(unique(.data$molecular_id), collapse = ","),
                weight_name = paste(unique(.data$weight_name), collapse = ",")) %>%
      ungroup() %>%
      select(all_of(gene_info_columns)) %>%
      as.data.frame()

    # combine gene_info from both molecular trait level and gene level
    gene_info <- rbind(molecular_trait_level_info, gene_level_info)
    gene_info <- unique(gene_info)
  }

  # get regions for each molecular trait or gene
  gene_region_info <- get_gene_regions(gene_info,
                                       region_info,
                                       ncore = ncore)

  # get boundary genes (n_regions > 1)
  boundary_genes <- gene_region_info[gene_region_info$n_regions > 1, ]
  boundary_genes <- boundary_genes[with(boundary_genes, order(chrom, p0, p1)), ]
  rownames(boundary_genes) <- NULL

  return(boundary_genes)
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
    loginfo("Group sizes before filtering \n {%s}: {%s}",
            names(gene_group_size), gene_group_size)
    groups_dropped <- names(gene_group_size)[gene_group_size < min_group_size]
    loginfo("Drop groups with group size < %d: %s", min_group_size, groups_dropped)
    z_gene <- z_gene[!z_gene$group %in% groups_dropped, , drop=FALSE]
    if (nrow(z_gene) == 0){
      stop("No genes left after group size filtering!")
    }
    gene_group_size <- table(z_gene$group)
    loginfo("Group sizes after filtering \n {%s}: {%s}",
            names(gene_group_size), gene_group_size)
  }
  return(z_gene)
}
