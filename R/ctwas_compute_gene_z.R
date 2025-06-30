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
#' @return a data frame of gene info.
#'
#' @export
get_gene_info <- function(weights){

  if (!inherits(weights,"list"))
    stop("'weights' should be a list.")

  gene_info <- data.frame(chrom = as.integer(sapply(weights, "[[", "chrom")),
                          p0 = as.integer(sapply(weights, "[[", "p0")),
                          p1 = as.integer(sapply(weights, "[[", "p1")),
                          id = names(weights),
                          molecular_id = sapply(weights, "[[", "molecular_id"))
  rownames(gene_info) <- NULL

  return(gene_info)
}

#' @title Map regions for each gene.
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
map_gene_regions <- function(gene_info,
                             region_info,
                             ncore = 1){

  gene_info <- gene_info[gene_info$chrom %in% unique(region_info$chrom), , drop=FALSE]

  if (!all(c("chrom", "p0", "p1") %in% colnames(gene_info))){
    stop("Columns chrom, p0, p1 are required in gene_info!")
  }

  loginfo("Map gene regions...")

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


#' @title Get cross-boundary genes.
#'
#' @param region_info a data frame of region definitions.
#'
#' @param weights a list of preprocessed weights.
#'
#' @param gene_ids a vector of selected gene IDs (z_gene$id).
#' If specified, limits to these genes. Default: use all genes in weights.
#'
#' @param mapping_table a data frame of mapping between molecular traits and genes,
#' with required columns: "molecular_id", "gene_name".
#'
#' @param show_mapping If TRUE, include the mapping between molecular traits and genes
#' in the result.
#'
#' @param ncore The number of cores used to parallelize over genes.
#'
#' @return a data frame of cross-boundary genes.
#'
#' @importFrom logging loginfo
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join select summarise group_by ungroup
#' @importFrom rlang .data
#'
#' @export
get_boundary_genes <- function(region_info,
                               weights,
                               gene_ids = NULL,
                               mapping_table = NULL,
                               show_mapping = TRUE,
                               ncore = 1){

  loginfo("Get boundary genes...")

  # limit to the selected genes
  if (!is.null(gene_ids)) {
    weights <- weights[gene_ids]
  }

  # get gene info from weights
  gene_info <- get_gene_info(weights)
  gene_info <- gene_info[, c("chrom", "p0", "p1", "id", "molecular_id")]

  if (!is.null(mapping_table)) {
    # get combined gene_info
    # map molecular traits to genes, allow for many-to-many matching, return all matches
    # a gene could have molecular traits, and a molecular trait could be linked to multiple genes
    loginfo("Map molecular traits to genes.")

    mapped_gene_info <- gene_info %>%
      left_join(mapping_table[,c("molecular_id", "gene_name")],
                by = "molecular_id",
                multiple = "all", relationship = "many-to-many")

    selected_columns <- c(colnames(gene_info), "gene_name")

    # get molecular trait level positions
    molecular_trait_level_info <- mapped_gene_info %>%
      group_by(.data$id) %>%
      summarise(chrom = .data$chrom[1],
                p0 = min(.data$p0),
                p1 = max(.data$p1),
                molecular_id = .data$molecular_id[1],
                gene_name = paste(unique(.data$gene_name), collapse = ",")) %>%
      ungroup() %>%
      select(all_of(selected_columns))

    # get gene level positions
    gene_level_info <- mapped_gene_info %>%
      group_by(.data$gene_name) %>%
      summarise(chrom = .data$chrom[1],
                p0 = min(.data$p0),
                p1 = max(.data$p1),
                id = paste(unique(.data$id), collapse = ","),
                molecular_id = paste(unique(.data$molecular_id), collapse = ",")) %>%
      ungroup() %>%
      select(all_of(selected_columns))

    if (!show_mapping) {
      # only show the combined molecular traits or genes
      # get molecular trait level positions
      molecular_trait_level_info$gene_name <- ""
      gene_level_info$id <- ""
      gene_level_info$molecular_id <- ""
    }

    # combine gene_info from both molecular trait level and gene level
    gene_info <- rbind(molecular_trait_level_info, gene_level_info)
    gene_info <- unique(gene_info)
    gene_info <- as.data.frame(gene_info)
  }

  # get regions for each molecular trait or gene
  gene_region_info <- map_gene_regions(gene_info,
                                       region_info,
                                       ncore = ncore)

  # get boundary genes (n_regions > 1)
  boundary_genes <- gene_region_info[gene_region_info$n_regions > 1, ]
  boundary_genes <- boundary_genes[with(boundary_genes, order(chrom, p0, p1)), ]
  rownames(boundary_genes) <- NULL

  return(boundary_genes)
}

#' @title Gets boundary genes and selects high PIP boundary genes
#'
#' @param region_info a data frame of region definitions.
#'
#' @param weights a list of preprocessed weights.
#'
#' @param gene_ids a vector of selected gene IDs (z_gene$id).
#' If specified, limits to these genes. Default: use all genes in weights.
#'
#' @param finemap_res a data frame of original finemapping result.
#'
#' @param susie_alpha_res a data frame of original susie alpha result.
#'
#' @param combine_PIPs if TRUE, select boundary genes after combining gene PIPs.
#'
#' @param mapping_table a data frame of mapping between molecular traits and genes,
#' with required columns: "molecular_id", "gene_name".
#'
#' @param pip_thresh PIP cutoff for selecting boundary genes to merge regions.
#'
#' @param filter_cs If TRUE, only select boundary genes in credible sets for region merge.
#'
#' @param ncore The number of cores used to parallelize computation over regions
#'
#' @return a list with boundary genes and selected boundary genes
#'
#' @importFrom logging loginfo
#'
#' @export
#'
select_boundary_genes <- function(region_info,
                                  weights,
                                  gene_ids,
                                  finemap_res,
                                  susie_alpha_res,
                                  combine_PIPs = TRUE,
                                  mapping_table = NULL,
                                  pip_thresh = 0.5,
                                  filter_cs = FALSE,
                                  ncore = 1) {

  # get boundary genes
  # get both molecular trait level and gene level positions when mapping_table is not NULL
  boundary_genes <- get_boundary_genes(region_info,
                                       weights,
                                       gene_ids = gene_ids,
                                       mapping_table = mapping_table,
                                       ncore = ncore)
  loginfo("%d boundary genes.", nrow(boundary_genes))

  # select boundary genes with high PIP and in CS
  finemap_gene_res <- finemap_res[finemap_res$group != "SNP",,drop=FALSE]
  high_PIP_finemap_gene_res <- finemap_gene_res[finemap_gene_res$susie_pip > pip_thresh,,drop=FALSE]

  if (filter_cs) {
    # limit to genes in credible sets
    high_PIP_finemap_gene_res <- high_PIP_finemap_gene_res[!is.na(high_PIP_finemap_gene_res$cs),,drop=FALSE]
  }

  high_PIP_ids <- unique(high_PIP_finemap_gene_res$id)
  selected_boundary_genes <- boundary_genes[which(boundary_genes$id %in% high_PIP_ids), , drop=FALSE]
  loginfo("Selected %d boundary genes with PIP > %s.", nrow(selected_boundary_genes), pip_thresh)

  # select boundary genes with high PIP after combining gene PIPs
  if (combine_PIPs) {
    if (!is.null(mapping_table)) {
      # combine gene-level PIPs
      combined_pip_res <- combine_gene_pips(susie_alpha_res,
                                            mapping_table = mapping_table,
                                            group_by = "gene_name",
                                            filter_cs = filter_cs)
      # select boundary genes with high combined PIP
      high_PIP_combined_pip_res <- combined_pip_res[combined_pip_res$combined_pip > pip_thresh,,drop=FALSE]
      high_PIP_gene_names <- unique(high_PIP_combined_pip_res$gene_name)
      selected_boundary_genes2 <- boundary_genes[which(boundary_genes$gene_name %in% high_PIP_gene_names), , drop=FALSE]
    } else {
      # combine molecular trait-level PIPs
      combined_pip_res <- combine_gene_pips(susie_alpha_res,
                                            group_by = "molecular_id",
                                            filter_cs = filter_cs)
      # select boundary genes with high combined PIP
      high_PIP_combined_pip_res <- combined_pip_res[combined_pip_res$combined_pip > pip_thresh,,drop=FALSE]
      high_PIP_molecular_ids <- unique(high_PIP_combined_pip_res$molecular_id)
      selected_boundary_genes2 <- boundary_genes[which(boundary_genes$molecular_id %in% high_PIP_molecular_ids), , drop=FALSE]
    }

    loginfo("Selected %d boundary genes with combined PIP > %s.", nrow(selected_boundary_genes2), pip_thresh)
    selected_boundary_genes <- unique(rbind(selected_boundary_genes, selected_boundary_genes2))
    loginfo("Selected %d boundary genes in total.", nrow(selected_boundary_genes))
  }

  return(list("boundary_genes" = boundary_genes,
              "selected_boundary_genes" = selected_boundary_genes))
}
