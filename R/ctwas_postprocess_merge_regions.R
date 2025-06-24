#' @title Runs cTWAS post-processing procedure for merging regions
#'
#' @param region_info a data frame of region definitions.
#'
#' @param region_data region_data to be fine-mapped.
#'
#' @param z_snp A data frame with columns: "id", "z", giving the z-scores for SNPs.
#'
#' @param z_gene A data frame with columns: "id", "z", giving the z-scores for genes.
#'
#' @param weights a list of preprocessed weights.
#'
#' @param LD_map a data frame with filenames of LD matrices for the regions.
#'
#' @param snp_map a list of data frames with SNP-to-region map for the reference.
#'
#' @param finemap_res a data frame of original finemapping result.
#'
#' @param susie_alpha_res a data frame of original susie alpha result.
#'
#' @param L the number of effects or a vector of number of effects for each region.
#'
#' @param group_prior a vector of prior inclusion probabilities for different groups.
#' If NULL, it will use uniform prior inclusion probabilities.
#'
#' @param group_prior_var a vector of prior variances for different groups.
#' If NULL, it will set prior variance = 50 as the default in \code{susie_rss}.
#'
#' @param combined_pip_res a data frame of combined gene PIP result from \code{combine_gene_pips()}.
#'
#' @param mapping_table a data frame of mapping between molecular traits and genes,
#' with required columns: "molecular_id", "gene_name".
#'
#' @param show_mapping If TRUE, include the mapping between molecular traits and genes
#' in the result.
#'
#' @param min_PIP PIP cutoff for selecting boundary genes to merge regions.
#'
#' @param filter_cs If TRUE, only select boundary genes in credible sets for region merge.
#'
#' @param maxSNP Inf or integer. Maximum number of SNPs in a region. Default is
#' Inf, no limit. This can be useful if there are many SNPs in a region and you don't
#' have enough memory to run the program.
#'
#' @param save_cor If TRUE, save correlation (R) matrices to \code{cor_dir}
#'
#' @param cor_dir a string, the directory to store correlation (R) matrices
#'
#' @param ncore The number of cores used to parallelize computation over regions
#'
#' @param verbose If TRUE, print detail messages.

#' @param logfile the log file, if NULL will print log info on screen
#'
#' @param ... Additional arguments of \code{finemap_regions}.
#'
#' @return a list with region merge results.
#'
#' @importFrom logging addHandler loginfo logwarn writeToFile
#'
#' @export
#'
postprocess_merge_regions <- function(region_info,
                                      region_data,
                                      z_snp,
                                      z_gene,
                                      weights,
                                      LD_map,
                                      snp_map,
                                      finemap_res,
                                      susie_alpha_res,
                                      L = 5,
                                      group_prior = NULL,
                                      group_prior_var = NULL,
                                      combined_pip_res = NULL,
                                      mapping_table = NULL,
                                      show_mapping = FALSE,
                                      min_PIP = 0.5,
                                      filter_cs = FALSE,
                                      maxSNP = Inf,
                                      save_cor = FALSE,
                                      cor_dir = NULL,
                                      ncore = 1,
                                      verbose = FALSE,
                                      logfile = NULL,
                                      ...){

  if (!is.null(logfile)) {
    addHandler(writeToFile, file=logfile, level='DEBUG')
  }

  loginfo("Run post-processing procedure for merging regions...")

  # get boundary genes (combining both molecular traits and genes)
  boundary_genes <- get_boundary_genes(region_info,
                                       weights,
                                       gene_ids = z_gene$id,
                                       mapping_table = mapping_table,
                                       ncore = ncore)

  loginfo("%d boundary genes (molecular traits).", nrow(boundary_genes))

  # select boundary genes with PIP > 0.5 and in CS
  finemap_gene_res <- finemap_res[finemap_res$group != "SNP",,drop=FALSE]
  high_PIP_finemap_gene_res <- finemap_gene_res[finemap_gene_res$susie_pip > min_PIP,,drop=FALSE]

  if (filter_cs) {
    # limit to genes in credible sets
    high_PIP_finemap_gene_res <- high_PIP_finemap_gene_res[!is.na(high_PIP_finemap_gene_res$cs),,drop=FALSE]
  }

  high_PIP_ids <- unique(high_PIP_finemap_gene_res$id)
  selected_boundary_genes <- boundary_genes[which(boundary_genes$id %in% high_PIP_ids), , drop=FALSE]
  if (verbose)
    loginfo("%d boundary molecular traits selected.", nrow(selected_boundary_genes))

  if (!is.null(combined_pip_res)) {
    if (is.null(boundary_genes$gene_name))
      stop("mapping_table is required when using combined_pip_res!")
    high_PIP_combined_pip_res <- combined_pip_res[combined_pip_res$combined_pip > min_PIP,,drop=FALSE]
    high_PIP_gene_names <- unique(high_PIP_combined_pip_res$gene_name)
    selected_boundary_genes2 <- boundary_genes[which(boundary_genes$gene_name %in% high_PIP_gene_names), , drop=FALSE]
    if (verbose)
      loginfo("%d boundary genes selected.", nrow(selected_boundary_genes2))
    selected_boundary_genes <- unique(rbind(selected_boundary_genes, selected_boundary_genes2))
  }
  loginfo("%d boundary genes (molecular traits) selected in total.", nrow(selected_boundary_genes))

  if (nrow(selected_boundary_genes) > 0) {
    # merge region data for selected boundary genes
    res <- merge_region_data(selected_boundary_genes,
                             region_data = region_data,
                             region_info = region_info,
                             LD_map = LD_map,
                             snp_map = snp_map,
                             z_snp = z_snp,
                             z_gene = z_gene,
                             maxSNP = maxSNP,
                             ncore = ncore)
    merged_region_data <- res$merged_region_data
    merged_region_info <- res$merged_region_info
    merged_LD_map <- res$merged_LD_map
    merged_snp_map <- res$merged_snp_map
    merged_region_id_map <- res$merged_region_id_map

    # rerun fine-mapping for the merged regions
    res <- finemap_regions(merged_region_data,
                           LD_map = merged_LD_map,
                           weights = weights,
                           group_prior = group_prior,
                           group_prior_var = group_prior_var,
                           L = L,
                           save_cor = save_cor,
                           cor_dir = cor_dir,
                           ncore = ncore,
                           verbose = verbose,
                           ...)
    merged_region_finemap_res <- res$finemap_res
    merged_region_susie_alpha_res <- res$susie_alpha_res

    # update fine-mapping results after region merging
    res <- update_merged_region_finemap_res(finemap_res,
                                            susie_alpha_res,
                                            merged_region_finemap_res,
                                            merged_region_susie_alpha_res,
                                            merged_region_id_map)
    updated_finemap_res <- res$finemap_res
    updated_susie_alpha_res <- res$susie_alpha_res

    # update region data after region merging
    res <- update_merged_region_data(region_data, merged_region_data,
                                     region_info, merged_region_info,
                                     LD_map, merged_LD_map,
                                     snp_map, merged_snp_map,
                                     merged_region_id_map)
    updated_region_data <- res$updated_region_data
    updated_region_info <- res$updated_region_info
    updated_LD_map <- res$updated_LD_map
    updated_snp_map <- res$updated_snp_map

    rm(res)

    return(list("boundary_genes" = boundary_genes,
                "selected_boundary_genes" = selected_boundary_genes,
                "merged_region_data" = merged_region_data,
                "merged_region_info" = merged_region_info,
                "merged_LD_map" = merged_LD_map,
                "merged_snp_map" = merged_snp_map,
                "merged_region_id_map" = merged_region_id_map,
                "merged_region_finemap_res" = merged_region_finemap_res,
                "merged_region_susie_alpha_res" = merged_region_susie_alpha_res,
                "updated_region_data" = updated_region_data,
                "updated_region_info" = updated_region_info,
                "updated_LD_map" = updated_LD_map,
                "updated_snp_map" = updated_snp_map,
                "updated_finemap_res" = updated_finemap_res,
                "updated_susie_alpha_res" = updated_susie_alpha_res))
  } else {
    loginfo("No boundary genes selected. Skip merging regions.")

    return(list("boundary_genes" = boundary_genes,
                "selected_boundary_genes" = selected_boundary_genes,
                "merged_region_data" = NULL,
                "merged_region_info" = NULL,
                "merged_LD_map" = NULL,
                "merged_snp_map" = NULL,
                "merged_region_id_map" = NULL,
                "merged_region_finemap_res" = NULL,
                "merged_region_susie_alpha_res" = NULL,
                "updated_region_data" = region_data,
                "updated_region_info" = region_info,
                "updated_LD_map" = LD_map,
                "updated_snp_map" = snp_map,
                "updated_finemap_res" = finemap_res,
                "updated_susie_alpha_res" = susie_alpha_res))
  }

}

#' @title Runs cTWAS post-processing procedure for merging regions without LD
#'
#' @param region_info a data frame of region definitions.
#'
#' @param region_data region_data to be fine-mapped.
#'
#' @param z_snp A data frame with columns: "id", "z", giving the z-scores for SNPs.
#'
#' @param z_gene A data frame with columns: "id", "z", giving the z-scores for genes.
#'
#' @param weights a list of preprocessed weights.
#'
#' @param snp_map a list of data frames with SNP-to-region map for the reference.
#'
#' @param finemap_res a data frame of original finemapping result.
#'
#' @param susie_alpha_res a data frame of original susie alpha result.
#'
#' @param group_prior a vector of prior inclusion probabilities for different groups.
#' If NULL, it will use uniform prior inclusion probabilities.
#'
#' @param group_prior_var a vector of prior variances for different groups.
#' If NULL, it will set prior variance = 50 as the default in \code{susie_rss}.
#'
#' @param combined_pip_res a data frame of combined gene PIP result from \code{combine_gene_pips()}.
#'
#' @param mapping_table a data frame of mapping between molecular traits and genes,
#' with required columns: "molecular_id", "gene_name".
#'
#' @param show_mapping If TRUE, include the mapping between molecular traits and genes
#' in the result.
#'
#' @param min_PIP PIP cutoff for selecting boundary genes to merge regions.
#'
#' @param filter_cs If TRUE, only select boundary genes in credible sets for region merge.
#'
#' @param maxSNP Inf or integer. Maximum number of SNPs in a region. Default is
#' Inf, no limit. This can be useful if there are many SNPs in a region and you don't
#' have enough memory to run the program.
#'
#' @param ncore The number of cores used to parallelize computation over regions
#'
#' @param verbose If TRUE, print detail messages.

#' @param logfile the log file, if NULL will print log info on screen
#'
#' @param ... Additional arguments of \code{finemap_regions}.
#'
#' @return a list with region merge results.
#'
#' @importFrom logging addHandler loginfo logwarn writeToFile
#'
#' @export
#'
postprocess_merge_regions_noLD <- function(region_info,
                                           region_data,
                                           z_snp,
                                           z_gene,
                                           weights,
                                           snp_map,
                                           finemap_res,
                                           susie_alpha_res,
                                           group_prior = NULL,
                                           group_prior_var = NULL,
                                           combined_pip_res = NULL,
                                           mapping_table = NULL,
                                           show_mapping = FALSE,
                                           min_PIP = 0.5,
                                           filter_cs = FALSE,
                                           maxSNP = Inf,
                                           ncore = 1,
                                           verbose = FALSE,
                                           logfile = NULL,
                                           ...){

  if (!is.null(logfile)) {
    addHandler(writeToFile, file=logfile, level='DEBUG')
  }

  loginfo("Run post-processing procedure for merging regions without LD...")

  # get boundary genes (combining both molecular traits and genes)
  boundary_genes <- get_boundary_genes(region_info,
                                       weights,
                                       gene_ids = z_gene$id,
                                       mapping_table = mapping_table,
                                       ncore = ncore)

  loginfo("%d boundary genes (molecular traits).", nrow(boundary_genes))

  # select boundary genes with PIP > 0.5 and in CS
  finemap_gene_res <- finemap_res[finemap_res$group != "SNP",,drop=FALSE]
  high_PIP_finemap_gene_res <- finemap_gene_res[finemap_gene_res$susie_pip > min_PIP,,drop=FALSE]

  if (filter_cs) {
    # limit to genes in credible sets
    high_PIP_finemap_gene_res <- high_PIP_finemap_gene_res[!is.na(high_PIP_finemap_gene_res$cs),,drop=FALSE]
  }

  high_PIP_ids <- unique(high_PIP_finemap_gene_res$id)
  selected_boundary_genes <- boundary_genes[which(boundary_genes$id %in% high_PIP_ids), , drop=FALSE]
  if (verbose)
    loginfo("%d boundary molecular traits selected.", nrow(selected_boundary_genes))

  if (!is.null(combined_pip_res)) {
    if (is.null(boundary_genes$gene_name))
      stop("mapping_table is required when using combined_pip_res!")
    high_PIP_combined_pip_res <- combined_pip_res[combined_pip_res$combined_pip > min_PIP,,drop=FALSE]
    high_PIP_gene_names <- unique(high_PIP_combined_pip_res$gene_name)
    selected_boundary_genes2 <- boundary_genes[which(boundary_genes$gene_name %in% high_PIP_gene_names), , drop=FALSE]
    if (verbose)
      loginfo("%d boundary genes selected.", nrow(selected_boundary_genes2))
    selected_boundary_genes <- unique(rbind(selected_boundary_genes, selected_boundary_genes2))
  }
  loginfo("%d boundary genes (molecular traits) selected in total.", nrow(selected_boundary_genes))

  if (nrow(selected_boundary_genes) > 0) {
    # merge region data for selected boundary genes
    res <- merge_region_data_noLD(selected_boundary_genes,
                                  region_data = region_data,
                                  region_info = region_info,
                                  snp_map = snp_map,
                                  z_snp = z_snp,
                                  z_gene = z_gene,
                                  maxSNP = maxSNP,
                                  ncore = ncore)
    merged_region_data <- res$merged_region_data
    merged_region_info <- res$merged_region_info
    merged_snp_map <- res$merged_snp_map
    merged_region_id_map <- res$merged_region_id_map

    # rerun fine-mapping for the merged regions
    res <- finemap_regions_noLD(merged_region_data,
                                group_prior = group_prior,
                                group_prior_var = group_prior_var,
                                ncore = ncore,
                                verbose = verbose,
                                ...)
    merged_region_finemap_res <- res$finemap_res
    merged_region_susie_alpha_res <- res$susie_alpha_res

    # update fine-mapping results after region merging
    res <- update_merged_region_finemap_res(finemap_res,
                                            susie_alpha_res,
                                            merged_region_finemap_res,
                                            merged_region_susie_alpha_res,
                                            merged_region_id_map)
    updated_finemap_res <- res$finemap_res
    updated_susie_alpha_res <- res$susie_alpha_res

    # update region data after region merging
    res <- update_merged_region_data_noLD(region_data, merged_region_data,
                                          region_info, merged_region_info,
                                          snp_map, merged_snp_map,
                                          merged_region_id_map)
    updated_region_data <- res$updated_region_data
    updated_region_info <- res$updated_region_info
    updated_snp_map <- res$updated_snp_map

    rm(res)

    return(list("boundary_genes" = boundary_genes,
                "selected_boundary_genes" = selected_boundary_genes,
                "merged_region_data" = merged_region_data,
                "merged_region_info" = merged_region_info,
                "merged_snp_map" = merged_snp_map,
                "merged_region_id_map" = merged_region_id_map,
                "merged_region_finemap_res" = merged_region_finemap_res,
                "merged_region_susie_alpha_res" = merged_region_susie_alpha_res,
                "updated_region_data" = updated_region_data,
                "updated_region_info" = updated_region_info,
                "updated_snp_map" = updated_snp_map,
                "updated_finemap_res" = updated_finemap_res,
                "updated_susie_alpha_res" = updated_susie_alpha_res))
  } else {
    loginfo("No boundary genes selected. Skip merging regions.")

    return(list("boundary_genes" = boundary_genes,
                "selected_boundary_genes" = selected_boundary_genes,
                "merged_region_data" = NULL,
                "merged_region_info" = NULL,
                "merged_snp_map" = NULL,
                "merged_region_id_map" = NULL,
                "merged_region_finemap_res" = NULL,
                "merged_region_susie_alpha_res" = NULL,
                "updated_region_data" = region_data,
                "updated_region_info" = region_info,
                "updated_snp_map" = snp_map,
                "updated_finemap_res" = finemap_res,
                "updated_susie_alpha_res" = susie_alpha_res))
  }

}
