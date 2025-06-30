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
#' @param combine_PIPs if TRUE, select boundary genes after combining gene PIPs.
#'
#' @param mapping_table a data frame of mapping between molecular traits and genes,
#' with required columns: "molecular_id", "gene_name".
#'
#' @param L the number of effects or a vector of number of effects for each region.
#'
#' @param group_prior a vector of prior inclusion probabilities for different groups.
#' If NULL, it will use uniform prior inclusion probabilities.
#'
#' @param group_prior_var a vector of prior variances for different groups.
#' If NULL, it will set prior variance = 50 as the default in \code{susie_rss}.
#'
#' @param pip_thresh PIP cutoff for selecting boundary genes to merge regions.
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
                                      combine_PIPs = TRUE,
                                      mapping_table = NULL,
                                      L = 5,
                                      group_prior = NULL,
                                      group_prior_var = NULL,
                                      pip_thresh = 0.5,
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

  # select boundary genes with high PIP and in CS
  res <- select_boundary_genes(region_info,
                               weights,
                               gene_ids = z_gene$id,
                               finemap_res = finemap_res,
                               susie_alpha_res = susie_alpha_res,
                               combine_PIPs = combine_PIPs,
                               mapping_table = mapping_table,
                               pip_thresh = pip_thresh,
                               filter_cs = filter_cs,
                               ncore = ncore)
  boundary_genes <- res$boundary_genes
  selected_boundary_genes <- res$selected_boundary_genes

  if (nrow(selected_boundary_genes) > 0) {
    # merge region data for selected boundary genes
    loginfo("Merge region data for selected boundary genes...")
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
    loginfo("Rerun fine-mapping for the merged regions...")
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
    loginfo("Update fine-mapping results...")
    res <- update_merged_region_finemap_res(finemap_res,
                                            susie_alpha_res,
                                            merged_region_finemap_res,
                                            merged_region_susie_alpha_res,
                                            merged_region_id_map)
    updated_finemap_res <- res$finemap_res
    updated_susie_alpha_res <- res$susie_alpha_res

    # update region data after region merging
    loginfo("Update region data...")
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
    loginfo("No boundary genes selected. Skipped merging regions.")

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
#' @param combine_PIPs if TRUE, select boundary genes after combining gene PIPs.
#'
#' @param mapping_table a data frame of mapping between molecular traits and genes,
#' with required columns: "molecular_id", "gene_name".
#'
#' @param group_prior a vector of prior inclusion probabilities for different groups.
#' If NULL, it will use uniform prior inclusion probabilities.
#'
#' @param group_prior_var a vector of prior variances for different groups.
#' If NULL, it will set prior variance = 50 as the default in \code{susie_rss}.
#'
#' @param pip_thresh PIP cutoff for selecting boundary genes to merge regions.
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
                                           combine_PIPs = TRUE,
                                           mapping_table = NULL,
                                           group_prior = NULL,
                                           group_prior_var = NULL,
                                           pip_thresh = 0.5,
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

  # select boundary genes with high PIP and in CS
  res <- select_boundary_genes(region_info,
                               weights,
                               gene_ids = z_gene$id,
                               finemap_res = finemap_res,
                               susie_alpha_res = susie_alpha_res,
                               combine_PIPs = combine_PIPs,
                               mapping_table = mapping_table,
                               pip_thresh = pip_thresh,
                               filter_cs = filter_cs,
                               ncore = ncore)
  boundary_genes <- res$boundary_genes
  selected_boundary_genes <- res$selected_boundary_genes

  if (nrow(selected_boundary_genes) > 0) {
    # merge region data for selected boundary genes
    loginfo("Merge region data for selected boundary genes...")
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
    loginfo("Rerun fine-mapping for the merged regions...")
    res <- finemap_regions_noLD(merged_region_data,
                                group_prior = group_prior,
                                group_prior_var = group_prior_var,
                                ncore = ncore,
                                verbose = verbose,
                                ...)
    merged_region_finemap_res <- res$finemap_res
    merged_region_susie_alpha_res <- res$susie_alpha_res

    # update fine-mapping results after region merging
    loginfo("Update fine-mapping results...")
    res <- update_merged_region_finemap_res(finemap_res,
                                            susie_alpha_res,
                                            merged_region_finemap_res,
                                            merged_region_susie_alpha_res,
                                            merged_region_id_map)
    updated_finemap_res <- res$finemap_res
    updated_susie_alpha_res <- res$susie_alpha_res

    # update region data after region merging
    loginfo("Update region data...")
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
    loginfo("No boundary genes selected. Skipped merging regions.")

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

