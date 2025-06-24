#' @title Runs cTWAS post-processing procedure for merging regions
#'
#' @param region_data region_data to be fine-mapped.
#'
#' @param z_snp A data frame with columns: "id", "z", giving the z-scores for SNPs.
#'
#' @param weights a list of preprocessed weights.
#'
#' @param LD_map a data frame with filenames of LD matrices for the regions.
#'
#' @param snp_map a list of data frames with SNP-to-region map for the reference.
#'
#' @param gwas_n integer, GWAS sample size.
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
#' @param min_nonSNP_PIP regions with total non-SNP PIP >= \code{min_nonSNP_PIP}
#' will be selected to run LD mismatch diagnosis.
#'
#' @param filter_cs If TRUE, computes region non-SNP PIP only in credible sets.
#'
#' @param p_diff_thresh p-value cutoff for identifying problematic SNPs
#' with significant difference between observed z-scores and estimated values.
#'
#' @param pip_thresh cutoff of gene PIP to select problematic genes.
#'
#' @param z_thresh cutoff of abs(z-scores) to select problematic genes.
#'
#' @param plot If TRUE, plot observed z-score vs the expected value.
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
postprocess_LD_mismatch <- function(region_data,
                                    z_snp,
                                    weights,
                                    LD_map,
                                    snp_map,
                                    gwas_n,
                                    finemap_res,
                                    susie_alpha_res,
                                    group_prior = NULL,
                                    group_prior_var = NULL,
                                    min_nonSNP_PIP = 0.8,
                                    filter_cs = FALSE,
                                    p_diff_thresh = 5e-6,
                                    pip_thresh = 0.5,
                                    z_thresh = NULL,
                                    plot = TRUE,
                                    maxSNP = Inf,
                                    ncore = 1,
                                    verbose = FALSE,
                                    logfile = NULL,
                                    ...){

  if (!is.null(logfile)) {
    addHandler(writeToFile, file=logfile, level='DEBUG')
  }

  loginfo("Run post-processing procedure for LD mismatch diagnosis...")

  # select regions
  nonSNP_PIPs <- compute_region_nonSNP_PIPs(finemap_res, filter_cs = filter_cs)
  selected_region_ids <- names(nonSNP_PIPs)[nonSNP_PIPs > min_nonSNP_PIP]
  loginfo("Selected %d regions with non-SNP PIP > %s.", length(selected_region_ids), min_nonSNP_PIP)

  if (length(selected_region_ids) > 0) {
    diagnose_LD_mismatch_res <- diagnose_LD_mismatch_susie(region_ids = selected_region_ids,
                                                            z_snp = z_snp,
                                                            LD_map = LD_map,
                                                            gwas_n = gwas_n,
                                                            p_diff_thresh = p_diff_thresh,
                                                            plot = plot,
                                                            ncore = ncore)
    problematic_snps <- diagnose_LD_mismatch_res$problematic_snps
    flipped_snps <- diagnose_LD_mismatch_res$flipped_snps
    condz_stats <- diagnose_LD_mismatch_res$condz_stats
    plots <- diagnose_LD_mismatch_res$plots
  } else {
    diagnose_LD_mismatch_res <- NULL
    problematic_snps <- NULL
    flipped_snps <- NULL
    condz_stats <- NULL
    plots <- NULL
  }

  # get problematic genes
  if (length(problematic_snps) > 0) {
    problematic_genes <- get_problematic_genes(problematic_snps,
                                               weights,
                                               finemap_res,
                                               pip_thresh = pip_thresh,
                                               z_thresh = z_thresh)
  } else {
    loginfo('No problematic SNPs. No need to update fine-mapping results.')
    problematic_genes <- NULL
  }

  # rerun the fine-mapping without LD for the regions with problematic genes
  if (length(problematic_genes) > 0) {
    problematic_region_ids <- unique(finemap_res$region_id[finemap_res$id %in% problematic_genes])
    rerun_region_data <- region_data[problematic_region_ids]
    thin <- min(sapply(rerun_region_data, "[[", "thin"))
    # expand regions with all SNPs before fine-mapping
    if (thin < 1) {
      loginfo("thin != 1 for region_data, expand regions to include all SNPs.")
      rerun_region_data <- expand_region_data(rerun_region_data,
                                              snp_map,
                                              z_snp,
                                              maxSNP = maxSNP,
                                              ncore = ncore)
    }

    # rerun the fine-mapping without LD
    res <- finemap_regions_noLD(rerun_region_data,
                                group_prior = group_prior,
                                group_prior_var = group_prior_var,
                                ncore = ncore,
                                verbose = verbose,
                                ...)
    rerun_finemap_res <- res$finemap_res
    rerun_susie_alpha_res <- res$susie_alpha_res

    # update the fine-mapping results for the problematic regions
    res <- update_finemap_res(finemap_res,
                              susie_alpha_res,
                              rerun_finemap_res,
                              rerun_susie_alpha_res,
                              updated_region_ids = problematic_region_ids)
    updated_finemap_res <- res$finemap_res
    updated_susie_alpha_res <- res$susie_alpha_res
    rm(res)
  } else {
    problematic_region_ids <- NULL
    rerun_finemap_res <- NULL
    rerun_susie_alpha_res <- NULL
    updated_finemap_res <- finemap_res
    updated_susie_alpha_res <- susie_alpha_res
  }

  return(list("selected_region_ids" = selected_region_ids,
              "problematic_snps" = problematic_snps,
              "flipped_snps" = flipped_snps,
              "condz_stats" = condz_stats,
              "plots" = plots,
              "problematic_genes" = problematic_genes,
              "problematic_region_ids" = problematic_region_ids,
              "rerun_finemap_res" = rerun_finemap_res,
              "rerun_susie_alpha_res" = rerun_susie_alpha_res,
              "updated_finemap_res" = updated_finemap_res,
              "updated_susie_alpha_res" = updated_susie_alpha_res))
}
