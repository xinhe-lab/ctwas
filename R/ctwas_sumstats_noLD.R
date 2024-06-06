#' @title Runs cTWAS analysis with "no LD" version
#'
#' @param z_snp A data frame with four columns: "id", "A1", "A2", "z".
#' giving the z scores for snps. "A1" is effect allele. "A2" is the other allele.
#'
#' @param weights a list of prediction weights
#'
#' @param region_info a data frame of region definition and associated LD file names
#'
#' @param snp_info a list or data frame of reference SNP info.
#'
#' @param z_gene A data frame with columns: "id", "z", giving the z-scores for genes.
#'
#' @param niter_prefit the number of iterations of the E-M algorithm to perform during the initial parameter estimation step
#'
#' @param niter the number of iterations of the E-M algorithm to perform during the complete parameter estimation step
#'
#' @param thin The proportion of SNPs to be used for estimating parameters and screening regions.
#'
#' @param init_group_prior a vector of initial values of prior inclusion probabilities for SNPs and genes.
#'
#' @param init_group_prior_var a vector of initial values of prior variances for SNPs and gene effects.
#'
#' @param group_prior_var_structure a string indicating the structure to put on the prior variance parameters.
#' "shared_type" allows all groups in one molecular QTL type to share the same variance parameter.
#' "shared_context" allows all groups in one context (tissue, cell type, condition) to share the same variance parameter.
#' "shared_nonSNP" allows all non-SNP groups to share the same variance parameter.
#' "shared_all" allows all groups to share the same variance parameter.
#' "independent" allows all groups to have their own separate variance parameters.
#'
#' @param maxSNP Inf or integer. Maximum number of SNPs in a region. Default is
#' Inf, no limit. This can be useful if there are many SNPs in a region and you don't
#' have enough memory to run the program. This applies to the finemapping step only.
#'
#' @param min_nonSNP_PIP Regions with non-SNP PIP >= \code{min_nonSNP_PIP}
#' will be selected to run finemapping using full SNPs.
#'
#' @param p_single_effect Regions with probability greater than \code{p_single_effect} of
#' having 1 or fewer effects will be used for parameter estimation
#'
#' @param use_null_weight TRUE/FALSE. If TRUE, allow for a probability of no effect in susie
#'
#' @param ncore The number of cores used to parallelize susie over regions
#'
#' @param outputdir The directory to store output. If specified, save outputs to the directory.
#'
#' @param outname The output name.
#'
#' @param logfile path to the log file, if NULL will print log info on screen.
#'
#' @param verbose TRUE/FALSE. If TRUE, print detailed messages
#'
#' @param ... Additional arguments of \code{susie_rss}.
#'
#' @importFrom logging addHandler loginfo writeToFile
#'
#' @return a list of estimated parameters, fine-mapping results, boundary genes,
#' z_gene, region_data, and screened_region_data
#'
#' @export
#'
ctwas_sumstats_noLD <- function(
    z_snp,
    weights,
    region_info,
    snp_info,
    z_gene = NULL,
    thin = 0.1,
    niter_prefit = 3,
    niter = 30,
    init_group_prior = NULL,
    init_group_prior_var = NULL,
    group_prior_var_structure = c("shared_type", "shared_context", "shared_nonSNP", "shared_all", "independent"),
    maxSNP = 20000,
    min_nonSNP_PIP = 0.5,
    p_single_effect = 0.8,
    use_null_weight = TRUE,
    ncore = 1,
    outputdir = NULL,
    outname = "ctwas_noLD",
    logfile = NULL,
    verbose = FALSE){

  if (!is.null(logfile)) {
    addHandler(writeToFile, file=logfile, level='DEBUG')
  }

  if (!is.null(outputdir)) {
    dir.create(outputdir, showWarnings=FALSE, recursive=TRUE)
  }

  # Compute gene z-scores
  if (is.null(z_gene)) {
    z_gene <- compute_gene_z(z_snp, weights, ncore = ncore)
    if (!is.null(outputdir)) {
      saveRDS(z_gene, file.path(outputdir, paste0(outname, ".z_gene.RDS")))
    }
  }

  # Get region_data, which contains SNPs and genes assigned to each region
  #. downsample SNPs if thin < 1
  #. assign SNP and gene IDs, and z-scores to each region
  #. find boundary genes and adjust region_data for boundary genes
  region_data_res <- assemble_region_data(region_info,
                                          z_snp,
                                          z_gene,
                                          weights,
                                          snp_info,
                                          thin = thin,
                                          maxSNP = maxSNP,
                                          trim_by = "random",
                                          adjust_boundary_genes = TRUE,
                                          ncore = ncore)
  region_data <- region_data_res$region_data
  boundary_genes <- region_data_res$boundary_genes
  if (!is.null(outputdir)) {
    saveRDS(region_data, file.path(outputdir, paste0(outname, ".region_data.thin", thin, ".RDS")))
    saveRDS(boundary_genes, file.path(outputdir, paste0(outname, ".boundary_genes.RDS")))
  }

  # Estimate parameters
  #. get region_data for all the regions
  #. run EM for two rounds with thinned SNPs using L = 1
  param <- est_param(region_data,
                     init_group_prior = init_group_prior,
                     init_group_prior_var = init_group_prior_var,
                     group_prior_var_structure = group_prior_var_structure,
                     niter_prefit = niter_prefit,
                     niter = niter,
                     p_single_effect = p_single_effect,
                     ncore = ncore,
                     verbose = verbose)

  group_prior <- param$group_prior
  group_prior_var <- param$group_prior_var
  if (!is.null(outputdir)) {
    saveRDS(param, file.path(outputdir, paste0(outname, ".param.RDS")))
  }

  # Screen regions
  #. fine-map all regions with thinned SNPs
  #. select regions with strong non-SNP signals
  screen_regions_res <- screen_regions(region_data,
                                       use_LD = FALSE,
                                       group_prior = group_prior,
                                       group_prior_var = group_prior_var,
                                       L = 1,
                                       min_nonSNP_PIP = min_nonSNP_PIP,
                                       ncore = ncore,
                                       verbose = verbose)
  screened_region_data <- screen_regions_res$screened_region_data
  region_nonSNP_PIP_df <- screen_regions_res$region_nonSNP_PIP_df

  # Expand screened region_data with all SNPs in the regions
  if (thin < 1){
    screened_region_data <- expand_region_data(screened_region_data,
                                               snp_info,
                                               z_snp,
                                               z_gene,
                                               trim_by = "z",
                                               maxSNP = maxSNP,
                                               ncore = ncore)
    # # update data in screened regions with screened_region_data (full SNPs)
    # region_data[screened_region_ids] <- screened_region_data
    # if (!is.null(outputdir)) {
    #   saveRDS(region_data, file.path(outputdir, paste0(outname, ".region_data.RDS")))
    # }
  }
  if (!is.null(outputdir)) {
    saveRDS(screened_region_data, file.path(outputdir, paste0(outname, ".screened_region_data.RDS")))
  }

  # Run fine-mapping for regions with strong gene signals using full SNPs
  #. save correlation matrices if save_cor is TRUE
  if (length(screened_region_data) > 0){
    finemap_res <- finemap_regions(screened_region_data,
                                   use_LD = FALSE,
                                   group_prior = group_prior,
                                   group_prior_var = group_prior_var,
                                   L = 1,
                                   use_null_weight = use_null_weight,
                                   include_cs_index = FALSE,
                                   ncore = ncore,
                                   verbose = verbose)
    if (!is.null(outputdir)) {
      saveRDS(finemap_res, file.path(outputdir, paste0(outname, ".finemap_res.RDS")))
    }
  } else {
    warning("No regions selected for finemapping.")
    finemap_res <- NULL
  }

  return(list("param" = param,
              "finemap_res" = finemap_res,
              "boundary_genes" = boundary_genes,
              "z_gene" = z_gene,
              "region_data" = region_data,
              "screened_region_data" = screened_region_data))

}

