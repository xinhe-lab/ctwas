#' Causal inference for TWAS using summary statistics
#'
#' @param z_snp A data frame with four columns: "id", "A1", "A2", "z".
#' giving the z scores for snps. "A1" is effect allele. "A2" is the other allele.
#'
#' @param weights a list of weights
#'
#' @param region_info a data frame of region definition and associated LD file names
#'
#' @param niter_prefit the number of iterations of the E-M algorithm to perform during the initial parameter estimation step
#'
#' @param niter the number of iterations of the E-M algorithm to perform during the complete parameter estimation step
#'
#' @param thin The proportion of SNPs to be used for estimating parameters and screening regions.
#'
#' @param L the number of effects for susie during the fine mapping steps
#'
#' @param init_group_prior a vector of initial values of prior inclusion probabilities for SNPs and genes.
#'
#' @param init_group_prior_var a vector of initial values of prior variances for SNPs and gene effects.
#'
#' @param group_prior_var_structure a string indicating the structure to put on the prior variance parameters.
#' "shared_type" allows all groups in one molecular QTL type to share the same variance parameter.
#' "shared_context" allows all groups in one molecular QTL type to share the same variance parameter.
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
#' @param coverage A number between 0 and 1 specifying the \dQuote{coverage} of the estimated confidence sets
#'
#' @param min_abs_corr Minimum absolute correlation allowed in a
#'.  credible set. The default, 0.5, corresponds to a squared
#'   correlation of 0.25, which is a commonly used threshold for
#'   genotype data in genetic studies.
#'
#' @param ncore The number of cores used to parallelize susie over regions
#'
#' @param save_cor TRUE/FALSE. If TRUE, save correlation (R) matrices to \code{outputdir}
#'
#' @param cor_dir a string, the directory to store correlation (R) matrices
#'
#' @param logfile the log file, if NULL will print log info on screen
#'
#' @param verbose TRUE/FALSE. If TRUE, print detailed messages
#'
#' @importFrom logging addHandler loginfo writeToFile
#'
#' @return a list of estimated parameters, fine-mapping results, boundary genes, z_gene, and region_data
#'
#' @export
#'
ctwas_sumstats <- function(
    z_snp,
    weights,
    region_info,
    thin = 0.1,
    niter_prefit = 3,
    niter = 30,
    L = 5,
    init_group_prior = NULL,
    init_group_prior_var = NULL,
    group_prior_var_structure = c("shared_type", "shared_context", "shared_nonSNP", "shared_all", "independent"),
    maxSNP = 20000,
    min_nonSNP_PIP = 0.5,
    p_single_effect = 0.8,
    use_null_weight = TRUE,
    coverage = 0.95,
    min_abs_corr = 0.5,
    ncore = 1,
    save_cor = FALSE,
    cor_dir = getwd(),
    logfile = NULL,
    verbose = FALSE,
    ...){

  if (!is.null(logfile)){
    addHandler(writeToFile, file=logfile, level='DEBUG')
  }

  # Compute gene z-scores
  z_gene <- compute_gene_z(z_snp, weights, ncore = ncore)

  # Get region_data, which contains SNPs and genes assigned to each region
  #. downsample SNPs if thin < 1
  #. assign SNP and gene IDs, and z-scores to each region
  #. find boundary genes and adjust region_data for boundary genes
  res <- assemble_region_data(region_info,
                              z_snp,
                              z_gene,
                              weights,
                              thin = thin,
                              maxSNP = maxSNP,
                              trim_by = "random",
                              adjust_boundary_genes = TRUE,
                              ncore = ncore)
  region_data <- res$region_data
  boundary_genes <- res$boundary_genes

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
                     ncore = ncore)
  group_prior <- param$group_prior
  group_prior_var <- param$group_prior_var

  # Screen regions
  #. fine-map all regions with thinned SNPs
  #. select regions with strong non-SNP signals
  region_nonSNP_PIP_df <- screen_regions(region_data,
                                         region_info,
                                         weights,
                                         group_prior = group_prior,
                                         group_prior_var = group_prior_var,
                                         L = L,
                                         minvar = 2,
                                         mingene = 1,
                                         min_nonSNP_PIP = min_nonSNP_PIP,
                                         ncore = ncore,
                                         verbose = verbose)
  screened_region_ids <- region_nonSNP_PIP_df$region_id
  screened_region_data <- region_data[screened_region_ids]

  # Expand screened region_data with all SNPs in the regions
  if (thin < 1){
    screened_region_data <- expand_region_data(screened_region_data,
                                               region_info,
                                               z_snp,
                                               z_gene,
                                               trim_by = "z",
                                               maxSNP = maxSNP,
                                               ncore = ncore)
    # update data in screened regions with screened_region_data (full SNPs)
    region_data[names(screened_region_data)] <- screened_region_data
  }

  # Run fine-mapping for regions with strong gene signals using full SNPs
  #. save correlation matrices if save_cor is TRUE
  if (length(screened_region_data) > 0){
    finemap_res <- finemap_regions(screened_region_data,
                                   region_info,
                                   weights,
                                   group_prior = group_prior,
                                   group_prior_var = group_prior_var,
                                   L = L,
                                   use_null_weight = use_null_weight,
                                   coverage = coverage,
                                   min_abs_corr = min_abs_corr,
                                   save_cor = save_cor,
                                   cor_dir = cor_dir,
                                   ncore = ncore,
                                   verbose = verbose)
  } else {
    finemap_res <- NULL
  }

  return(list("param" = param,
              "finemap_res" = finemap_res,
              "boundary_genes" = boundary_genes,
              "z_gene" = z_gene,
              "region_data" = region_data))

}

