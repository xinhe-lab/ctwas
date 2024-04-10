#' Causal inference for TWAS using summary statistics
#'
#' @param z_snp A data frame with four columns: "id", "A1", "A2", "z".
#' giving the z scores for snps. "A1" is effect allele. "A2" is the other allele.
#'
#' @param weights a list of weights
#'
#' @param region_info a data frame of region definition and associated LD file names
#'
#' @param niter1 the number of iterations of the E-M algorithm to perform during the initial parameter estimation step
#'
#' @param niter2 the number of iterations of the E-M algorithm to perform during the complete parameter estimation step
#'
#' @param thin The proportion of SNPs to be used for parameter estimation and initial screening regions.
#' Smaller \code{thin} parameters reduce runtime at the expense of accuracy.
#' The fine mapping step is rerun using full SNPs
#' for regions with strong gene signals; see \code{min_nonSNP_PIP}.
#'
#' @param L the number of effects for susie during the fine mapping steps
#'
#' @param init_group_prior a vector of initial values of prior inclusion probabilities for SNPs and genes.
#'
#' @param init_group_prior_var a vector of initial values of prior variances for SNPs and gene effects.
#'
#' @param group_prior_var_structure a string indicating the structure to put on the prior variance parameters.
#' "independent" is the default and allows all groups to have their own separate variance parameters.
#' "shared_all" allows all groups to share the same variance parameter.
#' "shared_QTLtype" allows all groups in one molecular QTL type to share the same variance parameter.
#'
#' @param max_snp_region Inf or integer. Maximum number of SNPs in a region. Default is
#' Inf, no limit. This can be useful if there are many SNPs in a region and you don't
#' have enough memory to run the program. This applies to the last rerun step
#' (using full SNPs and rerun susie for regions with strong gene signals) only.
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
#' @param verbose TRUE/FALSE. If TRUE, print detail messages
#'
#' @importFrom logging addHandler loginfo writeToFile
#'
#' @return a list of estimated parameters, fine-mapping results, and boundary genes
#'
#' @export
#'
ctwas_sumstats <- function(
    z_snp,
    weights,
    region_info,
    thin = 1,
    niter1 = 3,
    niter2 = 30,
    L = 5,
    init_group_prior = NULL,
    init_group_prior_var = NULL,
    group_prior_var_structure = c("independent","shared_all","shared_QTLtype"),
    max_snp_region = Inf,
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
  loginfo("Compute gene z-scores ...")
  z_gene <- compute_gene_z(z_snp, weights)

  # Get regionlist, which contains SNPs and genes assigned to each region
  #. downsample SNPs if thin < 1
  #. assign SNP and gene IDs, and z-scores to each region
  #. find boundary genes and adjust regionlist for boundary genes
  loginfo("Get regionlist ...")
  res <- get_regionlist(region_info,
                        z_snp,
                        z_gene,
                        weights,
                        maxSNP = max_snp_region,
                        trim_by = "random",
                        thin = thin,
                        minvar = 2,
                        adjust_boundary_genes = TRUE,
                        ncore = ncore)
  regionlist <- res$regionlist
  boundary_genes <- res$boundary_genes

  # Estimate parameters
  #. get regionlist for all the regions
  #. run EM for two rounds with thinned SNPs using L = 1
  loginfo("Estimating parameters...")
  param <- est_param(regionlist,
                     init_group_prior = init_group_prior,
                     init_group_prior_var = init_group_prior_var,
                     group_prior_var_structure = group_prior_var_structure,
                     niter1 = niter1,
                     niter2 = niter2,
                     p_single_effect = p_single_effect,
                     ncore = ncore)
  group_prior <- param$group_prior
  group_prior_var <- param$group_prior_var

  # Screen regions
  #. fine-map all regions with thinned SNPs
  #. select regions with strong non-SNP signals
  loginfo("Screening regions ...")
  screened_region_tags <- screen_regions(regionlist,
                                         region_info,
                                         weights,
                                         group_prior = group_prior,
                                         group_prior_var = group_prior_var,
                                         L = L,
                                         max_snp_region = max_snp_region,
                                         min_nonSNP_PIP = min_nonSNP_PIP,
                                         use_null_weight = use_null_weight,
                                         ncore = ncore)

  # Expand screened regionlist with all SNPs in the regions
  screened_regionlist <- regionlist[screened_region_tags]
  if (thin < 1){
    loginfo("Update regionlist with full SNPs for screened regions")
    screened_regionlist <- expand_regionlist(screened_regionlist,
                                             region_info,
                                             z_snp,
                                             z_gene,
                                             trim_by = "z",
                                             maxSNP = max_snp_region,
                                             ncore = ncore)
  }

  # Run fine-mapping for regions with strong gene signals using all SNPs
  #. save correlation matrices if save_cor is TRUE
  loginfo("Fine-mapping regions ...")
  finemap_res <- finemap_regions(screened_regionlist,
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

  ctwas_res <- list("param" = param,
                    "finemap_res" = finemap_res,
                    "boundary_genes" = boundary_genes)

  return(ctwas_res)

}

