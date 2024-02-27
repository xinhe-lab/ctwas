#' Causal inference for TWAS using summary statistics
#'
#' @param z_snp A data frame with four columns: "id", "A1", "A2", "z".
#' giving the z scores for snps. "A1" is effect allele. "A2" is the other allele.
#'
#' @param region_info a data frame of region definition and associated LD file names
#'
#' @param weights a string, pointing to a directory with the weights (.db file) in predictdb format.
#' A vector of multiple sets of weights in PredictDB format can also be specified;
#' genes will have their filename appended to their gene name to ensure IDs are unique.
#'
#' @param niter1 the number of iterations of the E-M algorithm to perform during the initial parameter estimation step
#'
#' @param niter2 the number of iterations of the E-M algorithm to perform during the complete parameter estimation step
#'
#' @param thin The proportion of SNPs to be used for parameter estimation and initial screening regions.
#' Smaller \code{thin} parameters reduce runtime at the expense of accuracy.
#' The fine mapping step is rerun using full SNPs
#' for regions with strong gene signals; see \code{rerun_gene_PIP}.
#'
#' @param L the number of effects for susie during the fine mapping steps
#'
#' @param group_prior a vector of two prior inclusion probabilities for SNPs and genes.
#'
#' @param group_prior_var a vector of two prior variances for SNPs and gene effects.
#'
#' @param group_prior_var_structure a string indicating the structure to put on the prior variance parameters.
#' "independent" is the default and allows all groups to have their own separate variance parameters.
#' "shared" allows all groups to share the same variance parameter.
#' "shared+snps" allows all groups to share the same variance parameter, and this variance parameter is also shared with SNPs.
#' "shared_type" allows all groups in one molecular QTL type to share the same variance parameter.
#'
#' @param max_snp_region Inf or integer. Maximum number of SNPs in a region. Default is
#' Inf, no limit. This can be useful if there are many SNPs in a region and you don't
#' have enough memory to run the program. This applies to the last rerun step
#' (using full SNPs and rerun susie for regions with strong gene signals) only.
#'
#' @param rerun_gene_PIP if thin <1, will rerun regions with the max gene PIP > \code{rerun_gene_PIP}
#' using full SNPs. if \code{rerun_gene_PIP} is 0, then
#' all blocks will rerun with full SNPs
#'
#' @param prob_single Regions with probability greater than \code{prob_single} of
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
#' @param save_R TRUE/FALSE. If TRUE, save correlation (R) matrices to \code{outputdir}
#'
#' @param outputdir a string, the directory to store correlation (R) matrices
#'
#' @param outname a string, the output name
#'
#' @param logfile the log file, if NULL will print log info on screen
#'
#' @importFrom logging addHandler loginfo writeToFile
#'
#' @return a list of estimated parameters, fine-mapping results,
#'. updated region info, and imputed gene scores
#'
#' @export
#'
ctwas_sumstats <- function(
    z_snp,
    weights,
    region_info,
    niter1 = 3,
    niter2 = 30,
    thin = 1,
    L = 5,
    group_prior = NULL,
    group_prior_var = NULL,
    group_prior_var_structure = c("independent","shared_all","shared+snps","shared_QTLtype"),
    max_snp_region = Inf,
    rerun_gene_PIP = 0.5,
    prob_single = 0.8,
    use_null_weight = T,
    coverage = 0.95,
    min_abs_corr = 0.5,
    ncore = 1,
    save_R = FALSE,
    outputdir = getwd(),
    outname = NULL,
    logfile = NULL){

  if (!is.null(logfile)){
    addHandler(writeToFile, file=logfile, level='DEBUG')
  }

  # Compute gene z-scores
  res <- compute_gene_z(z_snp = z_snp,
                        weights = weights,
                        region_info = region_info,
                        ncore = ncore)
  z_gene <- res$z_gene
  z_snp <- res$z_snp
  gene_info <- res$gene_info
  rm(res)

  # Estimate parameters
  #.  including computing correlation matrices for thinned SNPs
  res <- est_param(z_snp = z_snp,
                   region_info = region_info,
                   z_gene = z_gene,
                   gene_info = gene_info,
                   thin = thin,
                   group_prior = group_prior,
                   group_prior_var = group_prior_var,
                   group_prior_var_structure = group_prior_var_structure,
                   niter1 = niter1,
                   niter2 = niter2,
                   ncore = ncore)

  param <- res$param
  group_prior <- param$group_prior
  group_prior_var <- param$group_prior_var
  region_info <- res$region_info
  regionlist_thinned <- res$regionlist
  rm(res)

  # Screen regions
  #. index SNPs and genes for all the regions
  #. fine-map all regions with thinned SNPs with L = 1,
  #. select regions with strong gene signals
  res <- screen_regions(z_snp = z_snp,
                        z_gene = z_gene,
                        region_info = region_info,
                        gene_info = gene_info,
                        weights = weights,
                        regionlist = regionlist_thinned,
                        thin = thin,
                        max_snp_region = max_snp_region,
                        rerun_gene_PIP = rerun_gene_PIP,
                        L = 1,
                        group_prior = group_prior,
                        group_prior_var = group_prior_var,
                        use_null_weight = use_null_weight,
                        coverage = coverage,
                        min_abs_corr = min_abs_corr,
                        ncore = ncore)

  region_info <- res$region_info
  regionlist_allSNPs <- res$regionlist_allSNPs
  screened_region_tags <- res$screened_region_tags
  finemap_weak_res <- res$finemap_weak_res
  rm(res)

  # convert SNP prior
  group_prior["SNP"] <- group_prior["SNP"] * thin
  param$group_prior <- group_prior

  # Run fine-mapping for regions with strong gene signals using all SNPs
  #.  save correlation (R) matrices if save_R is TRUE
  for (region_tag in screened_region_tags) {
    finemap_strong_res <- finemap_region(z_snp = z_snp,
                                         z_gene = z_gene,
                                         region_info = region_info,
                                         gene_info = gene_info,
                                         regionlist = regionlist_allSNPs,
                                         region_tag = region_tag,
                                         L = L,
                                         group_prior = group_prior,
                                         group_prior_var = group_prior_var,
                                         use_null_weight = use_null_weight,
                                         coverage = coverage,
                                         min_abs_corr = min_abs_corr,
                                         save_R = save_R,
                                         outputdir = outputdir,
                                         outname = outname)
  }

  return(list("param" = param,
              "finemap_weak_res" = finemap_weak_res,
              "finemap_strong_res" = finemap_strong_res,
              "z_gene" = z_gene,
              "gene_info" = gene_info,
              "region_info" = region_info,
              "regionlist" = regionlist,
              "regionlist_allSNPs" = regionlist_allSNPs))

}

