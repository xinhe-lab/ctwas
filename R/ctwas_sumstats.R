#' Causal inference for TWAS using summary statistics
#'
#' @param z_snp A data frame with four columns: "id", "A1", "A2", "z".
#' giving the z scores for snps. "A1" is effect allele. "A2" is the other allele.
#'
#' @param region_info a data frame of region definition and associated LD file names
#'
#' @param weight_list a list of weights
#'
#' @param weight_info a data frame of weight information

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
#' @param min_nonSNP_PIP Regions with non-SNP PIP >= \code{min_nonSNP_PIP}
#' will be selected to run finemapping using full SNPs.
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
#' @param save_cor TRUE/FALSE. If TRUE, save correlation (R) matrices to \code{outputdir}
#'
#' @param cor_dir a string, the directory to store correlation (R) matrices
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
    region_info,
    weight_list,
    weight_info,
    niter1 = 3,
    niter2 = 30,
    thin = 1,
    L = 5,
    group_prior = NULL,
    group_prior_var = NULL,
    group_prior_var_structure = c("independent","shared_all","shared_QTLtype"),
    max_snp_region = Inf,
    min_nonSNP_PIP = 0.5,
    prob_single = 0.8,
    use_null_weight = T,
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
  res <- compute_gene_z(z_snp = z_snp,
                        region_info = region_info,
                        weight_list = weight_list,
                        weight_info = weight_info,
                        ncore = ncore)
  z_gene <- res$z_gene
  gene_info <- res$gene_info
  rm(res)

  # Estimate parameters
  #. get regionlist for all the regions
  #. run EM for two rounds with thinned SNPs using L = 1
  param <- est_param(z_snp = z_snp,
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
  group_prior <- param$group_prior
  group_prior_var <- param$group_prior_var
  regionlist <- param$regionlist
  weight_list <- param$weight_list
  boundary_genes <- param$boundary_genes
  rm(res)

  # Screen regions
  #. fine-map all regions with thinned SNPs
  #. select regions with strong non-SNP signals
  res <- screen_regions(z_snp = z_snp,
                        z_gene = z_gene,
                        gene_info = gene_info,
                        region_info = region_info,
                        regionlist = regionlist,
                        weight_list = weight_list,
                        thin = thin,
                        max_snp_region = max_snp_region,
                        min_nonSNP_PIP = min_nonSNP_PIP,
                        L = L,
                        group_prior = group_prior,
                        group_prior_var = group_prior_var,
                        use_null_weight = use_null_weight,
                        coverage = coverage,
                        min_abs_corr = min_abs_corr,
                        ncore = ncore)
  screened_regionlist <- res$screened_regionlist
  screened_region_tags <- res$screened_region_tags
  weak_region_finemap_res <- res$weak_region_finemap_res
  rm(res)

  # Run fine-mapping for regions with strong gene signals using all SNPs
  #. save correlation matrices if save_cor is TRUE

  # adjust group_prior parameter to account for thin argument
  group_prior["SNP"] <- group_prior["SNP"] * thin

  strong_region_finemap_res <- finemap_regions(z_snp = z_snp,
                                               z_gene = z_gene,
                                               gene_info = gene_info,
                                               regionlist = screened_regionlist,
                                               weight_list = weight_list,
                                               L = L,
                                               group_prior = group_prior,
                                               group_prior_var = group_prior_var,
                                               use_null_weight = use_null_weight,
                                               coverage = coverage,
                                               min_abs_corr = min_abs_corr,
                                               save_cor = save_cor,
                                               cor_dir = cor_dir,
                                               ncore = ncore,
                                               verbose = verbose)

  ctwas_res <- list("strong_region_finemap_res" = strong_region_finemap_res,
                    "weak_region_finemap_res" = weak_region_finemap_res,
                    "param" = param,
                    "z_gene" = z_gene,
                    "gene_info" = gene_info,
                    "region_info" = region_info,
                    "regionlist" = regionlist,
                    "weight_list" = weight_list,
                    "boundary_genes" = boundary_genes,
                    "screened_regionlist" = screened_regionlist,
                    "screened_region_tags" = screened_region_tags)
  return(ctwas_res)

}

