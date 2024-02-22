#' Estimate cTWAS parameters
#'
#' @param z_snp a data frame with four columns: "id", "A1", "A2", "z".
#' giving the z scores for SNPs. "A1" is effect allele. "A2" is the other allele.
#'
#' @param region_info a data frame of region definition and associated file names
#'
#' @param regionlist a list object indexing regions, variants and genes.
#'
#' @param z_gene a data frame with two columns: "id", "z". giving the z scores for genes.
#' Optionally, a "type" column can also be supplied; this is for using multiple sets of weights
#'
#' @param gene_info a data frame of gene information obtained from \code{compute_gene_z}
#'
#' @param weight a string, pointing to a directory with the FUSION/TWAS format of weights, or a .db file in predictdb format.
#' A vector of multiple sets of weights in PredictDB format can also be specified; genes will have their filename appended
#' to their gene name to ensure IDs are unique.

#' @param weight_format a string, the format of weight, PredictDB or FUSION
#'
#' @param method a string, blup/bslmm/lasso/top1/enet/best. This option is only used for FUSION weights.
#' "best" means the method giving the best cross #' validation R^2. Note that top1 uses only the weight
#' with largest effect.
#'
#' @param scale_by_ld_variance TRUE/FALSE. If TRUE, PredictDB weights are scaled by genotype variance, which is the default
#' behavior for PredictDB
#'
#' @param thin The proportion of SNPs to be used for the parameter estimation and
#' initial screening region steps.
#' Smaller \code{thin} parameters reduce runtime at the expense of accuracy.
#' The fine mapping step is rerun using full SNPs for regions with strong gene signals.
#'
#' @param prob_single Blocks with probability greater than \code{prob_single} of
#' having 1 or fewer effects will be used for parameter estimation
#'
#' @param niter1 the number of iterations of the E-M algorithm to perform during the initial parameter estimation step
#'
#' @param niter2 the number of iterations of the E-M algorithm to perform during the complete parameter estimation step
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
#' @param use_null_weight TRUE/FALSE. If TRUE, allow for a probability of no effect in susie
#'
#' @param coverage A number between 0 and 1 specifying the \dQuote{coverage} of the estimated confidence sets
#'
#' @param min_abs_corr Minimum absolute correlation allowed in a
#'   credible set. The default, 0.5, corresponds to a squared
#'   correlation of 0.25, which is a commonly used threshold for
#'   genotype data in genetic studies.
#'
#' @param ncore The number of cores used to parallelize computation over regions
#'
#' @param outputdir a string, the directory to store the correlation files
#'
#' @param outname a string, the output name for the correlation files
#'
#' @param logfile the log file, if NULL will print log info on screen
#'
#' @importFrom logging addHandler loginfo writeToFile
#'
#' @return a list with estimated parameters, and updated region_info with correlation file names for each region.
#'
#' @export
#'
est_param <- function(
    z_snp,
    region_info,
    regionlist = NULL,
    z_gene = NULL,
    gene_info = NULL,
    weight = NULL,
    weight_format = c("PredictDB", "FUSION"),
    method = c("lasso", "blup", "bslmm", "top1", "enet", "best"),
    scale_by_ld_variance = TRUE,
    thin = 1,
    prob_single = 0.8,
    niter1 = 3,
    niter2 = 30,
    L = 5,
    group_prior = NULL, # DO WE NEED to specify the priors here?
    group_prior_var = NULL,
    group_prior_var_structure = c("independent","shared_all","shared+snps","shared_QTLtype"),
    use_null_weight = TRUE,
    coverage = 0.95,
    min_abs_corr = 0.5,
    ncore = 1,
    outputdir = getwd(),
    outname = NULL,
    logfile = NULL){

  if (!is.null(logfile)){
    addHandler(writeToFile, file= logfile, level='DEBUG')
  }

  loginfo('Estimating cTWAS parameters ... ')

  group_prior_var_structure <- match.arg(group_prior_var_structure)

  # compute gene z-scores if not available
  if (is.null(z_gene)) {
    loginfo("Computing gene z-scores ...")
    compute_gene_z_res <- compute_gene_z(z_snp = z_snp,
                                         weight = weight,
                                         region_info = region_info,
                                         weight_format = weight_format,
                                         method = method,
                                         scale_by_ld_variance=scale_by_ld_variance,
                                         ncore=ncore)
    z_gene <- compute_gene_z_res$z_gene
    z_snp <- compute_gene_z_res$z_snp
    gene_info <- compute_gene_z_res$gene_info
    rm(compute_gene_z_res)
  }

  # combine z-scores of SNPs and genes
  zdf <- combine_z(z_gene, z_snp)

  if (thin <= 0 | thin > 1){
    stop("thin value needs to be in (0,1]")
  }

  # if correlation files not available, computing correlation matrices and generating regionlist
  if (is.null(region_info$cor_file)) {
    loginfo("Compute correlation matrices and generate regionlist with thin = %.2f", thin)
    compute_cor_res <- compute_cor(region_info = region_info,
                                   gene_info = gene_info,
                                   weight_list = weight,
                                   select = zdf$id,
                                   thin = thin,
                                   minvar = 2,
                                   outname = outname,
                                   outputdir = outputdir,
                                   merge = FALSE,
                                   ncore = ncore)

    regionlist <- compute_cor_res$regionlist
    region_info <- compute_cor_res$region_info # updated region_info containing correlation file names for each region
    # saveRDS(regionlist, file=file.path(outputdir, paste0(outname, ".regionlist.RDS")))
    rm(compute_cor_res)
  }

  # Run susieI with L = 1 for a few (niter1) iterations, getting rough estimates
  loginfo("Run susieI for %d iterations, getting rough estimates ...", niter1)

  if (!is.null(group_prior)){
    group_prior["SNP"] <- group_prior["SNP"]/thin
  }

  susieI_res <- ctwas_susieI_rss(zdf = zdf,
                                 regionlist = regionlist,
                                 region_info = region_info,
                                 niter = niter1,
                                 L = 1,
                                 group_prior = group_prior, # DO WE NEED to take the priors here?
                                 group_prior_var = group_prior_var,
                                 group_prior_var_structure = group_prior_var_structure,
                                 estimate_group_prior = TRUE,
                                 estimate_group_prior_var = TRUE,
                                 use_null_weight = use_null_weight,
                                 coverage = coverage,
                                 min_abs_corr = min_abs_corr,
                                 ncore = ncore)

  param <- susieI_res$param
  group_prior <- param$group_prior
  group_prior_var <- param$group_prior_var
  print("Roughly estimated group_prior: \n")
  print(group_prior)
  print("Roughly estimated group_prior_var: \n")
  print(group_prior_var)
  rm(susieI_res, param)

  # filter regions based on prob_single
  filtered_regionlist <- filter_regions(regionlist = regionlist,
                                        group_prior = group_prior,
                                        prob_single = prob_single,
                                        zdf = zdf)
  n_filtered_regions <- sum(unlist(lapply(filtered_regionlist, length)))
  loginfo("%d regions left after filtering with prob_single = %.2f",
          n_filtered_regions, prob_single)

  # Run susieI with L = 1 for more (niter2) iterations, getting rough estimates
  loginfo("Run susieI with L = 1 for %d iterations, getting accurate estimates ...", niter2)
  susieI_res <- ctwas_susieI_rss(zdf = zdf,
                                 regionlist = filtered_regionlist,
                                 region_info = region_info,
                                 niter = niter2,
                                 L = 1,
                                 group_prior = group_prior, # DO WE NEED to take the priors here?
                                 group_prior_var = group_prior_var,
                                 group_prior_var_structure = group_prior_var_structure,
                                 estimate_group_prior = TRUE,
                                 estimate_group_prior_var = TRUE,
                                 use_null_weight = use_null_weight,
                                 coverage = coverage,
                                 min_abs_corr = min_abs_corr,
                                 ncore = ncore)

  param <- susieI_res$param
  print("Estimated group_prior: \n")
  print(param$group_prior)
  print("Estimated group_prior_var: \n")
  print(param$group_prior_var)

  return(list("param" = param,
              "region_info" = region_info,
              "regionlist" = regionlist,
              "filtered_regionlist" = filtered_regionlist))

}

