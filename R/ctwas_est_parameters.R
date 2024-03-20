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
#' @param weight_list a list of weights
#'
#' @param thin The proportion of SNPs to be used for the parameter estimation and
#' initial screening region steps.
#' Smaller \code{thin} parameters reduce runtime at the expense of accuracy.
#' The fine mapping step is rerun using full SNPs for regions with strong gene signals.
#'
#' @param max_snp_region Inf or integer. Maximum number of SNPs in a region. Default is
#' Inf, no limit. This can be useful if there are many SNPs in a region and you don't
#' have enough memory to run the program. This applies to the last rerun step
#' (using full SNPs and rerun susie for regions with strong gene signals) only.
#'
#' @param prob_single Blocks with probability greater than \code{prob_single} of
#' having 1 or fewer effects will be used for parameter estimation
#'
#' @param niter1 the number of iterations of the E-M algorithm to perform during the initial parameter estimation step
#'
#' @param niter2 the number of iterations of the E-M algorithm to perform during the complete parameter estimation step
#'
#' @param group_prior a vector of two prior inclusion probabilities for SNPs and genes.
#'
#' @param group_prior_var a vector of two prior variances for SNPs and gene effects.
#'
#' @param group_prior_var_structure a string indicating the structure to put on the prior variance parameters.
#' "independent" is the default and allows all groups to have their own separate variance parameters.
#' "shared_all" allows all groups to share the same variance parameter.
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
#' @param max_iter Maximum number of IBSS iterations to perform.
#'
#' @param ncore The number of cores used to parallelize computation over regions
#'
#' @param logfile the log file, if NULL will print log info on screen
#'
#' @importFrom logging addHandler loginfo writeToFile
#'
#' @return a list with estimated parameters, regionlist, updated weight list and boundary genes.
#'
#' @export
#'
est_param <- function(
    z_snp,
    region_info,
    regionlist = NULL,
    z_gene = NULL,
    gene_info = NULL,
    weight_list = NULL,
    weight_info = NULL,
    thin = 1,
    max_snp_region = Inf,
    prob_single = 0.8,
    niter1 = 3,
    niter2 = 30,
    group_prior = NULL,
    group_prior_var = NULL,
    group_prior_var_structure = c("independent","shared_all","shared_QTLtype"),
    use_null_weight = TRUE,
    coverage = 0.95,
    min_abs_corr = 0.5,
    max_iter = 1,
    ncore = 1,
    logfile = NULL){

  if (!is.null(logfile)){
    addHandler(writeToFile, file= logfile, level='DEBUG')
  }

  loginfo('Estimating parameters ... ')

  group_prior_var_structure <- match.arg(group_prior_var_structure)

  if (thin <= 0 | thin > 1){
    stop("thin value needs to be in (0,1]")
  }

  # compute gene z-scores if not available
  if (is.null(z_gene)) {
    loginfo("Computing gene z-scores ...")
    res <- compute_gene_z(z_snp, region_info, weight_list, weight_info, ncore = ncore)
    z_gene <- res$z_gene
    gene_info <- res$gene_info
    rm(res)
  }

  # combine z-scores of SNPs and genes
  zdf <- combine_z(z_snp, z_gene)

  # get regionlist if not available
  if (is.null(regionlist)) {
    loginfo("Get regionlist with thin = %s", thin)
    res <- get_regionlist(region_info,
                          gene_info,
                          weight_list = weight_list,
                          select = zdf$id,
                          thin = thin,
                          maxSNP = max_snp_region,
                          minvar = 2,
                          adjust_boundary = TRUE)
    regionlist <- res$regionlist
    weight_list <- res$weight_list
    boundary_genes <- res$boundary_genes
    rm(res)
  }

  # Run EM for a few (niter1) iterations, getting rough estimates
  loginfo("Run EM for %d iterations, getting rough estimates ...", niter1)
  if (!is.null(group_prior)){
    group_prior["SNP"] <- group_prior["SNP"]/thin
  }

  EM1_res <- ctwas_EM(zdf,
                      regionlist,
                      region_info,
                      gene_info,
                      niter = niter1,
                      group_prior = group_prior,
                      group_prior_var = group_prior_var,
                      group_prior_var_structure = group_prior_var_structure,
                      use_null_weight = use_null_weight,
                      coverage = coverage,
                      min_abs_corr = min_abs_corr,
                      max_iter = max_iter,
                      ncore = ncore)
  group_prior <- EM1_res$group_prior
  group_prior_var <- EM1_res$group_prior_var
  loginfo("Roughly estimated group_prior {%s}: {%s}", names(group_prior), group_prior)
  loginfo("Roughly estimated group_prior_var {%s}: {%s}", names(group_prior_var), group_prior_var)
  group_size <- table(EM1_res$EM_susie_res$type)
  loginfo("group_size {%s}: {%s}", names(group_size), group_size)

  # filter regions based on prob_single
  filtered_regionlist <- filter_regions(regionlist, zdf, group_prior, prob_single = prob_single)

  # Run EM for more (niter2) iterations, getting rough estimates
  loginfo("Run EM for %d iterations on %d filtered regions, getting accurate estimates ...",
          niter2, length(filtered_regionlist))

  EM2_res <- ctwas_EM(zdf,
                      filtered_regionlist,
                      region_info,
                      gene_info,
                      niter = niter2,
                      group_prior = group_prior,
                      group_prior_var = group_prior_var,
                      group_prior_var_structure = group_prior_var_structure,
                      use_null_weight = use_null_weight,
                      coverage = coverage,
                      min_abs_corr = min_abs_corr,
                      max_iter = max_iter,
                      ncore = ncore)

  group_prior <- EM2_res$group_prior
  group_prior_var <- EM2_res$group_prior_var
  group_prior_rec <- EM2_res$group_prior_rec
  group_prior_var_rec <- EM2_res$group_prior_var_rec
  group_prior_var_structure <- EM2_res$group_prior_var_structure
  loginfo("Estimated group_prior {%s}: {%s}", names(group_prior), group_prior)
  loginfo("Estimated group_prior_var {%s}: {%s}", names(group_prior_var), group_prior_var)

  param <- list("group_prior" = group_prior,
                "group_prior_var" = group_prior_var,
                "group_prior_rec" = group_prior_rec,
                "group_prior_var_rec" = group_prior_var_rec,
                "group_prior_var_structure" = group_prior_var_structure,
                "group_size" = group_size,
                "thin" = thin,
                "regionlist" = regionlist,
                "filtered_regionlist" = filtered_regionlist)

  if (!is.null(weight_list)){
    param[["weight_list"]] <- weight_list
  }
  if (!is.null(boundary_genes)){
    param[["boundary_genes"]] <- boundary_genes
  }

  return(param)
}

