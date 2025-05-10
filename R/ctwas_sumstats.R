#' @title cTWAS analysis using summary statistics
#'
#' @param z_snp A data frame with four columns: "id", "A1", "A2", "z".
#' giving the z scores for SNPs. "A1" is effect allele. "A2" is the other allele.
#'
#' @param weights a list of pre-processed prediction weights.
#'
#' @param region_info a data frame of region definitions.
#'
#' @param LD_map a data frame with filenames of LD matrices and SNP information for the regions.
#'
#' @param snp_map a list of data frames with SNP-to-region map for the reference.
#'
#' @param z_gene A data frame with columns: "id", "z", giving the z-scores for genes.
#'
#' @param niter_prefit the number of iterations of the E-M algorithm to perform during the initial parameter estimation step.
#'
#' @param niter the maximum number of iterations of the E-M algorithm to perform during the complete parameter estimation step.
#'
#' @param thin The proportion of SNPs to be used for estimating parameters and screening regions.
#'
#' @param L the number of effects for susie during the fine mapping steps.
#'
#' @param init_group_prior a vector of initial values of prior inclusion probabilities for different groups.
#'
#' @param init_group_prior_var a vector of initial values of prior variances for different groups.
#'
#' @param group_prior_var_structure a string indicating the structure to put on the prior variance parameters.
#' "shared_all" allows all groups to share the same variance parameter.
#' "shared_type" allows all groups in one molecular QTL type to share the same variance parameter.
#' "shared_context" allows all groups in one context (tissue, cell type, condition) to share the same variance parameter.
#' "shared_nonSNP" allows all non-SNP groups to share the same variance parameter.
#' "independent" allows all groups to have their own separate variance parameters.
#' "fixed" sets prior variance parameters to values in \code{init_group_prior_var}.
#'
#' @param min_nonSNP_PIP Regions with non-SNP PIP >= \code{min_nonSNP_PIP}
#' will be selected to run finemapping using all SNPs.
#'
#' @param min_pval Keep regions with minimum p-values from z_snp and z_gene < \code{min_pval}.
#'
#' @param min_p_single_effect Regions with probability greater than \code{min_p_single_effect} of
#' having 1 or fewer effects will be used for parameter estimation.
#'
#' @param maxSNP Inf or integer. Maximum number of SNPs in a region. Default is
#' Inf, no limit. This can be useful if there are many SNPs in a region and you don't
#' have enough memory to run the program.
#'
#' @param min_var minimum number of variables (SNPs and genes) in a region
#' when estimating paramters and screening regions.
#'
#' @param min_gene minimum number of genes in a region
#' when estimating paramters and screening regions.
#'
#' @param min_group_size Minimum number of genes in a group.
#' Groups with number of genes < \code{min_group_size} will be removed for the analysis.
#'
#' @param null_method Method to compute null model, options: "ctwas", "susie" or "none".
#'
#' @param EM_tol A small, non-negative number specifying the convergence
#'   tolerance of log-likelihood for the EM iterations.
#'
#' @param coverage A number between 0 and 1 specifying the \dQuote{coverage} of
#' the estimated confidence sets.
#'
#' @param min_abs_corr Minimum absolute correlation allowed in a credible set.
#'
#' @param include_prior If TRUE, include priors in finemapping results.
#'
#' @param include_susie_result If TRUE, include the "susie" result object in finemapping results.
#'
#' @param LD_format file format for LD matrix. If "custom", use a user defined
#' \code{LD_loader_fun()} function to load LD matrix.
#'
#' @param LD_loader_fun a user defined function to load LD matrix when \code{LD_format = "custom"}.
#'
#' @param snpinfo_loader_fun a user defined function to load SNP information file,
#' if SNP information files are not in standard cTWAS reference format.
#'
#' @param force_compute_cor If TRUE, force computing correlation (R) matrices.
#'
#' @param save_cor If TRUE, save correlation (R) matrices to \code{cor_dir}.
#'
#' @param cor_dir The directory to store correlation (R) matrices.
#'
#' @param outputdir The directory to store output. If specified, save outputs to the directory.
#'
#' @param outname The output name.
#'
#' @param ncore The number of cores used to parallelize computing over regions.
#'
#' @param ncore_LD The number of cores used to parallelize computing correlation matrices,
#' in screening regions and fine-mapping steps with LD.
#'
#' @param seed seed for random sampling when thinning the SNPs in region data.
#'
#' @param logfile The log filename. If NULL, print log info on screen.
#'
#' @param verbose If TRUE, print detailed messages.
#'
#' @param ... Additional arguments of \code{susie_rss}.
#'
#' @importFrom logging addHandler loginfo logwarn writeToFile
#' @importFrom utils packageVersion
#'
#' @return a list, including z_gene, estimated parameters, region_data,
#' cross-boundary genes, screening region results, and fine-mapping results.
#'
#' @export
#'
ctwas_sumstats <- function(
    z_snp,
    weights,
    region_info,
    LD_map,
    snp_map,
    z_gene = NULL,
    thin = 1,
    niter_prefit = 3,
    niter = 50,
    L = 5,
    init_group_prior = NULL,
    init_group_prior_var = NULL,
    group_prior_var_structure = c("shared_all", "shared_type", "shared_context", "shared_nonSNP", "independent"),
    min_nonSNP_PIP = 0.5,
    min_pval = 5e-8,
    min_p_single_effect = 0.8,
    maxSNP = Inf,
    min_var = 2,
    min_gene = 1,
    min_group_size = 100,
    null_method = c("ctwas", "susie", "none"),
    EM_tol = 1e-4,
    coverage = 0.95,
    min_abs_corr = 0.1,
    include_prior = FALSE,
    include_susie_result = FALSE,
    LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
    LD_loader_fun = NULL,
    snpinfo_loader_fun = NULL,
    force_compute_cor = FALSE,
    save_cor = FALSE,
    cor_dir = NULL,
    outputdir = NULL,
    outname = "ctwas",
    ncore = 1,
    ncore_LD = max(ncore-1,1),
    seed = 99,
    logfile = NULL,
    verbose = FALSE,
    ...){

  if (!is.null(logfile)) {
    addHandler(writeToFile, file=logfile, level='DEBUG')
  }

  loginfo("Running cTWAS ...")
  loginfo("ctwas version: %s", packageVersion("ctwas"))

  # check inputs
  group_prior_var_structure <- match.arg(group_prior_var_structure)
  LD_format <- match.arg(LD_format)
  null_method <- match.arg(null_method)

  if (anyNA(z_snp))
    stop("z_snp contains missing values!")

  if (!inherits(snp_map,"list"))
    stop("'snp_map' should be a list.")

  if (!inherits(LD_map,"data.frame"))
    stop("'LD_map' should be a data frame")

  if (!inherits(weights,"list"))
    stop("'weights' should be a list.")

  if (any(sapply(weights, is.null)))
    stop("weights contain NULL, remove empty weights!")

  if (thin > 1 | thin <= 0)
    stop("thin needs to be in (0,1]")

  stopifnot(all(file.exists(LD_map$LD_file)))
  stopifnot(all(file.exists(LD_map$SNP_file)))

  if (!is.null(outputdir)) {
    dir.create(outputdir, showWarnings=FALSE, recursive=TRUE)
  }

  loginfo("ncore = %d", ncore)
  loginfo("ncore_LD = %d", ncore_LD)

  # check the number of SNPs in z_snp and weights
  count_snp_res <- count_z_snp_weights(z_snp, weights)
  loginfo("%d SNPs in z_snp", count_snp_res$n_snps_z)
  loginfo("%d SNPs in weights", count_snp_res$n_snps_weights)

  if (count_snp_res$n_snps_weights == 0)
    stop("Please check z_snp! No SNPs in z_snp outside weights! ")

  if (count_snp_res$N_snps_not_in_weights < count_snp_res$n_snps_weights)
    logwarn("SNPs outside weights are less than SNPs in weights!")

  # Compute gene z-scores
  if (is.null(z_gene)) {
    z_gene <- compute_gene_z(z_snp, weights, ncore = ncore)
    # filter groups with too few genes
    z_gene <- filter_z_gene_by_group_size(z_gene, min_group_size)
    if (!is.null(outputdir)) {
      saveRDS(z_gene, file.path(outputdir, paste0(outname, ".z_gene.RDS")))
    }
  } else {
    # filter groups with too few genes
    z_gene <- filter_z_gene_by_group_size(z_gene, min_group_size)
  }

  if (anyNA(z_gene))
    stop("z_gene contains missing values!")

  # Get region_data, which contains SNPs and genes assigned to each region
  #. downsample SNPs if thin < 1
  #. assign SNP and gene IDs, and z-scores to each region
  #. find boundary genes and adjust region_data for boundary genes
  region_data_res <- assemble_region_data(region_info,
                                          z_snp,
                                          z_gene,
                                          weights,
                                          snp_map,
                                          thin = thin,
                                          maxSNP = maxSNP,
                                          min_group_size = min_group_size,
                                          trim_by = "random",
                                          adjust_boundary_genes = TRUE,
                                          ncore = ncore,
                                          seed = seed)
  region_data <- region_data_res$region_data
  boundary_genes <- region_data_res$boundary_genes
  rm(region_data_res)
  if (!is.null(outputdir)) {
    saveRDS(region_data, file.path(outputdir, paste0(outname, ".region_data.thin", thin, ".RDS")))
    saveRDS(boundary_genes, file.path(outputdir, paste0(outname, ".boundary_genes.RDS")))
  }

  # Estimate parameters
  #. get region_data for all the regions
  #. run EM for two rounds using the SER model (L = 1)
  param <- est_param(region_data,
                     init_group_prior = init_group_prior,
                     init_group_prior_var = init_group_prior_var,
                     group_prior_var_structure = group_prior_var_structure,
                     niter_prefit = niter_prefit,
                     niter = niter,
                     min_var = min_var,
                     min_gene = min_gene,
                     min_group_size = min_group_size,
                     min_p_single_effect = min_p_single_effect,
                     null_method = null_method,
                     EM_tol = EM_tol,
                     ncore = ncore,
                     verbose = verbose)
  group_prior <- param$group_prior
  group_prior_var <- param$group_prior_var
  if (!is.null(outputdir)) {
    saveRDS(param, file.path(outputdir, paste0(outname, ".param.RDS")))
  }

  # expand regions with all SNPs before screening and fine-mapping
  if (thin < 1) {
    all_region_data <- expand_region_data(region_data,
                                          snp_map,
                                          z_snp,
                                          maxSNP = maxSNP,
                                          ncore = ncore)
  } else {
    all_region_data <- region_data
  }

  # Screen regions
  #. fine-map all regions without LD (using SER model, L = 1)
  #. select regions with strong non-SNP signals
  screen_res <- screen_regions(all_region_data,
                               group_prior = group_prior,
                               group_prior_var = group_prior_var,
                               min_var = min_var,
                               min_gene = min_gene,
                               min_nonSNP_PIP = min_nonSNP_PIP,
                               min_pval = min_pval,
                               null_method = null_method,
                               ncore = ncore,
                               verbose = verbose)
  screened_region_data <- screen_res$screened_region_data
  if (!is.null(outputdir)) {
    saveRDS(screen_res, file.path(outputdir, paste0(outname, ".screen_res.RDS")))
  }

  # Run fine-mapping for regions with strong gene signals using all SNPs
  #. save correlation matrices if save_cor is TRUE
  if (length(screened_region_data) > 0){
    res <- finemap_regions(screened_region_data,
                           LD_map = LD_map,
                           weights = weights,
                           group_prior = group_prior,
                           group_prior_var = group_prior_var,
                           L = L,
                           null_method = null_method,
                           coverage = coverage,
                           min_abs_corr = min_abs_corr,
                           include_prior = include_prior,
                           include_susie_result = include_susie_result,
                           force_compute_cor = force_compute_cor,
                           save_cor = save_cor,
                           cor_dir = cor_dir,
                           LD_format = LD_format,
                           LD_loader_fun = LD_loader_fun,
                           snpinfo_loader_fun = snpinfo_loader_fun,
                           ncore = ncore_LD,
                           verbose = verbose,
                           ...)
    finemap_res <- res$finemap_res
    susie_alpha_res <- res$susie_alpha_res
    susie_res <- res$susie_res
    if (!is.null(outputdir)) {
      saveRDS(finemap_res, file.path(outputdir, paste0(outname, ".finemap_res.RDS")))
      if (!is.null(susie_alpha_res))
        saveRDS(susie_alpha_res, file.path(outputdir, paste0(outname, ".susie_alpha_res.RDS")))
      if (!is.null(susie_res))
        saveRDS(susie_res, file.path(outputdir, paste0(outname, ".susie_res.RDS")))
    }
  } else {
    loginfo("No regions selected for fine-mapping.")
    finemap_res <- NULL
    susie_alpha_res <- NULL
    susie_res <- NULL
  }

  return(list("z_gene" = z_gene,
              "param" = param,
              "finemap_res" = finemap_res,
              "susie_alpha_res" = susie_alpha_res,
              "susie_res" = susie_res,
              "region_data" = region_data,
              "boundary_genes" = boundary_genes,
              "screen_res" = screen_res))

}

