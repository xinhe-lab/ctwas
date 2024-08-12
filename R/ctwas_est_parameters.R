#' @title Estimates cTWAS parameters using EM
#'
#' @param region_data a list object indexing regions, variants and genes.
#'
#' @param init_group_prior a vector of initial values of prior inclusion probabilities for SNPs and genes.
#'
#' @param init_group_prior_var a vector of initial values of prior variances for SNPs and gene effects.
#'
#' @param niter_prefit the number of iterations of the E-M algorithm to perform
#' during the initial parameter estimation step
#'
#' @param niter the number of iterations of the E-M algorithm to perform during
#' the complete parameter estimation step
#'
#' @param min_p_single_effect Regions with probability >= \code{min_p_single_effect}
#' of having at most one causal effect will be selected for the final EM step.
#'
#' @param use_null_weight If TRUE, allow for a probability of no effect in susie
#'
#' @param min_snps minimum number of SNPs in a region.
#'
#' @param min_genes minimum number of genes in a region.
#'
#' @param ncore The number of cores used to parallelize computation over regions
#'
#' @param logfile The log filename. If NULL, print log info on screen.
#'
#' @param verbose If TRUE, print detail messages
#'
#' @param ... Additional arguments of \code{susie_rss}.
#'
#' @importFrom logging addHandler loginfo writeToFile
#'
#' @return a list with estimated parameters
#'
#' @export
#'
est_param <- function(
    region_data,
    init_group_prior = NULL,
    init_group_prior_var = NULL,
    niter_prefit = 3,
    niter = 30,
    min_p_single_effect = 0.8,
    use_null_weight = TRUE,
    min_snps = 2,
    min_genes = 1,
    ncore = 1,
    logfile = NULL,
    verbose = FALSE,
    ...){

  if (!is.null(logfile)){
    addHandler(writeToFile, file= logfile, level='DEBUG')
  }

  loginfo('Estimating parameters ... ')

  # check inputs
  if (!inherits(region_data, "list")){
    stop("'region_data' should be a list.")
  }

  # extract thin value from region_data
  thin <- unique(sapply(region_data, "[[", "thin"))
  if (length(thin) > 1) {
    thin <- min(thin)
    loginfo("thin has more than one value in region_data, use the minimum thin value.")
  }
  loginfo("thin = %s", thin)

  # adjust to account for thin argument
  if (!is.null(init_group_prior)){
    init_group_prior["SNP"] <- init_group_prior["SNP"]/thin
  }

  region_ids <- names(region_data)
  n_gids <- sapply(region_data, function(x){length(x$gid)})
  n_sids <- sapply(region_data, function(x){length(x$sid)})
  p_single_effect_df <- data.frame(region_id = region_ids,
                                   p_single_effect = NA)

  # skip regions with fewer than min_snps SNPs
  if (min_snps > 0) {
    skip_region_ids <- region_ids[n_sids < min_snps]
    if (length(skip_region_ids) > 0){
      loginfo("Skip %d regions with number of SNPs < %d.", length(skip_region_ids), min_snps)
      region_data[skip_region_ids] <- NULL
    }
  }

  # skip regions with fewer than min_genes genes
  if (min_genes > 0) {
    skip_region_ids <- region_ids[n_gids < min_genes]
    if (length(skip_region_ids) > 0){
      loginfo("Skip %d regions with number of genes < %d.", length(skip_region_ids), min_genes)
      region_data[skip_region_ids] <- NULL
    }
  }

  # Run EM for a few (niter_prefit) iterations, getting rough estimates
  loginfo("Run EM (prefit) for %d iterations, getting rough estimates ...", niter_prefit)
  EM_prefit_res <- fit_EM(region_data,
                          niter = niter_prefit,
                          init_group_prior = init_group_prior,
                          init_group_prior_var = init_group_prior_var,
                          use_null_weight = use_null_weight,
                          ncore = ncore,
                          verbose = verbose,
                          ...)
  adjusted_EM_prefit_group_prior <- EM_prefit_res$group_prior
  adjusted_EM_prefit_group_prior["SNP"] <- EM_prefit_res$group_prior["SNP"] * thin
  loginfo("Roughly estimated group_prior {%s}: {%s}",
          names(EM_prefit_res$group_prior), format(adjusted_EM_prefit_group_prior, digits = 4))
  loginfo("Roughly estimated group_prior_var {%s}: {%s}",
          names(EM_prefit_res$group_prior_var), format(EM_prefit_res$group_prior_var, digits = 4))
  group_size <- EM_prefit_res$group_size

  # Select regions with single effect
  p_single_effect <- compute_region_p_single_effect(region_data, EM_prefit_res$group_prior)
  selected_region_ids <- names(p_single_effect)[p_single_effect >= min_p_single_effect]
  loginfo("Selected %d regions with p(single effect) >= %s", length(selected_region_ids), min_p_single_effect)
  selected_region_data <- region_data[selected_region_ids]
  idx <- match(names(p_single_effect), p_single_effect_df$region_id)
  p_single_effect_df$p_single_effect[idx] <- p_single_effect

  # Run EM for more (niter) iterations, getting rough estimates
  loginfo("Run EM for %d iterations, getting accurate estimates ...", niter)
  EM_res <- fit_EM(selected_region_data,
                   niter = niter,
                   init_group_prior = EM_prefit_res$group_prior,
                   init_group_prior_var = EM_prefit_res$group_prior_var,
                   use_null_weight = use_null_weight,
                   ncore = ncore,
                   verbose = verbose,
                   ...)
  group_prior <- EM_res$group_prior
  group_prior_var <- EM_res$group_prior_var

  # record estimated parameters from all iterations
  group_prior_iters <- EM_res$group_prior_iters
  group_prior_var_iters <- EM_res$group_prior_var_iters

  # adjust parameters to account for thin
  if (thin < 1){
    group_prior["SNP"] <- group_prior["SNP"] * thin
    group_prior_iters["SNP",] <- group_prior_iters["SNP",] * thin
    group_size["SNP"] <- group_size["SNP"] / thin
  }
  group_size <- group_size[names(group_prior)]

  loginfo("Estimated group_prior {%s}: {%s}", names(group_prior), format(group_prior, digits = 4))
  loginfo("Estimated group_prior_var {%s}: {%s}", names(group_prior_var), format(group_prior_var, digits = 4))
  loginfo("group_size {%s}: {%s}", names(group_size), group_size)

  param <- list("group_prior" = group_prior,
                "group_prior_var" = group_prior_var,
                "group_prior_iters" = group_prior_iters,
                "group_prior_var_iters" = group_prior_var_iters,
                "group_size" = group_size,
                "p_single_effect" = p_single_effect_df)

  return(param)
}

