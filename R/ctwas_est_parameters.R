#' @title Estimates cTWAS parameters using EM
#'
#' @param region_data a list object indexing regions, variants and genes.
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
#' @param niter_prefit the number of iterations of the E-M algorithm to perform during the initial parameter estimation step.
#'
#' @param niter the number of iterations of the E-M algorithm to perform during the complete parameter estimation step.
#'
#' @param min_p_single_effect Regions with probability greater than \code{min_p_single_effect} of
#' having at most one causal effect will be used selected for the complete parameter estimation step.
#'
#' @param use_null_weight If TRUE, allow for a probability of no effect in susie.
#'
#' @param min_var minimum number of variables (SNPs and genes) in a region.
#'
#' @param min_gene minimum number of genes in a region.
#'
#' @param min_group_size Minimum number of variables in a group.
#'
#' @param ncore The number of cores used to parallelize computation over regions.
#'
#' @param logfile The log filename. If NULL, print log info on screen.
#'
#' @param verbose If TRUE, print detail messages.
#'
#' @param ... Additional arguments of \code{susie_rss}.
#'
#' @importFrom logging addHandler loginfo logwarn writeToFile
#'
#' @return a list with estimated parameters
#'
#' @export
#'
est_param <- function(
    region_data,
    init_group_prior = NULL,
    init_group_prior_var = NULL,
    group_prior_var_structure = c("shared_type", "shared_context", "shared_nonSNP", "shared_all", "independent"),
    niter_prefit = 3,
    niter = 30,
    min_p_single_effect = 0.8,
    use_null_weight = TRUE,
    min_var = 2,
    min_gene = 1,
    min_group_size = 100,
    ncore = 1,
    logfile = NULL,
    verbose = FALSE,
    ...){

  if (!is.null(logfile)){
    addHandler(writeToFile, file=logfile, level='DEBUG')
  }

  loginfo('Estimating parameters ... ')

  # check inputs
  group_prior_var_structure <- match.arg(group_prior_var_structure)

  if (!inherits(region_data,"list"))
    stop("'region_data' should be a list.")

  if (anyDuplicated(names(region_data)))
    logwarn("Duplicated names of region_data found! Please use unique names for region_data!")

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

  # filter groups with number of variables < min_group_size
  group_size_thinned <- get_group_size_from_region_data(region_data)
  if (any(group_size_thinned < min_group_size)){
    group_size_drop <- group_size_thinned[group_size_thinned < min_group_size]
    logwarn("Groups with group size < %d:\n {%s}: {%s}", min_group_size, names(group_size_drop), group_size_drop)
    stop("Parameters may not be reliable for groups with too few variables! Remove those groups from z_gene!")
  }

  region_ids <- names(region_data)
  n_gids <- sapply(region_data, function(x){length(x$gid)})
  n_sids <- sapply(region_data, function(x){length(x$sid)})
  p_single_effect_df <- data.frame(region_id = region_ids,
                                   p_single_effect = NA)

  # skip regions with fewer than min_var variables
  if (min_var > 0) {
    skip_region_ids <- region_ids[(n_sids + n_gids) < min_var]
    if (length(skip_region_ids) > 0){
      loginfo("Skip %d regions with number of variables < %d", length(skip_region_ids), min_var)
      region_data <- region_data[!names(region_data) %in% skip_region_ids]
    }
  }

  # skip regions with fewer than min_gene genes
  if (min_gene > 0) {
    skip_region_ids <- region_ids[n_gids < min_gene]
    if (length(skip_region_ids) > 0){
      loginfo("Skip %d regions with number of genes < %d", length(skip_region_ids), min_gene)
      region_data <- region_data[!names(region_data) %in% skip_region_ids]
    }
  }

  # Run EM for a few (niter_prefit) iterations, getting rough estimates
  loginfo("Run EM (prefit) for %d iterations, getting rough estimates ...", niter_prefit)
  loginfo("Using data in %d regions", length(region_data))
  if (length(region_data) == 0){
    stop("No regions selected!")
  }
  EM_prefit_res <- fit_EM(region_data,
                          niter = niter_prefit,
                          init_group_prior = init_group_prior,
                          init_group_prior_var = init_group_prior_var,
                          group_prior_var_structure = group_prior_var_structure,
                          use_null_weight = use_null_weight,
                          ncore = ncore,
                          verbose = verbose,
                          ...)
  adjusted_EM_prefit_group_prior <- EM_prefit_res$group_prior
  group_size <- EM_prefit_res$group_size
  # adjust thin
  if (thin != 1){
    adjusted_EM_prefit_group_prior["SNP"] <- EM_prefit_res$group_prior["SNP"] * thin
    group_size["SNP"] <- group_size["SNP"]/thin
  }
  group_size <- group_size[names(EM_prefit_res$group_prior)]
  loginfo("group_size {%s}: {%s}", names(group_size), group_size)

  loginfo("Roughly estimated group_prior {%s}: {%s}",
          names(EM_prefit_res$group_prior), format(adjusted_EM_prefit_group_prior, digits = 4))
  loginfo("Roughly estimated group_prior_var {%s}: {%s}",
          names(EM_prefit_res$group_prior_var), format(EM_prefit_res$group_prior_var, digits = 4))

  # Select regions with single effect
  p_single_effect <- compute_region_p_single_effect(region_data, EM_prefit_res$group_prior)
  selected_region_ids <- names(p_single_effect)[p_single_effect >= min_p_single_effect]
  loginfo("Selected %d regions with p(single effect) >= %s", length(selected_region_ids), min_p_single_effect)
  selected_region_data <- region_data[selected_region_ids]
  idx <- match(names(p_single_effect), p_single_effect_df$region_id)
  p_single_effect_df$p_single_effect[idx] <- p_single_effect

  # Run EM for more (niter) iterations, getting rough estimates
  loginfo("Run EM for %d iterations, getting accurate estimates ...", niter)
  loginfo("Using data in %d regions", length(selected_region_data))
  if (length(selected_region_data) == 0){
    stop("No regions selected!")
  }
  EM_res <- fit_EM(selected_region_data,
                   niter = niter,
                   init_group_prior = EM_prefit_res$group_prior,
                   init_group_prior_var = EM_prefit_res$group_prior_var,
                   group_prior_var_structure = group_prior_var_structure,
                   use_null_weight = use_null_weight,
                   ncore = ncore,
                   verbose = verbose,
                   ...)
  group_prior <- EM_res$group_prior
  group_prior_var <- EM_res$group_prior_var
  group_prior_var_structure <- EM_res$group_prior_var_structure

  # record estimated parameters from all iterations
  group_prior_iters <- EM_res$group_prior_iters
  group_prior_var_iters <- EM_res$group_prior_var_iters

  # adjust parameters to account for thin
  if (thin != 1){
    group_prior["SNP"] <- group_prior["SNP"] * thin
    group_prior_iters["SNP",] <- group_prior_iters["SNP",] * thin
  }

  loginfo("Estimated group_prior {%s}: {%s}", names(group_prior), format(group_prior, digits = 4))
  loginfo("Estimated group_prior_var {%s}: {%s}", names(group_prior_var), format(group_prior_var, digits = 4))

  if (anyNA(group_prior)) {
    stop("Estimated group_prior contains NAs!")
  }

  if (anyNA(group_prior_var)) {
    stop("Estimated group_prior_var contains NAs!")
  }

  rownames(p_single_effect_df) <- NULL

  param <- list("group_prior" = group_prior,
                "group_prior_var" = group_prior_var,
                "group_prior_iters" = group_prior_iters,
                "group_prior_var_iters" = group_prior_var_iters,
                "group_prior_var_structure" = group_prior_var_structure,
                "group_size" = group_size,
                "p_single_effect" = p_single_effect_df)

  return(param)
}
