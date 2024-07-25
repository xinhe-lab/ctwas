#' @title Estimates cTWAS parameters using EM
#'
#' @param region_data a list object indexing regions, variants and genes.
#'
#' @param init_group_prior a vector of initial values of prior inclusion probabilities for SNPs and genes.
#'
#' @param init_group_prior_var a vector of initial values of prior variances for SNPs and gene effects.
#'
#' @param niter_prefit the number of iterations of the E-M algorithm to perform during the initial parameter estimation step
#'
#' @param niter the number of iterations of the E-M algorithm to perform during the complete parameter estimation step
#'
#' @param p_single_effect Regions with probability greater than \code{p_single_effect} of
#' having at most one causal effect will be used selected for the complete parameter estimation step
#'
#' @param use_null_weight TRUE/FALSE. If TRUE, allow for a probability of no effect in susie
#'
#' @param minvar minimum number of variables (snps and genes) in a region
#'
#' @param mingene minimum number of genes in a region
#'
#' @param ncore The number of cores used to parallelize computation over regions
#'
#' @param logfile the log file, if NULL will print log info on screen
#'
#' @param verbose TRUE/FALSE. If TRUE, print detail messages
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
    p_single_effect = 0.8,
    use_null_weight = TRUE,
    minvar = 2,
    mingene = 1,
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

  # remove regions with fewer than minvar variables
  if (minvar > 0) {
    n.var <- sapply(region_data, function(x){length(x$gid) + length(x$sid)})
    drop.idx <- which(n.var < minvar)
    if (length(drop.idx) > 0){
      loginfo("Remove %d regions with number of variables < %d.", length(drop.idx), minvar)
      region_data[drop.idx] <- NULL
    }
  }

  # remove regions with fewer than mingene genes
  if (mingene > 0) {
    n.gid <- sapply(region_data, function(x){length(x$gid)})
    drop.idx <- which(n.gid < mingene)
    if (length(drop.idx) > 0){
      loginfo("Remove %d regions with number of genes < %d.", length(drop.idx), mingene)
      region_data[drop.idx] <- NULL
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
  loginfo("Roughly estimated group_prior {%s}: {%s} (thin = %s)", names(EM_prefit_res$group_prior), format(EM_prefit_res$group_prior, digits = 4), thin)
  loginfo("Roughly estimated group_prior_var {%s}: {%s}", names(EM_prefit_res$group_prior_var), format(EM_prefit_res$group_prior_var, digits = 4))
  group_size <- EM_prefit_res$group_size

  # Select regions with single effect
  region_single_effect_df <- compute_region_p_single_effect(region_data, EM_prefit_res$group_prior)
  region_single_effect_df <- region_single_effect_df[region_single_effect_df$p_single_effect >= p_single_effect, ]
  selected_region_ids <- region_single_effect_df$region_id
  loginfo("Selected %d regions with P(single effect) >= %s", length(selected_region_ids), p_single_effect)
  selected_region_data <- region_data[selected_region_ids]

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
  group_prior_var_structure <- EM_res$group_prior_var_structure
  loginfo("Estimated group_prior {%s}: {%s} (thin = %s)", names(group_prior), format(group_prior, digits = 4), thin)
  loginfo("Estimated group_prior_var {%s}: {%s}", names(group_prior_var), format(group_prior_var, digits = 4))

  # record estimated parameters from all iterations
  group_prior_iters <- EM_res$group_prior_iters
  group_prior_var_iters <- EM_res$group_prior_var_iters

  # adjust parameters to account for thin
  if (thin != 1){
    loginfo("Adjust parameters to account for thin (thin = %s)", thin)
    group_prior["SNP"] <- group_prior["SNP"] * thin
    group_prior_iters["SNP",] <- group_prior_iters["SNP",] * thin
    group_size["SNP"] <- group_size["SNP"] / thin
  }
  group_size <- group_size[names(group_prior)]
  loginfo("Estimated group_prior {%s}: {%s}", names(group_prior), format(group_prior, digits = 4))
  loginfo("Estimated group_prior_var {%s}: {%s}", names(group_prior_var), format(group_prior_var, digits = 4))
  loginfo("Estimated group_size {%s}: {%s}", names(group_size), group_size)

  param <- list("group_prior" = group_prior,
                "group_prior_var" = group_prior_var,
                "group_prior_iters" = group_prior_iters,
                "group_prior_var_iters" = group_prior_var_iters,
                "group_size" = group_size)

  return(param)
}


# Select single effect regions
compute_region_p_single_effect <- function(region_data, group_prior){
  region_ids <- names(region_data)
  if (length(region_ids) == 0)
    stop("No region_ids in region_data!")

  p_single_effect <- sapply(region_ids, function(region_id){
    res <- extract_region_data(region_data, region_id)
    gs_group <- res$gs_group
    group_size <- table(gs_group)[names(group_prior)]
    group_size[is.na(group_size)] <- 0
    p1 <- prod((1-group_prior)^group_size) * (1 + sum(group_size*(group_prior/(1-group_prior))))
    p1
  })
  return(data.frame(region_id = region_ids, p_single_effect = p_single_effect))

}

