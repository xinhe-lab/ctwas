#' @title Estimates cTWAS parameters using EM with cTWAS SER model
#'
#' @param region_data a list object indexing regions, variants and genes.
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
#' @param niter_prefit the maximum number of iterations of the E-M algorithm to perform during the initial parameter estimation step.
#'
#' @param niter the maximum number of iterations of the E-M algorithm to perform during the complete parameter estimation step.
#'
#' @param min_var minimum number of variables (SNPs and genes) in a region.
#'
#' @param min_gene minimum number of genes in a region.
#'
#' @param min_group_size Minimum number of variables in a group.
#'
#' @param min_p_single_effect Regions with probability greater than \code{min_p_single_effect} of
#' having at most one causal effect will be used selected for the complete parameter estimation step.
#'
#' @param null_method Method to compute null model, options: "ctwas", "susie" or "none".
#'
#' @param EM_tol A small, non-negative number specifying the convergence
#'   tolerance of log-likelihood for the EM iterations.
#'
#' @param force_run_niter If TRUE, run all the \code{niter} EM iterations.
#'
#' @param enrichment_test Method to test enrichment,
#' options: "G" (G-test), "fisher" (Fisher's exact test).
#' Only used when \code{run_enrichment_test = TRUE}.
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
    group_prior_var_structure = c("shared_all", "shared_type", "shared_context", "shared_nonSNP", "independent", "fixed"),
    niter_prefit = 3,
    niter = 50,
    min_var = 2,
    min_gene = 1,
    min_group_size = 100,
    min_p_single_effect = 0.8,
    null_method = c("ctwas", "susie", "none"),
    enrichment_test = c("G", "fisher"),
    EM_tol = 1e-4,
    force_run_niter = FALSE,
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
  null_method <- match.arg(null_method)
  enrichment_test <- match.arg(enrichment_test)

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
                                   # n_gids = n_gids,
                                   # n_sids = n_sids,
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
  loginfo("Run EM (prefit) iterations, getting rough estimates ...")
  loginfo("group_prior_var_structure = '%s'", group_prior_var_structure)

  loginfo("Using data in %d regions", length(region_data))
  if (length(region_data) == 0){
    stop("No regions selected!")
  }

  EM_prefit_res <- fit_EM(region_data,
                          niter = niter_prefit,
                          init_group_prior = init_group_prior,
                          init_group_prior_var = init_group_prior_var,
                          group_prior_var_structure = group_prior_var_structure,
                          null_method = null_method,
                          EM_tol = EM_tol,
                          force_run_niter = force_run_niter,
                          warn_converge_fail = FALSE,
                          ncore = ncore,
                          verbose = verbose)
  adjusted_EM_prefit_group_prior <- EM_prefit_res$group_prior

  # adjust thin
  if (thin != 1){
    adjusted_EM_prefit_group_prior["SNP"] <- EM_prefit_res$group_prior["SNP"] * thin
  }

  loginfo("Roughly estimated group_prior {%s}: {%s}",
          names(EM_prefit_res$group_prior), format(adjusted_EM_prefit_group_prior, digits = 4))
  loginfo("Roughly estimated group_prior_var {%s}: {%s}",
          names(EM_prefit_res$group_prior_var), format(EM_prefit_res$group_prior_var, digits = 4))

  # Select regions with single effect
  p_single_effect <- compute_region_p_single_effect(region_data, EM_prefit_res$group_prior)
  selected_region_ids <- names(p_single_effect)[p_single_effect > min_p_single_effect]
  loginfo("Selected %d regions with p(single effect) > %s", length(selected_region_ids), min_p_single_effect)
  selected_region_data <- region_data[selected_region_ids]
  idx <- match(names(p_single_effect), p_single_effect_df$region_id)
  p_single_effect_df$p_single_effect[idx] <- p_single_effect
  rownames(p_single_effect_df) <- NULL

  # Run EM for more (niter) iterations, getting rough estimates
  loginfo("Run EM iterations, getting accurate estimates ...")
  loginfo("Using data in %d regions", length(selected_region_data))
  if (length(selected_region_data) == 0){
    stop("No regions selected!")
  }

  EM_res <- fit_EM(selected_region_data,
                   niter = niter,
                   init_group_prior = EM_prefit_res$group_prior,
                   init_group_prior_var = EM_prefit_res$group_prior_var,
                   group_prior_var_structure = group_prior_var_structure,
                   null_method = null_method,
                   EM_tol = EM_tol,
                   force_run_niter = force_run_niter,
                   warn_converge_fail = TRUE,
                   ncore = ncore,
                   verbose = verbose)
  group_prior <- EM_res$group_prior
  group_prior_var <- EM_res$group_prior_var
  group_prior_var_structure <- EM_res$group_prior_var_structure

  # record estimated parameters from all iterations
  group_prior_iters <- EM_res$group_prior_iters
  group_prior_var_iters <- EM_res$group_prior_var_iters

  if (anyNA(group_prior)) {
    stop("Estimated group_prior contains NAs!")
  }

  if (anyNA(group_prior_var)) {
    stop("Estimated group_prior_var contains NAs!")
  }

  group_size <- get_group_size_from_region_data(region_data)
  group_size <- group_size[names(group_prior)]

  # adjust parameters to account for thin
  if (thin != 1){
    group_prior["SNP"] <- group_prior["SNP"] * thin
    group_prior_iters["SNP",] <- group_prior_iters["SNP",] * thin
    group_size["SNP"] <- group_size["SNP"]/thin
  }

  loginfo("Estimated group_prior {%s}: {%s}", names(group_prior), format(group_prior, digits = 4))
  loginfo("Estimated group_prior_var {%s}: {%s}", names(group_prior_var), format(group_prior_var, digits = 4))
  loginfo("group_size {%s}: {%s}", names(group_size), group_size)

  param <- list("group_prior" = group_prior,
                "group_prior_var" = group_prior_var,
                "group_prior_iters" = group_prior_iters,
                "group_prior_var_iters" = group_prior_var_iters,
                "group_prior_var_structure" = group_prior_var_structure,
                "group_size" = group_size,
                "p_single_effect" = p_single_effect_df)

  # compute enrichment, s.e. and p-values
  enrichment_res <- compute_enrichment_test(group_prior = group_prior,
                                            group_size = group_size,
                                            enrichment_test = enrichment_test)

  loginfo("Estimated enrichment (log scale) {%s}: {%s}",
          names(enrichment_res$enrichment),
          format(enrichment_res$enrichment, digits = 4))

  param$enrichment = enrichment_res$enrichment
  param$enrichment_se = enrichment_res$se
  param$enrichment_pval = enrichment_res$p.value

  # include log-likelihood from all iterations
  param$loglik_iters <- EM_res$loglik_iters
  param$converged <- EM_res$converged
  param$niter <- EM_res$niter

  return(param)
}


#' @title Computes enrichment (log-scale), standard error and p-value
#'
#' @param group_prior a vector of prior inclusion probabilities for different groups.
#'
#' @param group_size a vector of number of variables in different groups.
#'
#' @param enrichment_test Method to test enrichment,
#'"G": G-test, "fisher": Fisher's exact test.
#'
#' @return Estimated enrichment, S.E. and p-value from G-test or Fisher's exact test.
#'
#' @importFrom AMR g.test
#' @importFrom stats fisher.test
#'
#' @export
#'
compute_enrichment_test <- function(group_prior,
                                    group_size,
                                    enrichment_test = c("G", "fisher")){

  enrichment_test <- match.arg(enrichment_test)

  # compute sum of PIPs for each group
  groups <- names(group_prior)
  group_size <- group_size[groups]
  group_pip <- group_prior * group_size
  names(group_pip) <- groups

  gene_groups <- setdiff(groups, "SNP")

  enrichment <- rep(NA, length(gene_groups))
  names(enrichment) <- gene_groups
  enrichment.se <- rep(NA, length(gene_groups))
  names(enrichment.se) <- gene_groups
  enrichment.pval <- rep(NA, length(gene_groups))
  names(enrichment.pval) <- gene_groups

  for (group in gene_groups) {
    k0 <- group_size["SNP"]
    k1 <- group_size[group]

    r0 <- group_pip["SNP"]
    r1 <- group_pip[group]

    # enrichment: relative risk estimate in a 2 × 2 contingency table
    enrichment[group] <- log((r1/k1) / (r0/k0))

    # standard error: standard error of a relative risk
    enrichment.se[group] <- sqrt(1/r1 + 1/r0 - 1/k1 - 1/k0)

    obs <- rbind(c(r1, k1), c(r0, k0))
    colnames(obs) <- c("r", "k")
    rownames(obs) <- c(group, "SNP")

    if (enrichment_test == "G"){
      # evaluate statistical significance of enrichment using G-test
      test_res <- g.test(obs)
      enrichment.pval[group] <- test_res$p.value
    } else if (enrichment_test == "fisher"){
      # evaluate statistical significance of enrichment using Fisher's exact test
      test_res <- fisher.test(round(obs))
      enrichment.pval[group] <- test_res$p.value
    }
  }

  return(list("enrichment" = enrichment,
              "se" = enrichment.se,
              "p.value" = enrichment.pval))
}
