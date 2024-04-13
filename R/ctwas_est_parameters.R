#' Estimate cTWAS parameters
#'
#' @param region_data a list object indexing regions, variants and genes.
#'
#' @param init_group_prior a vector of initial values of prior inclusion probabilities for SNPs and genes.
#'
#' @param init_group_prior_var a vector of initial values of prior variances for SNPs and gene effects.
#'
#' @param group_prior_var_structure a string indicating the structure to put on the prior variance parameters.
#' "independent" is the default and allows all groups to have their own separate variance parameters.
#' "shared_all" allows all groups to share the same variance parameter.
#' "shared_type" allows all groups in one molecular QTL type to share the same variance parameter.
#'
#' @param niter1 the number of iterations of the E-M algorithm to perform during the initial parameter estimation step
#'
#' @param niter2 the number of iterations of the E-M algorithm to perform during the complete parameter estimation step
#'
#' @param p_single_effect Regions with probability greater than \code{p_single_effect} of
#' having at most one causal effect will be used selected for the complete parameter estimation step
#'
#' @param ncore The number of cores used to parallelize computation over regions
#'
#' @param logfile the log file, if NULL will print log info on screen
#'
#' @param verbose TRUE/FALSE. If TRUE, print detail messages
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
    group_prior_var_structure = c("independent","shared_all","shared_type"),
    niter1 = 3,
    niter2 = 30,
    p_single_effect = 0.8,
    ncore = 1,
    logfile = NULL,
    verbose = FALSE){

  if (!is.null(logfile)){
    addHandler(writeToFile, file= logfile, level='DEBUG')
  }

  loginfo('Estimating parameters ... ')

  group_prior_var_structure <- match.arg(group_prior_var_structure)

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

  # Run EM for a few (niter1) iterations, getting rough estimates
  loginfo("Run EM for %d iterations on %d regions, getting rough estimates ...",
          niter1, length(region_data))
  EM1_res <- EM(region_data,
                niter = niter1,
                init_group_prior = init_group_prior,
                init_group_prior_var = init_group_prior_var,
                group_prior_var_structure = group_prior_var_structure,
                ncore = ncore,
                verbose = verbose)
  loginfo("Roughly estimated group_prior {%s}: {%s}", names(EM1_res$group_prior), EM1_res$group_prior)
  loginfo("Roughly estimated group_prior_var {%s}: {%s}", names(EM1_res$group_prior_var), EM1_res$group_prior_var)

  # Select regions with single effect
  region_p1_df <- compute_region_p_single_effect(region_data, EM1_res$group_prior)
  selected_region_ids <- region_p1_df[region_p1_df$p1 >= p_single_effect, "region_id"]
  loginfo("Selected %d regions with P(single effect) >= %s", length(selected_region_ids), p_single_effect)
  selected_region_data <- region_data[selected_region_ids]

  # Run EM for more (niter2) iterations, getting rough estimates
  loginfo("Run EM for %d iterations on %d regions, getting accurate estimates ...",
          niter2, length(selected_region_data))

  EM2_res <- EM(selected_region_data,
                niter = niter2,
                init_group_prior = EM1_res$group_prior,
                init_group_prior_var = EM1_res$group_prior_var,
                group_prior_var_structure = group_prior_var_structure,
                ncore = ncore,
                verbose = verbose)
  group_prior <- EM2_res$group_prior
  group_prior_var <- EM2_res$group_prior_var
  group_prior_var_structure <- EM2_res$group_prior_var_structure
  loginfo("Estimated group_prior {%s}: {%s}", names(group_prior), group_prior)
  loginfo("Estimated group_prior_var {%s}: {%s}", names(group_prior_var), group_prior_var)

  # estimated parameters from all iterations
  group_prior_iters <- EM2_res$group_prior_iters
  group_prior_var_iters <- EM2_res$group_prior_var_iters

  # adjust parameters to account for thin argument
  group_prior["SNP"] <- group_prior["SNP"] * thin
  group_prior_iters["SNP",] <- group_prior_iters["SNP",] * thin

  group_size <- EM1_res$group_size
  group_size["SNP"] <- group_size["SNP"]/thin

  param <- list("group_prior" = group_prior,
                "group_prior_var" = group_prior_var,
                "group_prior_iters" = group_prior_iters,
                "group_prior_var_iters" = group_prior_var_iters,
                "group_prior_var_structure" = group_prior_var_structure,
                "group_size" = group_size)

  return(param)
}


