#' Estimate cTWAS parameters
#'
#' @param regionlist a list object indexing regions, variants and genes.
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
#' @param thin The proportion of SNPs to be used for the parameter estimation and
#' initial screening region steps.
#' Smaller \code{thin} parameters reduce runtime at the expense of accuracy.
#' The fine mapping step is rerun using full SNPs for regions with strong gene signals.
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
    regionlist,
    init_group_prior = NULL,
    init_group_prior_var = NULL,
    group_prior_var_structure = c("independent","shared_all","shared_type"),
    thin = 1,
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

  if (thin <= 0 | thin > 1){
    stop("thin value needs to be in (0,1]")
  }

  if (!is.null(init_group_prior)){
    init_group_prior["SNP"] <- init_group_prior["SNP"]/thin # adjust to account for thin argument
  }

  # Run EM for a few (niter1) iterations, getting rough estimates
  loginfo("Run EM for %d iterations on %d regions, getting rough estimates ...",
          niter1, length(regionlist))
  EM1_res <- EM(regionlist,
                niter = niter1,
                init_group_prior = init_group_prior,
                init_group_prior_var = init_group_prior_var,
                group_prior_var_structure = group_prior_var_structure,
                ncore = ncore,
                verbose = verbose)
  loginfo("Roughly estimated group_prior {%s}: {%s}", names(EM1_res$group_prior), EM1_res$group_prior)
  loginfo("Roughly estimated group_prior_var {%s}: {%s}", names(EM1_res$group_prior_var), EM1_res$group_prior_var)

  # Select regions with single effect
  selected_region_tags <- select_single_effect_regions(regionlist, EM1_res$group_prior, p_single_effect)
  loginfo("Selected %d regions with P(single effect) >= %s", length(selected_region_tags), p_single_effect)
  selected_regionlist <- regionlist[selected_region_tags]

  # Run EM for more (niter2) iterations, getting rough estimates
  loginfo("Run EM for %d iterations on %d regions, getting accurate estimates ...",
          niter2, length(selected_regionlist))

  EM2_res <- EM(selected_regionlist,
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

  # estimated prior records (all iterations)
  group_prior_rec <- cbind(EM1_res$group_prior_rec, EM2_res$group_prior_rec)
  colnames(group_prior_rec) <- c(paste0("EM1_iter", 1:ncol(EM1_res$group_prior_rec)), paste0("EM2_iter", 1:ncol(EM2_res$group_prior_rec)))

  group_prior_var_rec <- cbind(EM1_res$group_prior_var_rec, EM2_res$group_prior_var_rec)
  colnames(group_prior_var_rec) <- c(paste0("EM1_iter", 1:ncol(EM1_res$group_prior_var_rec)), paste0("EM2_iter", 1:ncol(EM2_res$group_prior_var_rec)))

  # adjust parameters to account for thin argument
  group_prior["SNP"] <- group_prior["SNP"] * thin
  group_prior_rec["SNP",] <- group_prior_rec["SNP",] * thin
  group_size <- EM1_res$group_size
  group_size["SNP"] <- group_size["SNP"]/thin

  param <- list("group_prior" = group_prior,
                "group_prior_var" = group_prior_var,
                "group_prior_rec" = group_prior_rec,
                "group_prior_var_rec" = group_prior_var_rec,
                "group_prior_var_structure" = group_prior_var_structure,
                "group_size" = group_size)

  return(param)
}


