#' @title Runs cTWAS fine-mapping for regions
#'
#' @param region_data region_data to be finemapped
#'
#' @param LD_map a data frame with filenames of LD matrices for the regions.
#'
#' @param weights a list of preprocessed weights.
#'
#' @param L the number of effects or a vector of number of effects for each region.
#'
#' @param group_prior a vector of two prior inclusion probabilities for SNPs and genes.
#' If NULL, it will use uniform prior inclusion probabilities.
#'
#' @param group_prior_var a vector of two prior variances for SNPs and gene effects.
#' If NULL, it will set prior variance = 50 as the default in \code{susie_rss}.
#'
#' @param use_null_weight If TRUE, allow for a probability of no effect in susie
#'
#' @param coverage A number between 0 and 1 specifying the \dQuote{coverage} of the estimated confidence sets
#'
#' @param min_abs_corr Minimum absolute correlation allowed in a credible set.
#'
#' @param include_cs If TRUE, add credible sets (CS) to finemapping results.
#'
#' @param get_susie_alpha If TRUE, get susie alpha matrix from finemapping results.
#'
#' @param snps_only If TRUE, use only SNPs in the region data.
#'
#' @param force_compute_cor If TRUE, force computing correlation (R) matrices
#'
#' @param save_cor If TRUE, save correlation (R) matrices to \code{cor_dir}
#'
#' @param cor_dir a string, the directory to store correlation (R) matrices
#'
#' @param LD_format file format for LD matrix. If "custom", use a user defined
#' \code{LD_loader_fun()} function to load LD matrix.
#'
#' @param LD_loader_fun a user defined function to load LD matrix when \code{LD_format = "custom"}.
#'
#' @param snpinfo_loader_fun a user defined function to load SNP information file,
#' if SNP information files are not in standard cTWAS reference format.
#'
#' @param ncore The number of cores used to parallelize computation over regions
#'
#' @param verbose If TRUE, print detail messages
#'
#' @param logfile the log file, if NULL will print log info on screen
#'
#' @param ... Additional arguments of \code{susie_rss}.
#'
#' @return a list with cTWAS finemapping results.
#'
#' @importFrom logging addHandler loginfo writeToFile
#' @importFrom parallel mclapply
#'
#' @export
#'
finemap_regions <- function(region_data,
                            LD_map,
                            weights,
                            L = 5,
                            group_prior = NULL,
                            group_prior_var = NULL,
                            use_null_weight = TRUE,
                            coverage = 0.95,
                            min_abs_corr = 0.1,
                            include_cs = TRUE,
                            get_susie_alpha = TRUE,
                            snps_only = FALSE,
                            force_compute_cor = FALSE,
                            save_cor = FALSE,
                            cor_dir = NULL,
                            LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
                            LD_loader_fun = NULL,
                            snpinfo_loader_fun = NULL,
                            ncore = 1,
                            verbose = FALSE,
                            logfile = NULL,
                            ...){

  if (!is.null(logfile)){
    addHandler(writeToFile, file=logfile, level='DEBUG')
  }

  loginfo("Fine-mapping %d regions ...", length(region_data))

  # check inputs
  LD_format <- match.arg(LD_format)

  if (!inherits(region_data,"list"))
    stop("'region_data' should be a list.")

  if (anyDuplicated(names(region_data)))
    logwarn("Duplicated names of region_data found! Please use unique names for region_data!")

  if (missing(LD_map))
    stop("'LD_map' is required when running with LD!")

  if (!inherits(LD_map,"data.frame"))
    stop("'LD_map' should be a data frame!")

  if (missing(weights))
    stop("'weights' is required when running with LD!")

  if (!inherits(weights,"list"))
    stop("'weights' should be a list!")

  if (any(sapply(weights, is.null)))
    stop("'weights' contain NULL, remove empty weights!")

  if (length(L) == 1) {
    L <- rep(L, length(region_data))
    names(L) <- names(region_data)
  } else if (length(L) == length(region_data)) {
    if (!all.equal(names(L), names(region_data)))
      stop("The names of L do not match with region_data!")
  } else {
    stop("L needs to an integer or a vector of the same length as region_data!")
  }

  # check if all groups have group_prior and group_prior_var
  groups <- unique(unlist(lapply(region_data, "[[", "groups")))
  if (!is.null(group_prior)) {
    groups_without_prior <- setdiff(groups, names(group_prior))
    if (length(groups_without_prior) > 0) {
      stop(paste("Missing group_prior for group:", groups_without_prior, "!"))
    }
  }

  if (!is.null(group_prior_var)) {
    groups_without_prior_var <- setdiff(groups, names(group_prior))
    if (length(groups_without_prior_var) > 0) {
      stop(paste("Missing group_prior_var for group:", groups_without_prior_var, "!"))
    }
  }

  if (verbose) {
    if (is.null(group_prior)) {
      loginfo("Use uniform prior")
    }
    loginfo("coverage = %s", coverage)
    loginfo("min_abs_corr = %s", min_abs_corr)
  }

  region_ids <- names(region_data)
  res <- mclapply_check(region_ids, function(region_id){
    finemap_single_region(region_data = region_data,
                          region_id = region_id,
                          LD_map = LD_map,
                          weights = weights,
                          L = L[region_id],
                          group_prior = group_prior,
                          group_prior_var = group_prior_var,
                          use_null_weight = use_null_weight,
                          coverage = coverage,
                          min_abs_corr = min_abs_corr,
                          include_cs = include_cs,
                          get_susie_alpha = get_susie_alpha,
                          snps_only = snps_only,
                          force_compute_cor = force_compute_cor,
                          save_cor = save_cor,
                          cor_dir = cor_dir,
                          LD_format = LD_format,
                          LD_loader_fun = LD_loader_fun,
                          snpinfo_loader_fun = snpinfo_loader_fun,
                          verbose = verbose,
                          ...)
  }, mc.cores = ncore, stop_if_missing = TRUE)

  finemap_res <- do.call(rbind, lapply(res, "[[", 1))
  rownames(finemap_res) <- NULL

  if (get_susie_alpha) {
    susie_alpha_res <- do.call(rbind, lapply(res, "[[", 2))
    rownames(susie_alpha_res) <- NULL
  } else {
    susie_alpha_res <- NULL
  }

  return(list("finemap_res" = finemap_res,
              "susie_alpha_res" = susie_alpha_res))
}

#' @title Runs cTWAS fine-mapping for regions without LD (L = 1)
#'
#' @param region_data region_data to be finemapped
#'
#' @param group_prior a vector of two prior inclusion probabilities for SNPs and genes.
#' If NULL, it will use uniform prior inclusion probabilities.
#'
#' @param group_prior_var a vector of two prior variances for SNPs and gene effects.
#' If NULL, it will set prior variance = 50 as the default in \code{susie_rss}.
#'
#' @param use_null_weight If TRUE, allow for a probability of no effect in susie
#'
#' @param get_susie_alpha If TRUE, get susie alpha matrix from finemapping results.
#'
#' @param snps_only If TRUE, use only SNPs in the region data.
#'
#' @param ncore The number of cores used to parallelize computation over regions
#'
#' @param verbose If TRUE, print detail messages
#'
#' @param logfile the log file, if NULL will print log info on screen
#'
#' @param ... Additional arguments of \code{susie_rss}.
#'
#' @return a list with cTWAS finemapping results.
#'
#' @importFrom logging addHandler loginfo logwarn writeToFile
#' @importFrom parallel mclapply
#'
#' @export
#'
finemap_regions_noLD <- function(region_data,
                                 group_prior = NULL,
                                 group_prior_var = NULL,
                                 use_null_weight = TRUE,
                                 get_susie_alpha = TRUE,
                                 snps_only = FALSE,
                                 ncore = 1,
                                 verbose = FALSE,
                                 logfile = NULL,
                                 ...){

  if (!is.null(logfile)){
    addHandler(writeToFile, file=logfile, level='DEBUG')
  }

  loginfo("Fine-mapping %d regions without LD ...", length(region_data))

  if (!inherits(region_data,"list"))
    stop("'region_data' should be a list!")

  if (anyDuplicated(names(region_data)))
    logwarn("Duplicated names of region_data found! Please use unique names for region_data!")

  # check if all groups have group_prior and group_prior_var
  groups <- unique(unlist(lapply(region_data, "[[", "groups")))
  if (!is.null(group_prior)) {
    groups_without_prior <- setdiff(groups, names(group_prior))
    if (length(groups_without_prior) > 0) {
      stop(paste("Missing group_prior for group:", groups_without_prior, "!"))
    }
  }

  if (!is.null(group_prior_var)) {
    groups_without_prior_var <- setdiff(groups, names(group_prior))
    if (length(groups_without_prior_var) > 0) {
      stop(paste("Missing group_prior_var for group:", groups_without_prior_var, "!"))
    }
  }

  if (verbose) {
    if (is.null(group_prior)) {
      loginfo("Use uniform prior")
    }
  }

  region_ids <- names(region_data)
  res <- mclapply_check(region_ids, function(region_id){
    finemap_single_region_noLD(region_data = region_data,
                               region_id = region_id,
                               group_prior = group_prior,
                               group_prior_var = group_prior_var,
                               use_null_weight = use_null_weight,
                               get_susie_alpha = get_susie_alpha,
                               snps_only = snps_only,
                               verbose = verbose,
                               ...)
  }, mc.cores = ncore, stop_if_missing = TRUE)

  finemap_res <- do.call(rbind, lapply(res, "[[", 1))
  rownames(finemap_res) <- NULL

  if (get_susie_alpha) {
    susie_alpha_res <- do.call(rbind, lapply(res, "[[", 2))
    rownames(susie_alpha_res) <- NULL
  } else {
    susie_alpha_res <- NULL
  }

  return(list("finemap_res" = finemap_res,
              "susie_alpha_res" = susie_alpha_res))
}

# Runs cTWAS finemapping for a single region
finemap_single_region <- function(region_data,
                                  region_id,
                                  LD_map,
                                  weights,
                                  L = 5,
                                  group_prior = NULL,
                                  group_prior_var = NULL,
                                  use_null_weight = TRUE,
                                  coverage = 0.95,
                                  min_abs_corr = 0.1,
                                  include_cs = TRUE,
                                  get_susie_alpha = TRUE,
                                  snps_only = FALSE,
                                  force_compute_cor = FALSE,
                                  save_cor = FALSE,
                                  cor_dir = NULL,
                                  LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
                                  LD_loader_fun = NULL,
                                  snpinfo_loader_fun = NULL,
                                  verbose = FALSE,
                                  ...){

  if (verbose){
    loginfo("Fine-mapping region %s", region_id)
  }

  # check inputs
  LD_format <- match.arg(LD_format)

  if (!inherits(region_data,"list"))
    stop("'region_data' should be a list.")

  if (!inherits(LD_map,"data.frame"))
    stop("'LD_map' should be a data frame")

  if (!inherits(weights,"list"))
    stop("'weights' should be a list.")

  if (any(sapply(weights, is.null)))
    stop("'weights' contain NULL, remove empty weights!")

  # load input data for the region
  regiondata <- extract_region_data(region_data, region_id,
                                    snps_only = snps_only)
  gids <- regiondata[["gid"]]
  sids <- regiondata[["sid"]]
  z <- regiondata[["z"]]
  gs_group <- regiondata[["gs_group"]]
  g_type <- regiondata[["g_type"]]
  g_context <- regiondata[["g_context"]]
  g_group <- regiondata[["g_group"]]
  groups <- regiondata$groups
  rm(regiondata)

  if (verbose){
    loginfo("%d genes, %d SNPs in the region", length(gids), length(sids))
  }

  if (length(z) < 2) {
    stop(paste(length(z), "variables in the region. At least two variables in a region are needed to run susie"))
  }

  if (!is.null(group_prior)) {
    groups_without_prior <- setdiff(groups, names(group_prior))
    if (length(groups_without_prior) > 0) {
      stop(paste("Missing group_prior for group:", groups_without_prior, "!"))
    }
  }
  if (!is.null(group_prior_var)) {
    groups_without_prior_var <- setdiff(groups, names(group_prior))
    if (length(groups_without_prior_var) > 0) {
      stop(paste("Missing group_prior_var for group:", groups_without_prior_var, "!"))
    }
  }

  res <- initiate_group_priors(group_prior[groups], group_prior_var[groups], groups)
  pi_prior <- res$pi_prior
  V_prior <- res$V_prior
  rm(res)

  # set prior and prior variance values for the region
  res <- set_region_susie_priors(pi_prior, V_prior, gs_group, L = L, use_null_weight = use_null_weight)
  prior <- res$prior
  V <- res$V
  null_weight <- res$null_weight
  rm(res)

  cor_res <- get_region_cor(region_id,
                            sids = sids,
                            gids = gids,
                            LD_map = LD_map,
                            weights = weights,
                            force_compute_cor = force_compute_cor,
                            save_cor = save_cor,
                            cor_dir = cor_dir,
                            LD_format = LD_format,
                            LD_loader_fun = LD_loader_fun,
                            snpinfo_loader_fun = snpinfo_loader_fun)

  # gene first then SNPs
  R <- rbind(cbind(cor_res$R_gene, t(cor_res$R_snp_gene)),
             cbind(cor_res$R_snp_gene, cor_res$R_snp))
  rm(cor_res)

  if (anyNA(R))
    stop("R matrix contains missing values!")

  if (length(z) != nrow(R))
    stop("R matrix dimension does not match with z!")

  # run susie for this region
  # in susie, prior_variance is under standardized scale (if performed)
  if (verbose){
    loginfo("run susie for region %s with L = %d", region_id, L)
  }
  susie_res <- ctwas_susie_rss(z = z,
                               R = R,
                               prior_weights = prior,
                               prior_variance = V,
                               L = L,
                               null_weight = null_weight,
                               coverage = coverage,
                               min_abs_corr = min_abs_corr,
                               ...)

  susie_res_df <- anno_susie(susie_res,
                             gids = gids,
                             sids = sids,
                             g_type = g_type,
                             g_context = g_context,
                             g_group = g_group,
                             region_id = region_id,
                             z = z,
                             include_cs = include_cs)

  if (get_susie_alpha) {
    # extract alpha matrix from susie result
    susie_alpha_df <- get_susie_alpha_res(susie_res, susie_res_df, keep_genes_only = TRUE)
  } else {
    susie_alpha_df <- NULL
  }

  return(list("susie_res_df" = susie_res_df,
              "susie_alpha_df" = susie_alpha_df))

}


# Runs cTWAS finemapping for a single region without LD
finemap_single_region_noLD <- function(region_data,
                                       region_id,
                                       group_prior = NULL,
                                       group_prior_var = NULL,
                                       use_null_weight = TRUE,
                                       coverage = 0.95,
                                       include_cs = TRUE,
                                       get_susie_alpha = TRUE,
                                       snps_only = FALSE,
                                       verbose = FALSE,
                                       ...){

  if (verbose){
    loginfo("Fine-mapping region %s without LD", region_id)
  }

  if (!inherits(region_data,"list"))
    stop("'region_data' should be a list.")

  regiondata <- extract_region_data(region_data, region_id,
                                    snps_only = snps_only)
  gids <- regiondata[["gid"]]
  sids <- regiondata[["sid"]]
  z <- regiondata[["z"]]
  gs_group <- regiondata[["gs_group"]]
  g_type <- regiondata[["g_type"]]
  g_context <- regiondata[["g_context"]]
  g_group <- regiondata[["g_group"]]
  groups <- regiondata$groups
  rm(regiondata)

  if (verbose){
    loginfo("%d genes, %d SNPs in the region", length(gids), length(sids))
  }

  if (length(z) < 2) {
    stop(paste(length(z), "variables in the region. At least two variables in a region are needed to run susie"))
  }

  if (!is.null(group_prior)) {
    groups_without_prior <- setdiff(groups, names(group_prior))
    if (length(groups_without_prior) > 0) {
      stop(paste("Missing group_prior for group:", groups_without_prior, "!"))
    }
  }
  if (!is.null(group_prior_var)) {
    groups_without_prior_var <- setdiff(groups, names(group_prior))
    if (length(groups_without_prior_var) > 0) {
      stop(paste("Missing group_prior_var for group:", groups_without_prior_var, "!"))
    }
  }

  res <- initiate_group_priors(group_prior[groups], group_prior_var[groups], groups)
  pi_prior <- res$pi_prior
  V_prior <- res$V_prior
  rm(res)

  # set prior and prior variance values for the region
  res <- set_region_susie_priors(pi_prior, V_prior, gs_group, L = 1, use_null_weight = use_null_weight)
  prior <- res$prior
  V <- res$V
  null_weight <- res$null_weight
  rm(res)

  # use an identity matrix as R in no-LD version
  R <- diag(length(z))

  if (anyNA(R))
    stop("R matrix contains missing values!")

  if (length(z) != nrow(R))
    stop("R matrix dimension does not match with z!")

  # run susie for this region
  # in susie, prior_variance is under standardized scale (if performed)
  if (verbose){
    loginfo("run susie for region %s with L = 1", region_id)
  }
  susie_res <- ctwas_susie_rss(z = z,
                               R = R,
                               prior_weights = prior,
                               prior_variance = V,
                               L = 1,
                               null_weight = null_weight,
                               coverage = coverage,
                               min_abs_corr = 0,
                               ...)

  susie_res_df <- anno_susie(susie_res,
                             gids = gids,
                             sids = sids,
                             g_type = g_type,
                             g_context = g_context,
                             g_group = g_group,
                             region_id = region_id,
                             z = z,
                             include_cs = include_cs)

  if (get_susie_alpha) {
    # extract alpha matrix from susie result
    susie_alpha_df <- get_susie_alpha_res(susie_res, susie_res_df, keep_genes_only = TRUE)
  } else {
    susie_alpha_df <- NULL
  }

  return(list("susie_res_df" = susie_res_df,
              "susie_alpha_df" = susie_alpha_df))
}

# finemap a single region with L = 1 without LD, used in EM
fast_finemap_single_region_L1_noLD <- function(region_data,
                                               region_id,
                                               pi_prior,
                                               V_prior,
                                               use_null_weight = TRUE,
                                               ...){
  # load region data
  regiondata <- extract_region_data(region_data, region_id)
  gids <- regiondata[["gid"]]
  sids <- regiondata[["sid"]]
  z <- regiondata[["z"]]
  gs_group <- regiondata[["gs_group"]]
  g_type <- regiondata[["g_type"]]
  g_context <- regiondata[["g_context"]]
  g_group <- regiondata[["g_group"]]
  rm(regiondata)

  if (length(z) < 2) {
    stop(paste(length(z), "variables in the region", region_id, "\n",
               "At least two variables in a region are needed to run susie"))
  }

  # update priors, prior variances and null_weight
  res <- set_region_susie_priors(pi_prior, V_prior, gs_group, L = 1, use_null_weight = use_null_weight)
  prior <- res$prior
  V <- res$V
  null_weight <- res$null_weight
  rm(res)

  # Use an identity matrix as LD, R does not matter for susie when L = 1
  R <- diag(length(z))

  # in susie, prior_variance is under standardized scale (if performed)
  susie_res <- ctwas_susie_rss(z = z,
                               R = R,
                               prior_weights = prior,
                               prior_variance = V,
                               L = 1,
                               null_weight = null_weight,
                               max_iter = 1,
                               warn_converge_fail = FALSE,
                               ...)

  # annotate susie result
  susie_res_df <- anno_susie(susie_res,
                             gids = gids,
                             sids = sids,
                             g_type = g_type,
                             g_context = g_context,
                             g_group = g_group,
                             region_id = region_id,
                             include_cs = FALSE)

  return(susie_res_df)
}
