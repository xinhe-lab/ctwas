#' @title Runs cTWAS finemapping for a single region
#'
#' @param region_data a list object with data for the regions
#'
#' @param region_id a character string of region id to be finemapped
#'
#' @param use_LD TRUE/FALSE. If TRUE, use LD for finemapping. Otherwise, use "no-LD" version.
#'
#' @param LD_info a list of paths to LD matrices for each of the regions. Required when \code{use_LD = TRUE}.
#'
#' @param snp_info a list of SNP info data frames for LD reference. Required when \code{use_LD = TRUE}.
#'
#' @param weights a list of preprocessed weights
#'
#' @param L the number of effects for susie during the fine mapping
#'
#' @param group_prior a vector of two prior inclusion probabilities for SNPs and genes.
#'
#' @param group_prior_var a vector of two prior variances for SNPs and gene effects.
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
#' @param force_compute_cor TRUE/FALSE. If TRUE, force computing correlation (R) matrices
#'
#' @param save_cor TRUE/FALSE. If TRUE, save correlation (R) matrices to \code{cor_dir}
#'
#' @param cor_dir a string, the directory to store correlation (R) matrices
#'
#' @param LD_format file format for LD matrix. If "custom", use a user defined
#' \code{LD_loader_fun()} function to load LD matrix.
#'
#' @param LD_loader_fun a user defined function to load LD matrix when \code{LD_format = "custom"}.
#'
#' @param include_cs_index TRUE/FALSE. If TRUE, add cs_index to finemapping results.
#'
#' @param verbose TRUE/FALSE. If TRUE, print detail messages
#'
#' @param ... Additional arguments of \code{susie_rss}.
#'
#' @return a data frame of finemapping results.
#'
#' @importFrom logging loginfo
#' @importFrom Matrix bdiag
#'
#' @export
#'
finemap_region <- function(region_data,
                           region_id,
                           use_LD = TRUE,
                           LD_info = NULL,
                           snp_info = NULL,
                           weights = NULL,
                           L = 5,
                           group_prior = NULL,
                           group_prior_var = NULL,
                           use_null_weight = TRUE,
                           coverage = 0.95,
                           min_abs_corr = 0.5,
                           force_compute_cor = FALSE,
                           save_cor = FALSE,
                           cor_dir = NULL,
                           LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
                           LD_loader_fun,
                           include_cs_index = TRUE,
                           verbose = FALSE,
                           ...){

  if (verbose){
    loginfo("Fine-mapping region %s with L = %d", region_id, L)
  }

  # check weights
  if (!is.null(weights)){
    if (!inherits(weights,"list")){
      stop("'weights' should be a list.")
    }
  }

  LD_format <- match.arg(LD_format)

  # load input data for the region
  regiondata <- region_data[[region_id]]
  sids <- regiondata[["sid"]]
  gids <- regiondata[["gid"]]
  z <- regiondata[["z"]]
  gs_group <- regiondata[["gs_group"]]
  g_type <- regiondata[["g_type"]]
  g_context <- regiondata[["g_context"]]
  g_group <- regiondata[["g_group"]]

  # set pi_prior and V_prior based on group_prior and group_prior_var
  if(!is.null(group_prior)){
    groups <- names(group_prior)
  }else{
    groups <- unique(unlist(lapply(region_data, "[[", "gs_group")))
  }
  res <- initiate_group_priors(group_prior, group_prior_var, groups)
  pi_prior <- res$pi_prior
  V_prior <- res$V_prior
  rm(res)

  # set prior and prior variance values for the region
  res <- set_region_susie_priors(pi_prior, V_prior, gs_group, L = L, use_null_weight = use_null_weight)
  prior <- res$prior
  V <- res$V
  null_weight <- res$null_weight
  rm(res)

  if (!use_LD) {
    if (verbose){
      loginfo("Fine-mapping using no-LD version ...")
    }
    # use an identity matrix as R in no-LD version
    R <- diag(length(z))
    # do not include cs_index in no-LD version
    include_cs_index <- FALSE
  } else {
    if (verbose){
      loginfo("Fine-mapping using LD version ...")
    }
    cor_res <- get_region_cor(region_id,
                              region_data = region_data,
                              LD_info = LD_info,
                              snp_info = snp_info,
                              weights = weights,
                              force_compute_cor = force_compute_cor,
                              save_cor = save_cor,
                              cor_dir = cor_dir,
                              LD_format = LD_format,
                              LD_loader_fun = LD_loader_fun,
                              verbose = verbose)
    # gene first then SNPs
    R <- rbind(cbind(cor_res$R_gene, t(cor_res$R_snp_gene)),
               cbind(cor_res$R_snp_gene, cor_res$R_snp))
    rm(cor_res)
  }

  if (anyNA(R))
    stop("R matrix contains missing values!")

  if (length(z) != nrow(R))
    stop("R matrix dimension does not match with z!")

  # run susie for this region
  # in susie, prior_variance is under standardized scale (if performed)
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
                             region_id = region_id,
                             z = z,
                             g_type = g_type,
                             g_context = g_context,
                             g_group = g_group,
                             include_cs_index = include_cs_index)

  return(susie_res_df)

}

#' @title Runs cTWAS finemapping for multiple regions
#'
#' @param region_data region_data to be finemapped
#'
#' @param use_LD TRUE/FALSE. If TRUE, use LD for finemapping. Otherwise, use "no-LD" version.
#'
#' @param LD_info a list of paths to LD matrices for each of the regions. Required when \code{use_LD = TRUE}.
#'
#' @param snp_info a list of SNP info data frames for LD reference. Required when \code{use_LD = TRUE}.
#'
#' @param weights a list of weights for each gene
#'
#' @param L the number of effects or a vector of number of effects for each region.
#'
#' @param group_prior a vector of two prior inclusion probabilities for SNPs and genes.
#'
#' @param group_prior_var a vector of two prior variances for SNPs and gene effects.
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
#' @param force_compute_cor TRUE/FALSE. If TRUE, force computing correlation (R) matrices
#'
#' @param save_cor TRUE/FALSE. If TRUE, save correlation (R) matrices to \code{cor_dir}
#'
#' @param cor_dir a string, the directory to store correlation (R) matrices
#'
#' @param LD_format file format for LD matrix. If "custom", use a user defined
#' \code{LD_loader_fun()} function to load LD matrix.
#'
#' @param LD_loader_fun a user defined function to load LD matrix when \code{LD_format = "custom"}.
#'
#' @param include_cs_index TRUE/FALSE. If TRUE, add cs_index to finemapping results.
#'
#' @param verbose TRUE/FALSE. If TRUE, print detail messages
#'
#' @param logfile the log file, if NULL will print log info on screen
#'
#' @param ... Additional arguments of \code{susie_rss}.
#'
#' @return a data frame of cTWAS finemapping results.
#'
#' @importFrom logging addHandler loginfo writeToFile
#' @importFrom parallel mclapply
#'
#' @export
#'
finemap_regions <- function(region_data,
                            use_LD = TRUE,
                            LD_info = NULL,
                            snp_info = NULL,
                            weights = NULL,
                            L = 5,
                            group_prior = NULL,
                            group_prior_var = NULL,
                            use_null_weight = TRUE,
                            coverage = 0.95,
                            min_abs_corr = 0.5,
                            ncore = 1,
                            force_compute_cor = FALSE,
                            save_cor = FALSE,
                            cor_dir = NULL,
                            LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
                            LD_loader_fun,
                            include_cs_index = TRUE,
                            verbose = FALSE,
                            logfile = NULL,
                            ...){

  if (!is.null(logfile)){
    addHandler(writeToFile, file= logfile, level='DEBUG')
  }

  loginfo('Fine-mapping %d regions ...', length(region_data))

  # check weights
  if (!is.null(weights)){
    if (!inherits(weights,"list")){
      stop("'weights' should be a list.")
    }
  }

  LD_format <- match.arg(LD_format)

  if (!use_LD) {
    if (L != 1){
      loginfo("L has to be 1 for no-LD version. Set L = 1")
      L <- 1
    }
  }

  region_ids <- names(region_data)

  finemap_region_res_list <- mclapply(region_ids, function(region_id){

    if (length(L) == 1) {
      region_L <- L
    } else if (length(L) > 1 && length(L) == length(region_data)) {
      region_L <- L[region_id]
    } else{
      stop("L needs to an integer or a vector of the same length as region_data")
    }

    finemap_region(region_data = region_data,
                   region_id = region_id,
                   use_LD = use_LD,
                   LD_info = LD_info,
                   snp_info = snp_info,
                   weights = weights,
                   L = region_L,
                   group_prior = group_prior,
                   group_prior_var = group_prior_var,
                   use_null_weight = use_null_weight,
                   coverage = coverage,
                   min_abs_corr = min_abs_corr,
                   force_compute_cor = force_compute_cor,
                   save_cor = save_cor,
                   cor_dir = cor_dir,
                   LD_format = LD_format,
                   LD_loader_fun = LD_loader_fun,
                   include_cs_index = include_cs_index,
                   verbose = verbose,
                   ...)

  }, mc.cores = ncore)
  check_mc_res(finemap_region_res_list)

  finemap_res <- do.call(rbind, finemap_region_res_list)
  rownames(finemap_res) <- NULL

  return(finemap_res)
}
