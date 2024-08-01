#' @title Runs cTWAS fine-mapping for regions
#'
#' @param region_data region_data to be finemapped
#'
#' @param use_LD If TRUE, use LD for finemapping. Otherwise, use "no-LD" version.
#'
#' @param LD_map a data frame with filenames of LD matrices for the regions. Required when \code{use_LD = TRUE}.
#'
#' @param snp_map a list of SNP-to-region map for the reference. Required when \code{use_LD = TRUE}.
#'
#' @param weights a list of preprocessed weights.
#'
#' @param L the number of effects or a vector of number of effects for each region.
#'
#' @param group_prior a vector of two prior inclusion probabilities for SNPs and genes.
#'
#' @param group_prior_var a vector of two prior variances for SNPs and gene effects.
#'
#' @param use_null_weight If TRUE, allow for a probability of no effect in susie
#'
#' @param coverage A number between 0 and 1 specifying the \dQuote{coverage} of the estimated confidence sets
#'
#' @param min_abs_corr Minimum absolute correlation allowed in a
#'   credible set. The default, 0.5, corresponds to a squared
#'   correlation of 0.25, which is a commonly used threshold for
#'   genotype data in genetic studies.
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
#' @param include_cs_index If TRUE, add cs_index to finemapping results.
#'
#' @param ncore The number of cores used to parallelize computation over regions
#'
#' @param verbose If TRUE, print detail messages
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
                            LD_map,
                            snp_map,
                            weights,
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
                            ncore = 1,
                            verbose = FALSE,
                            logfile = NULL,
                            ...){

  if (!is.null(logfile)){
    addHandler(writeToFile, file=logfile, level='DEBUG')
  }

  if (use_LD){
    loginfo("Fine-mapping %d regions using LD version ...", length(region_data))
  } else {
    loginfo("Fine-mapping %d regions using no-LD version ...", length(region_data))
  }

  # check inputs
  LD_format <- match.arg(LD_format)

  if (!inherits(region_data,"list"))
    stop("'region_data' should be a list.")

  if (use_LD) {
    if (missing(LD_map) || missing(snp_map) || missing(weights))
      stop("'LD_map', 'snp_map' and 'weights' are required when use_LD = TRUE")

    if (!inherits(LD_map,"data.frame"))
      stop("'LD_map' should be a data frame")

    if (!inherits(snp_map,"list"))
      stop("'snp_map' should be a list.")

    if (!inherits(weights,"list"))
      stop("'weights' should be a list.")

    if (any(sapply(weights, is.null)))
      stop("'weights' contain NULL, remove empty weights!")
  }

  # set L = 1 in no-LD version
  if (!use_LD){
    L <- 1
  }

  if (length(L) == 1) {
    L <- rep(L, length(region_data))
    names(L) <- names(region_data)
  } else if (length(L) == length(region_data)) {
    if (!all.equal(names(L), names(region_data)))
      stop("the names of L do not match with region_data!")
  } else {
    stop("L needs to an integer or a vector of the same length as region_data!")
  }

  if (verbose) {
    if (is.null(group_prior)) {
      loginfo("Use uniform prior")
    } else {
      loginfo("group_prior {%s}: {%s}", names(group_prior), format(group_prior, digits = 4))
      loginfo("group_prior_var {%s}: {%s}", names(group_prior_var), format(group_prior_var, digits = 4))
    }
  }

  region_ids <- names(region_data)

  finemap_res_list <- mclapply_check(region_ids, function(region_id){

    finemap_region(region_data = region_data,
                   region_id = region_id,
                   use_LD = use_LD,
                   LD_map = LD_map,
                   snp_map = snp_map,
                   weights = weights,
                   L = L[region_id],
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

  }, mc.cores = ncore, stop_if_missing = TRUE)

  finemap_res <- do.call(rbind, finemap_res_list)
  rownames(finemap_res) <- NULL

  return(finemap_res)
}

#' @title Runs cTWAS finemapping for a single region
#'
#' @param region_data a list object with data for the regions
#'
#' @param region_id a character string of region id to be finemapped
#'
#' @param use_LD If TRUE, use LD for finemapping. Otherwise, use "no-LD" version.
#'
#' @param LD_map a data frame with filenames of LD matrices for the regions. Required when \code{use_LD = TRUE}.
#'
#' @param snp_map a list of SNP-to-region map for the reference. Required when \code{use_LD = TRUE}.
#'
#' @param weights a list of preprocessed weights. Required when \code{use_LD = TRUE}.
#'
#' @param L the number of effects for susie during the fine mapping
#'
#' @param group_prior a vector of two prior inclusion probabilities for SNPs and genes.
#'
#' @param group_prior_var a vector of two prior variances for SNPs and gene effects.
#'
#' @param use_null_weight If TRUE, allow for a probability of no effect in susie
#'
#' @param coverage A number between 0 and 1 specifying the \dQuote{coverage} of the estimated confidence sets
#'
#' @param min_abs_corr Minimum absolute correlation allowed in a
#'   credible set. The default, 0.5, corresponds to a squared
#'   correlation of 0.25, which is a commonly used threshold for
#'   genotype data in genetic studies.
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
#' @param include_cs_index If TRUE, add cs_index to finemapping results.
#'
#' @param verbose If TRUE, print detail messages
#'
#' @param ... Additional arguments of \code{susie_rss}.
#'
#' @return a data frame of finemapping results.
#'
#' @importFrom logging loginfo
#' @importFrom Matrix bdiag
#'
#' @keywords internal
#'
finemap_region <- function(region_data,
                           region_id,
                           use_LD = TRUE,
                           LD_map,
                           snp_map,
                           weights,
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

  # check inputs
  LD_format <- match.arg(LD_format)

  if (!inherits(region_data,"list")){
    stop("'region_data' should be a list.")
  }

  if (use_LD) {
    if (missing(LD_map) || missing(snp_map) || missing(weights))
      stop("'LD_map', 'snp_map' and 'weights' are required when use_LD = TRUE")

    if (!inherits(LD_map,"data.frame"))
      stop("'LD_map' should be a data frame")

    if (!inherits(snp_map,"list"))
      stop("'snp_map' should be a list.")

    if (!inherits(weights,"list"))
      stop("'weights' should be a list.")

    if (any(sapply(weights, is.null)))
      stop("'weights' contain NULL, remove empty weights!")
  }

  # Set L = 1 in no-LD version
  if (!use_LD){
    L <- 1
  }

  # load input data for the region
  regiondata <- extract_region_data(region_data, region_id)
  gids <- regiondata$gid
  sids <- regiondata$sid
  z <- regiondata$z
  gs_group <- regiondata$gs_group
  groups <- regiondata$groups
  rm(regiondata)

  if (verbose){
    loginfo("%d genes, %d SNPs in the region", length(gids), length(sids))
  }

  # set pi_prior and V_prior based on group_prior and group_prior_var
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

  if (!use_LD) {
    # use an identity matrix as R in no-LD version
    R <- diag(length(z))
    # do not include cs_index in no-LD version
    include_cs_index <- FALSE
  } else {
    cor_res <- get_region_cor(region_id,
                              region_data = region_data,
                              LD_map = LD_map,
                              snp_map = snp_map,
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
                             include_cs_index = include_cs_index)

  return(susie_res_df)

}

