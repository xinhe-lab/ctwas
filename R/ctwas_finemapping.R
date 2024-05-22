#' run cTWAS finemapping for a single region
#'
#' @param region_data a list object with data for the regions
#'
#' @param region_id a character string of region id to be finemapped
#'
#' @param weights a list of preprocessed weights
#'
#' @param use_LD TRUE/FALSE. If TRUE, use LD for finemapping. Otherwise, use "no-LD" version.
#'
#' @param LD_info a list of paths to LD matrices for each of the regions. Required when \code{use_LD = TRUE}.
#'
#' @param snp_info a list of SNP info data frames for LD reference. Required when \code{use_LD = TRUE}.
#'
#' @param L the number of effects for susie during the fine mapping steps
#'
#' @param group_prior a vector of two prior inclusion probabilities for SNPs and genes.
#'
#' @param group_prior_var a vector of two prior variances for SNPs and gene effects.
#'
#' @param groups a vector of group names for group_prior and group_prior_var.
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
#' @param annotate_susie_result TRUE/FALSE. If TRUE, add gene and SNP information and cs_index to
#' the data frame of finemapping results.
#'
#' @param verbose TRUE/FALSE. If TRUE, print detail messages
#'
#' @param ... Additional arguments of \code{susie_rss}.
#'
#' @importFrom logging loginfo
#'
#' @return a data frame of finemapping results.
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
                           groups = NULL,
                           use_null_weight = TRUE,
                           coverage = 0.95,
                           min_abs_corr = 0.5,
                           force_compute_cor = FALSE,
                           save_cor = FALSE,
                           cor_dir = NULL,
                           include_cs_index = TRUE,
                           verbose = FALSE,
                           ...){

  if (verbose){
    loginfo("Fine-mapping region %s with L = %d", region_id, L)
  }

  # check weights
  if (!is.null(weights)){
    if (!is.list(weights)){
      stop("'weights' should be a list.")
    }
  }

  # load input data for the region
  regiondata <- region_data[[region_id]]
  sid <- regiondata[["sid"]]
  gid <- regiondata[["gid"]]
  z <- regiondata[["z"]]
  gs_group <- regiondata[["gs_group"]]
  g_type <- regiondata[["g_type"]]
  g_context <- regiondata[["g_context"]]
  g_group <- regiondata[["g_group"]]

  # set pi_prior and V_prior based on group_prior and group_prior_var
  if (is.null(groups)){
    if(!is.null(group_prior)){
      groups <- names(group_prior)
    }else{
      stop("'groups' is required when group_prior is null")
    }
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
    if (L != 1){
      warning("L has to be 1 for no-LD version. Set L = 1")
      L <- 1
    }
    R <- diag(length(z))
    include_cs_index <- FALSE
  } else {
    if (verbose){
      loginfo("Fine-mapping using LD version ...")
    }
    # compute correlation matrices
    if (!is.null(cor_dir)) {
      if (!dir.exists(cor_dir))
        dir.create(cor_dir, recursive = TRUE)
      R_sg_file <- file.path(cor_dir, paste0("region.", region_id, ".R_snp_gene.RDS"))
      R_g_file <- file.path(cor_dir, paste0("region.", region_id, ".R_gene.RDS"))
      R_s_file <- file.path(cor_dir, paste0("region.", region_id, ".R_snp.RDS"))
    }
    cor_files_exist <- isTRUE(!is.null(cor_dir) && file.exists(R_sg_file) && file.exists(R_g_file) && file.exists(R_s_file))
    if (cor_files_exist && !force_compute_cor) {
      if (verbose){
        loginfo("Load correlation matrices for region %s", region_id)
      }
      # load precomputed correlation matrices
      R_snp_gene <- load_LD(R_sg_file)
      R_gene <- load_LD(R_g_file)
      R_snp <- load_LD(R_s_file)
      # gene first then SNPs
      R <- rbind(cbind(R_gene, t(R_snp_gene)),
                 cbind(R_snp_gene, R_snp))
      rm(R_gene, R_snp_gene, R_snp)
    } else {
      # if no precomputed correlation matrices, or force_compute_cor = TRUE,
      # compute correlation matrices
      if (verbose){
        loginfo("Compute correlation matrices for region %s", region_id)
      }
      # load LD matrix of the region
      if (is.null(LD_info) || is.null(snp_info)) {
        stop("LD_info and snp_info are required for computing correlation matrices")
      }
      LD_matrix_files <- LD_info[[region_id]]$LD_matrix
      stopifnot(all(file.exists(LD_matrix_files)))
      if (length(LD_matrix_files)==1) {
        R_snp <- load_LD(LD_matrix_files)
      } else {
        R_snp <- lapply(LD_matrix_files, load_LD)
        R_snp <- suppressWarnings(as.matrix(Matrix::bdiag(R_snp)))
      }
      # load SNP info of the region
      snpinfo <- snp_info[[region_id]]

      # Compute correlation matrices
      res <- compute_region_cor(sid, gid, R_snp, snpinfo$id, weights)
      R_snp <- res$R_snp
      R_snp_gene <- res$R_snp_gene
      R_gene <- res$R_gene
      rm(res)
      # save correlation matrices
      if (isTRUE(save_cor && !is.null(cor_dir))) {
        saveRDS(R_snp_gene, file=R_sg_file)
        saveRDS(R_gene, file=R_g_file)
        saveRDS(R_snp, file=R_s_file)
      }
      # gene first then SNPs
      R <- rbind(cbind(R_gene, t(R_snp_gene)),
                 cbind(R_snp_gene, R_snp))
      rm(R_gene, R_snp_gene, R_snp)
    }
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
                             gid = gid,
                             sid = sid,
                             region_id = region_id,
                             z = z,
                             g_type = g_type,
                             g_context = g_context,
                             g_group = g_group,
                             include_cs_index = include_cs_index)

  return(susie_res_df)

}

#' run cTWAS finemapping for multiple regions
#'
#' @param region_data region_data to be finemapped
#'
#' @param weights a list of weights for each gene
#'
#' @param use_LD TRUE/FALSE. If TRUE, use LD for finemapping. Otherwise, use "no-LD" version.
#'
#' @param LD_info a list of paths to LD matrices for each of the regions. Required when \code{use_LD = TRUE}.
#'
#' @param snp_info a list of SNP info data frames for LD reference. Required when \code{use_LD = TRUE}.
#'
#' @param L the number of effects for susie during the fine mapping steps
#'
#' @param group_prior a vector of two prior inclusion probabilities for SNPs and genes.
#'
#' @param group_prior_var a vector of two prior variances for SNPs and gene effects.
#'
#' @param groups a vector of group names for group_prior and group_prior_var.
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
#' @param verbose TRUE/FALSE. If TRUE, print detail messages
#'
#' @param logfile the log file, if NULL will print log info on screen
#'
#' @param ... Additional arguments of \code{susie_rss}.
#'
#' @importFrom logging addHandler loginfo writeToFile
#'
#' @return finemapping results.
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
                            groups = NULL,
                            use_null_weight = TRUE,
                            coverage = 0.95,
                            min_abs_corr = 0.5,
                            ncore = 1,
                            force_compute_cor = FALSE,
                            save_cor = FALSE,
                            cor_dir = NULL,
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
    if (!is.list(weights)){
      stop("'weights' should be a list.")
    }
  }
  if (!use_LD) {
    if (L != 1){
      warning("L has to be 1 for no-LD version. Set L = 1")
      L <- 1
    }
  }

  if (is.null(groups)){
    if(!is.null(group_prior)){
      groups <- names(group_prior)
    }else{
      groups <- unique(unlist(lapply(region_data, "[[", "gs_group")))
    }
  }

  region_ids <- names(region_data)

  finemap_region_res_list <- parallel::mclapply(region_ids, function(region_id){
    finemap_region(region_data = region_data,
                   region_id = region_id,
                   use_LD = use_LD,
                   LD_info = LD_info,
                   snp_info = snp_info,
                   weights = weights,
                   L = L,
                   group_prior = group_prior,
                   group_prior_var = group_prior_var,
                   groups = groups,
                   use_null_weight = use_null_weight,
                   coverage = coverage,
                   min_abs_corr = min_abs_corr,
                   force_compute_cor = force_compute_cor,
                   save_cor = save_cor,
                   cor_dir = cor_dir,
                   include_cs_index = include_cs_index,
                   verbose = verbose,
                   ...)
  }, mc.cores = ncore)

  finemap_res <- do.call(rbind, finemap_region_res_list)

  return(finemap_res)
}
