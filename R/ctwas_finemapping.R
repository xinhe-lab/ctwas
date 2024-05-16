#' run cTWAS finemapping for a single region
#'
#' @param region_data a list object with data for the regions
#'
#' @param region_id a character string of region id to be finemapped
#'
#' @param region_info a data frame of region definition and associated LD file names
#'
#' @param weights a list of weights for each gene
#'
#' @param snp_info a data frame, SNP info for LD reference,
#'  with columns "chrom", "id", "pos", and "region_id".
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
                           region_info,
                           weights,
                           snp_info,
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
                           annotate_susie_result = TRUE,
                           verbose = FALSE,
                           ...){
  # check input data
  if (!is.list(weights)){
    stop("'weights' should be a list.")
  }

  if (verbose){
    loginfo("Fine-mapping region %s with L = %d", region_id, L)
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
  } else if (L == 1 && !force_compute_cor) {
    # if L = 1, do not need LD
    R <- diag(length(z))
  } else {
    # if no precomputed correlation matrices, or L > 1, or force_compute_cor = TRUE,
    # compute correlation matrices
    if (verbose){
      loginfo("Compute correlation matrices for region %s", region_id)
    }
    # load LD matrix of the region
    regioninfo <- region_info[which(region_info$region_id == region_id), ]
    if(!is.character(regioninfo$LD_matrix) || is.null(regioninfo$LD_matrix)){
      stop("LD_matrix in region_info is required for computing correlation matrices")
    }
    LD_matrix_files <- unlist(strsplit(regioninfo$LD_matrix, split = ";"))
    stopifnot(all(file.exists(LD_matrix_files)))
    if (length(LD_matrix_files)==1) {
      R_snp <- load_LD(LD_matrix_files)
    } else {
      R_snp <- lapply(LD_matrix_files, load_LD)
      R_snp <- suppressWarnings(as.matrix(Matrix::bdiag(R_snp)))
    }
    # load SNP info of the region
    if (!is.character(regioninfo$SNP_info) || is.null(regioninfo$SNP_info)){
      snpinfo <- snp_info[which(snp_info$region_id == region_id), ]
    } else{
      snp_info_files <- unlist(strsplit(regioninfo$SNP_info, split = ";"))
      stopifnot(all(file.exists(snp_info_files)))
      snpinfo <- read_snp_info_files(snp_info_files)
    }

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

  # annotate susie result
  if (isTRUE(annotate_susie_result)) {
    geneinfo <- get_gene_info(weights[gid])
    snpinfo <- snp_info[which(snp_info$region_id == region_id), ]

    susie_res_df <- anno_susie(susie_res,
                               gid = gid,
                               sid = sid,
                               region_id = region_id,
                               z = z,
                               g_type = g_type,
                               g_context = g_context,
                               g_group = g_group,
                               geneinfo = geneinfo,
                               snpinfo = snpinfo,
                               include_cs_index = TRUE)

  } else {
    # skip annotating gene and SNP info, and cs_index
    susie_res_df <- anno_susie(susie_res,
                               gid = gid,
                               sid = sid,
                               region_id = region_id,
                               z = z,
                               g_type = g_type,
                               g_context = g_context,
                               g_group = g_group,
                               include_cs_index = FALSE)
  }

  return(susie_res_df)

}

#' run cTWAS finemapping for multiple regions
#'
#' @param region_data region_data to be finemapped
#'
#' @param region_info a data frame of region definition and associated LD file names
#'
#' @param weights a list of weights for each gene
#'
#' @param snp_info a data frame, SNP info for LD reference,
#'  with columns "chrom", "id", "pos", and "region_id".
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
                            region_info,
                            weights,
                            snp_info,
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
                            annotate_susie_result = TRUE,
                            verbose = FALSE,
                            logfile = NULL,
                            ...){

  if (!is.null(logfile)){
    addHandler(writeToFile, file= logfile, level='DEBUG')
  }
  # check input data
  if (!is.list(weights)){
    stop("'weights' should be a list.")
  }

  loginfo('Fine-mapping %d regions ...', length(region_data))

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
                   region_info = region_info,
                   weights = weights,
                   snp_info = snp_info,
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
                   annotate_susie_result = annotate_susie_result,
                   verbose = verbose,
                   ...)
  }, mc.cores = ncore)

  finemap_res <- do.call(rbind, finemap_region_res_list)

  return(finemap_res)
}
