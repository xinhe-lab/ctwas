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
#' @param max_iter Maximum number of IBSS iterations to perform.
#'
#' @param force_compute_cor TRUE/FALSE. If TRUE, force computing correlation (R) matrices
#'
#' @param save_cor TRUE/FALSE. If TRUE, save correlation (R) matrices
#'
#' @param cor_dir a string, the directory to store correlation (R) matrices
#'
#' @param annotate_susie_result TRUE/FALSE. If TRUE, add gene and SNP information and cs_index to
#' the data frame of finemapping results.
#'
#' @param verbose TRUE/FALSE. If TRUE, print detail messages
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
                           L = 5,
                           group_prior = NULL,
                           group_prior_var = NULL,
                           groups = NULL,
                           use_null_weight = TRUE,
                           coverage = 0.95,
                           min_abs_corr = 0.5,
                           max_iter = 100,
                           force_compute_cor = FALSE,
                           save_cor = FALSE,
                           cor_dir = getwd(),
                           annotate_susie_result = TRUE,
                           verbose = FALSE,
                           ...){

  # check input data
  if (!is.list(region_data)){
    stop("'region_data' should be a list.")
  }

  if (!is.data.frame(region_info)){
    stop("'region_info' should be a data frame.")
  }

  if (!is.list(weights)){
    stop("'weights' should be a list.")
  }

  if (!region_id %in% names(region_data)){
    stop("'region_data' does not contain 'region_id'.")
  }

  if (!region_id %in% region_info$region_id){
    stop("'region_info' does not contain 'region_id'.")
  }

  if (verbose){
    loginfo("Finemapping region %s with L = %d", region_id, L)
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

  # select region info for the region ids to finemap
  regioninfo <- region_info[region_info$region_id %in% region_id, ]

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

  # get LD matrix files and SNP info files
  LD_matrix_files <- unlist(strsplit(regioninfo$LD_matrix, split = ";"))
  stopifnot(all(file.exists(LD_matrix_files)))

  SNP_info_files <- unlist(strsplit(regioninfo$SNP_info, split = ";"))
  stopifnot(all(file.exists(SNP_info_files)))

  # load SNP info for the region
  ld_snpinfo <- read_LD_SNP_files(SNP_info_files)

  # compute correlation matrices
  R_sg_file <- file.path(cor_dir, paste0("region.", region_id, ".R_snp_gene.RDS"))
  R_g_file <- file.path(cor_dir, paste0("region.", region_id, ".R_gene.RDS"))
  R_s_file <- file.path(cor_dir, paste0("region.", region_id, ".R_snp.RDS"))

  if (isTRUE(force_compute_cor)) {
    # force compute correlation matrix
    if (length(LD_matrix_files) == 1){
      R_snp <- load_LD(LD_matrix_files)
    } else {
      R_snp <- lapply(LD_matrix_files, load_LD)
      R_snp <- suppressWarnings(as.matrix(Matrix::bdiag(R_snp)))
    }
    res <- compute_region_cor(sid, gid, R_snp, weights, ld_snpinfo)
    R_snp <- res$R_snp
    R_snp_gene <- res$R_snp_gene
    R_gene <- res$R_gene
    rm(res)
    if (isTRUE(save_cor)) {
      if (!dir.exists(cor_dir))
        dir.create(cor_dir, recursive = TRUE)
      saveRDS(R_snp_gene, file=R_sg_file)
      saveRDS(R_gene, file=R_g_file)
      saveRDS(R_snp, file=R_s_file)
    }
    # gene first then SNPs
    R <- rbind(cbind(R_gene, t(R_snp_gene)),
               cbind(R_snp_gene, R_snp))
    rm(R_gene, R_snp_gene, R_snp)
  } else {
    if (all(file.exists(c(R_sg_file, R_g_file, R_s_file)))) {
      # load precomputed correlation matrices
      R_snp_gene <- load_LD(R_sg_file)
      R_gene <- load_LD(R_g_file)
      R_snp <- load_LD(R_s_file)
      # gene first then SNPs
      R <- rbind(cbind(R_gene, t(R_snp_gene)),
                 cbind(R_snp_gene, R_snp))
      rm(R_gene, R_snp_gene, R_snp)
    } else if (L == 1){
      # R does not matter for susie when L = 1
      # loginfo("L = 1, skip computing correlation matrices")
      R <- diag(length(z))
    } else {
      # compute correlation matrix if L > 1
      if (length(LD_matrix_files)==1){
        R_snp <- load_LD(LD_matrix_files)
      } else {
        R_snp <- lapply(LD_matrix_files, load_LD)
        R_snp <- suppressWarnings(as.matrix(Matrix::bdiag(R_snp)))
      }
      res <- compute_region_cor(sid, gid, R_snp, weights, ld_snpinfo)
      R_snp <- res$R_snp
      R_snp_gene <- res$R_snp_gene
      R_gene <- res$R_gene
      rm(res)
      if (isTRUE(save_cor)) {
        if (!dir.exists(cor_dir))
          dir.create(cor_dir, recursive = TRUE)
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
                               max_iter = max_iter,
                               ...)

  # annotate susie result
  if (isTRUE(annotate_susie_result)) {

    geneinfo <- get_gene_info(weights[gid])

    susie_res_df <- anno_susie(susie_res,
                               gid = gid,
                               sid = sid,
                               region_id = region_id,
                               g_type = g_type,
                               g_context = g_context,
                               g_group = g_group,
                               geneinfo = geneinfo,
                               snpinfo = ld_snpinfo,
                               include_cs_index = TRUE)
    # add z-scores to finemapping result
    susie_res_df$z <- z

  } else {
    # skip annotating gene and SNP info, and cs_index
    susie_res_df <- anno_susie(susie_res,
                               gid = gid,
                               sid = sid,
                               region_id = region_id,
                               g_type = g_type,
                               g_context = g_context,
                               g_group = g_group,
                               geneinfo = NULL,
                               snpinfo = NULL,
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
#' @param max_iter Maximum number of IBSS iterations to perform.
#'
#' @param ncore The number of cores used to parallelize computation over regions
#'
#' @param force_compute_cor TRUE/FALSE. If TRUE, force computing correlation (R) matrices
#'
#' @param save_cor TRUE/FALSE. If TRUE, save correlation (R) matrices
#'
#' @param cor_dir a string, the directory to store correlation (R) matrices
#'
#' @param verbose TRUE/FALSE. If TRUE, print detail messages
#'
#' @param logfile the log file, if NULL will print log info on screen
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
                            L = 5,
                            group_prior = NULL,
                            group_prior_var = NULL,
                            groups = NULL,
                            use_null_weight = TRUE,
                            coverage = 0.95,
                            min_abs_corr = 0.5,
                            max_iter = 100,
                            ncore = 1,
                            force_compute_cor = FALSE,
                            save_cor = FALSE,
                            cor_dir = getwd(),
                            annotate_susie_result = TRUE,
                            verbose = FALSE,
                            logfile = NULL,
                            ...){

  if (!is.null(logfile)){
    addHandler(writeToFile, file= logfile, level='DEBUG')
  }

  # check input data
  if (!is.list(region_data)){
    stop("'region_data' should be a list.")
  }

  if (!is.data.frame(region_info)){
    stop("'region_info' should be a data frame.")
  }

  if (!is.list(weights)){
    stop("'weights' should be a list.")
  }

  if (!all(names(region_data) %in% region_info$region_id)){
    stop("Some 'region_id' in 'region_data' were missed in 'region_info'.")
  }

  loginfo('Finemapping %d regions ...', length(region_data))

  if (is.null(groups)){
    if(!is.null(group_prior)){
      groups <- names(group_prior)
    }else{
      groups <- unique(unlist(lapply(region_data, "[[", "gs_group")))
    }
  }

  cl <- parallel::makeCluster(ncore, outfile = "")
  doParallel::registerDoParallel(cl)

  corelist <- region2core(region_data, ncore)

  finemap_res <- foreach (core = 1:length(corelist), .combine = "rbind", .packages = "ctwas") %dopar% {

    finemap_res.core.list <- list()
    # run finemapping for each region
    region_ids.core <- corelist[[core]]
    for (region_id in region_ids.core) {
      finemap_res.core.list[[region_id]] <- finemap_region(region_data = region_data,
                                                           region_id = region_id,
                                                           region_info = region_info,
                                                           weights = weights,
                                                           L = L,
                                                           group_prior = group_prior,
                                                           group_prior_var = group_prior_var,
                                                           groups = groups,
                                                           use_null_weight = use_null_weight,
                                                           coverage = coverage,
                                                           min_abs_corr = min_abs_corr,
                                                           max_iter = max_iter,
                                                           force_compute_cor = force_compute_cor,
                                                           save_cor = save_cor,
                                                           cor_dir = cor_dir,
                                                           annotate_susie_result = annotate_susie_result,
                                                           verbose = verbose,
                                                           ...)
    }
    finemap_res.core <- do.call(rbind, finemap_res.core.list)
    finemap_res.core
  }
  parallel::stopCluster(cl)
  rownames(finemap_res) <- NULL

  return(finemap_res)
}
