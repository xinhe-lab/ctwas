#' run cTWAS finemapping for a single region
#'
#' @param region_data a list object indexing regions, variants and genes.
#'
#' @param region_id a character string of region ids to be finemapped
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

  if (verbose){
    loginfo("Finemapping region %s with L = %d", region_id, L)
  }

  regioninfo <- region_info[region_info$region_id == region_id, ]

  # get susie input data
  sid <- region_data[[region_id]][["sid"]]
  gid <- region_data[[region_id]][["gid"]]
  z <- region_data[[region_id]][["z"]]
  gs_group <- region_data[[region_id]][["gs_group"]]
  g_type <- region_data[[region_id]][["g_type"]]
  g_context <- region_data[[region_id]][["g_context"]]
  g_group <- region_data[[region_id]][["g_group"]]

  # set pi_prior and V_prior based on group_prior and group_prior_var
  groups <- unique(gs_group)
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
  if (length(region_id) > 1){
    region_id <- paste(region_id, collapse = "_")
  }
  R_sg_file <- file.path(cor_dir, paste0("region.", region_id, ".R_snp_gene.RDS"))
  R_g_file <- file.path(cor_dir, paste0("region.", region_id,  ".R_gene.RDS"))
  R_s_file <- file.path(cor_dir, paste0("region.", region_id, ".R_snp.RDS"))

  if (isTRUE(force_compute_cor)) {
    # force compute correlation matrix
    if (length(regioninfo$LD_matrix)==1){
      R_snp <- load_LD(regioninfo$LD_matrix)
    } else {
      R_snp <- lapply(regioninfo$LD_matrix, load_LD)
      R_snp <- suppressWarnings(as.matrix(Matrix::bdiag(R_snp)))
    }
    ld_snpinfo <- read_LD_SNP_files(regioninfo$SNP_info)
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
      if (length(regioninfo$LD_matrix)==1){
        R_snp <- load_LD(regioninfo$LD_matrix)
      } else {
        R_snp <- lapply(regioninfo$LD_matrix, load_LD)
        R_snp <- suppressWarnings(as.matrix(Matrix::bdiag(R_snp)))
      }
      ld_snpinfo <- read_LD_SNP_files(regioninfo$SNP_info)
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
    gene_info <- get_gene_info(weights)
    ld_snpinfo <- read_LD_SNP_files(regioninfo$SNP_info)
    susie_res_df <- anno_susie(susie_res,
                               gid = gid,
                               sid = sid,
                               region_id = region_id,
                               g_type = g_type,
                               g_context = g_context,
                               g_group = g_group,
                               geneinfo = gene_info,
                               snpinfo = ld_snpinfo,
                               include_cs_index = TRUE)
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
#' @importFrom logging loginfo
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

  loginfo('Finemapping %d regions ...', length(region_data))

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
