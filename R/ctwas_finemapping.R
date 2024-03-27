#' run cTWAS finemapping for a single region
#'
#' @param z_snp A data frame with four columns: "id", "A1", "A2", "z".
#' giving the z scores for snps. "A1" is effect allele. "A2" is the other allele.
#'
#' @param z_gene A data frame with two columns: "id", "z". giving the z scores for genes.
#' Optionally, a "type" column can also be supplied; this is for using multiple sets of weights
#'
#' @param sid SNP IDs in the region
#'
#' @param gid gene IDs in the region
#'
#' @param region_tag a character string of region tags to be finemapped
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
#' @param verbose TRUE/FALSE. If TRUE, print detail messages
#'
#' @importFrom logging loginfo
#'
#' @return finemapping results.
#'
#' @export
#'
finemap_region <- function(z_snp,
                           z_gene,
                           sid,
                           gid,
                           region_tag,
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
                           verbose = FALSE,
                           ...){

  if (verbose){
    loginfo("Finemapping region %s with L = %d", region_tag, L)
  }

  # prepare susie input data
  zdf <- combine_z(z_snp, z_gene)

  # set pi_prior and V_prior based on group_prior and group_prior_var
  types <- unique(zdf$type)
  if (is.null(group_prior)){
    group_prior <- structure(as.numeric(rep(NA,length(types))), names=types)
  }
  if (is.null(group_prior_var)){
    group_prior_var <- structure(as.numeric(rep(NA,length(types))), names=types)
  }
  pi_prior <- list()
  V_prior <- list()
  for (type in types){
    pi_prior[[type]] <- unname(group_prior[type])
    V_prior[[type]] <- unname(group_prior_var[type])
  }
  pi_prior <- unlist(pi_prior)
  V_prior <- unlist(V_prior)

  susie_input <- assemble_susie_input(sid, gid, zdf, pi_prior, V_prior,
                                      L = L, use_null_weight = TRUE)
  sid <- susie_input$sid
  gid <- susie_input$gid
  z <- susie_input$z
  g_type <- susie_input$g_type
  g_QTLtype <- susie_input$g_QTLtype
  gs_type <- susie_input$gs_type
  prior <- susie_input$prior
  V <- susie_input$V
  null_weight <- susie_input$null_weight

  # get gene info from weights
  gene_info <- get_gene_info(weights)

  # LD and SNP information for this region
  regioninfo <- region_info[region_info$region_tag %in% region_tag, ]
  LD_matrix_file <- regioninfo$LD_matrix
  ld_snpinfo <- read_LD_SNP_files(regioninfo$SNP_info)

  # compute correlation matrices
  if (length(region_tag) > 1){
    region_tag <- paste(region_tag, collapse = "_")
  }
  R_sg_file <- file.path(cor_dir, paste0("region.", region_tag, ".R_snp_gene.RDS"))
  R_g_file <- file.path(cor_dir, paste0("region.", region_tag,  ".R_gene.RDS"))
  R_s_file <- file.path(cor_dir, paste0("region.", region_tag, ".R_snp.RDS"))

  if (isTRUE(force_compute_cor)) {
    # force compute correlation matrix
    if (verbose){
      loginfo("Compute correlation matrices ...")
    }
    if (length(LD_matrix_file)==1){
      R_snp <- load_LD(LD_matrix_file)
    } else {
      R_snp <- lapply(LD_matrix_file, load_LD)
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
      if (length(LD_matrix_file)==1){
        R_snp <- load_LD(LD_matrix_file)
      } else {
        R_snp <- lapply(LD_matrix_file, load_LD)
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
  susie_res_df <- anno_susie(susie_res,
                             gid = gid,
                             sid = sid,
                             g_type = g_type,
                             g_QTLtype = g_QTLtype,
                             region_tag = region_tag,
                             geneinfo = gene_info,
                             snpinfo = ld_snpinfo)

  return(susie_res_df)

}

#' run cTWAS finemapping for multiple regions
#'
#' @param z_snp A data frame with four columns: "id", "A1", "A2", "z".
#' giving the z scores for snps. "A1" is effect allele. "A2" is the other allele.
#'
#' @param z_gene A data frame with two columns: "id", "z". giving the z scores for genes.
#' Optionally, a "type" column can also be supplied; this is for using multiple sets of weights
#'
#' @param regionlist regionlist to be finemapped
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
finemap_regions <- function(z_snp,
                            z_gene,
                            regionlist,
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
                            verbose = FALSE,
                            logfile = NULL,
                            ...){

  if (!is.null(logfile)){
    addHandler(writeToFile, file= logfile, level='DEBUG')
  }

  loginfo('Finemapping %d regions ...', length(regionlist))

  cl <- parallel::makeCluster(ncore, outfile = "")
  doParallel::registerDoParallel(cl)

  corelist <- region2core(regionlist, ncore)

  finemap_res <- foreach (core = 1:length(corelist), .combine = "rbind", .packages = "ctwas") %dopar% {

    finemap_res.core.list <- list()
    # run finemapping for each region
    region_tags.core <- corelist[[core]]
    for (region_tag in region_tags.core) {
      sid <- regionlist[[region_tag]][["sid"]]
      gid <- regionlist[[region_tag]][["gid"]]
      finemap_res.core.list[[region_tag]] <- finemap_region(z_snp,
                                                            z_gene,
                                                            sid,
                                                            gid,
                                                            region_tag,
                                                            region_info,
                                                            weights,
                                                            L = L,
                                                            group_prior = group_prior,
                                                            group_prior_var = group_prior_var,
                                                            use_null_weight = use_null_weight,
                                                            coverage = coverage,
                                                            min_abs_corr = min_abs_corr,
                                                            max_iter = max_iter,
                                                            save_cor = save_cor,
                                                            cor_dir = cor_dir,
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
