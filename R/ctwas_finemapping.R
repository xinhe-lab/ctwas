#' run cTWAS finemapping for a single region
#'
#' @param z_snp A data frame with four columns: "id", "A1", "A2", "z".
#' giving the z scores for snps. "A1" is effect allele. "A2" is the other allele.
#'
#' @param z_gene A data frame with two columns: "id", "z". giving the z scores for genes.
#' Optionally, a "type" column can also be supplied; this is for using multiple sets of weights
#'
#' @param gene_info a data frame of gene information obtained from \code{compute_gene_z}
#'
#' @param regionlist a list object indexing regions, variants and genes.
#'
#' @param region_tag a character string of region tag to be finemapped
#'
#' @param weight_list a list of weights for each gene
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
#' @importFrom logging loginfo
#'
#' @return finemapping results.
#'
#' @export
#'
finemap_region <- function(z_snp,
                           z_gene,
                           gene_info,
                           regionlist,
                           region_tag,
                           region_info = NULL,
                           weight_list = NULL,
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
                           ...){

  loginfo("Run finemapping with L = %d for region %s", L, region_tag)

  # get regionlist with full SNPs if not available
  if (missing(regionlist)) {
    loginfo("Get regionlist for region %s", region_tag)
    regioninfo <- region_info[region_info$region_tag == region_tag, ]
    res <- get_regionlist(regioninfo, gene_info, adjust_boundary = FALSE)
    regionlist <- res$regionlist
    # boundary_genes <- res$boundary_genes
    rm(res)
  }

  # combine z-scores of SNPs and genes
  z_snp$type <- "SNP"
  z_snp$QTLtype <- "SNP"
  if (is.null(z_gene$type)){
    z_gene$type <- "gene"
  }
  if (is.null(z_gene$QTLtype)){
    z_gene$QTLtype <- "gene"
  }
  zdf <- rbind(z_snp[, c("id", "z", "type", "QTLtype")],
               z_gene[, c("id", "z", "type", "QTLtype")])

  types <- unique(zdf$type)

  # priors
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

  # run susie for this region
  gid <- regionlist[[region_tag]][["gid"]]
  sid <- regionlist[[region_tag]][["sid"]]
  g_type <- zdf$type[match(gid, zdf$id)]
  s_type <- zdf$type[match(sid, zdf$id)]
  gs_type <- c(g_type, s_type)
  # g_QTLtype <- zdf$QTLtype[match(gid, zdf$id)]

  p <- length(gid) + length(sid)

  if (any(is.na(pi_prior))){
    prior <- rep(1/p, p)
  } else {
    prior <- unname(pi_prior[gs_type])
  }

  if (any(is.na(V_prior))){
    V <- matrix(rep(50, L * p), nrow = L)
    # following the default in susieR::susie_rss
  } else{
    V <- unname(V_prior[gs_type])
    V <- matrix(rep(V, each = L), nrow=L)
  }

  if (isTRUE(use_null_weight)){
    null_weight <- max(0, 1 - sum(prior))
    prior <- prior/(1-null_weight)
  } else {
    null_weight <- NULL
  }

  z.g <- zdf[match(gid, zdf$id), ][["z"]]
  z.s <- zdf[match(sid, zdf$id), ][["z"]]
  z <- c(z.g, z.s)

  # SNP information in this region
  ld_snpinfo <- do.call(rbind, lapply(regionlist[[region_tag]][["SNP_info"]],read_LD_SNP_file))
  # sidx <- match(sid, ld_snpinfo$id)

  # compute correlation matrices
  R_sg_file <- file.path(cor_dir, paste0(region_tag, ".R_snp_gene.RDS"))
  R_g_file <- file.path(cor_dir, paste0(region_tag,  ".R_gene.RDS"))
  R_s_file <- file.path(cor_dir, paste0(region_tag, ".R_snp.RDS"))

  if (isTRUE(force_compute_cor)) {
    # force compute correlation matrix
    res <- compute_region_cor(regionlist, region_tag, weight_list)
    R_snp <- res$R_snp
    R_snp_gene <- res$R_snp_gene
    R_gene <- res$R_gene
    # R_snp_gene <- R_snp_gene[sidx, , drop = F]
    rm(res)

    if (isTRUE(save_cor)) {
      loginfo("Save correlation matrices to %s", cor_dir)
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
      loginfo("Load precomputed correlation matrices from %s", cor_dir)
      R_snp_gene <- read_LD(R_sg_file)
      R_gene <- read_LD(R_g_file)
      R_snp <- read_LD(R_s_file)
      # R_snp_gene <- R_snp_gene[sidx, , drop = F]

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
      res <- compute_region_cor(regionlist, region_tag, weight_list)
      R_snp <- res$R_snp
      R_snp_gene <- res$R_snp_gene
      R_gene <- res$R_gene
      # R_snp_gene <- R_snp_gene[sidx, , drop = F]
      rm(res)

      if (isTRUE(save_cor)) {
        loginfo("Save correlation matrices to %s", cor_dir)
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

  # run susie
  # in susie, prior_variance is under standardized scale (if performed)
  loginfo("run susie_rss ...")
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

  # annotate susie result with SNP and gene information
  loginfo("annotate susie result ...")
  susie_res_df <- anno_susie(susie_res,
                             geneinfo = gene_info,
                             snpinfo = ld_snpinfo,
                             gid = gid,
                             sid = sid,
                             zdf = zdf,
                             region_tag = region_tag)

  return(susie_res_df)

}


#' Run cTWAS version of susie_rss for a single region
ctwas_susie_rss <- function(z,
                            R,
                            prior_weights = NULL,
                            prior_variance = NULL,
                            L = 5,
                            z_ld_weight = 0,
                            null_weight = NULL,
                            coverage = 0.95,
                            min_abs_corr = 0.5,
                            max_iter = 100,
                            ...){

  if (missing(R)) {
    if (L == 1){
      # R does not matter for susie when L = 1
      R <- diag(length(z))
    } else {
      stop("R (correlation matrix) is required when L > 1")
    }
  }

  # in susie, prior_variance is under standardized scale (if performed)
  susie_res <- susie_rss(z,
                         R,
                         prior_weights = prior_weights,
                         prior_variance = prior_variance,
                         estimate_prior_variance = F,
                         L = L,
                         z_ld_weight = z_ld_weight,
                         null_weight = null_weight,
                         coverage = coverage,
                         min_abs_corr = min_abs_corr,
                         max_iter = max_iter,
                         ...)

  return(susie_res)
}

# annotate susie results with SNP and gene information
anno_susie <- function(susie_res,
                       geneinfo,
                       snpinfo,
                       gid,
                       sid,
                       zdf,
                       region_tag) {

  gidx <- match(gid, geneinfo$id)
  sidx <- match(sid, snpinfo$id)
  g_type <- zdf$type[match(gid, zdf$id)]
  g_QTLtype <- zdf$QTLtype[match(gid, zdf$id)]

  if (length(geneinfo) != 0) {
    gene_anno <- data.frame(geneinfo[gidx,  c("chrom", "id", "p0")], type = g_type, QTLtype = g_QTLtype)
    colnames(gene_anno) <-  c("chrom", "id", "pos", "type", "QTLtype")
  } else {
    gene_anno <- NULL
  }

  snp_anno <- data.frame(snpinfo[sidx, c("chrom", "id", "pos")], type = "SNP", QTLtype = "SNP")
  colnames(snp_anno) <-  c("chrom", "id", "pos", "type", "QTLtype")

  anno <- as.data.frame(rbind(gene_anno, snp_anno))
  susie_res_df <- cbind(anno, region_tag = region_tag, susie_pip = susie_res$pip)

  p <- length(gid) + length(sid)
  susie_res_df$mu2 <- colSums(susie_res$mu2[, seq(1, p)[1:p!=susie_res$null_index], drop = F]) #WARN: not sure for L>1

  susie_res_df$cs_index <- 0
  if (!is.null(susie_res$sets$cs)){
    for (cs_i in susie_res$sets$cs_index){
      X.idx <- susie_res$sets$cs[[paste0("L", cs_i)]]
      X.idx <- X.idx[X.idx != susie_res$null_index] # susie_rss' bug
      susie_res_df$cs_index[X.idx] <- cs_i
      #TODO: note this ignores the fact that some variants can belong to multiple CS
    }
  }

  return(susie_res_df)
}

