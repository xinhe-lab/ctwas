#' run cTWAS finemapping for a single region
#'
#' @param z_snp A data frame with four columns: "id", "A1", "A2", "z".
#' giving the z scores for snps. "A1" is effect allele. "A2" is the other allele.
#'
#' @param z_gene A data frame with two columns: "id", "z". giving the z scores for genes.
#' Optionally, a "type" column can also be supplied; this is for using multiple sets of weights
#'
#' @param region_info a data frame of region definition and associated file names
#'
#' @param gene_info a data frame of gene information obtained from \code{compute_gene_z}
#'
#' @param weight a string, pointing to a directory with the weights (.db file) in predictdb format.
#' A vector of multiple sets of weights in PredictDB format can also be specified;
#' genes will have their filename appended to their gene name to ensure IDs are unique.
#'
#' @param regionlist the list of regions with SNP and gene indices.
#'
#' @param region_tag a character string of region tag to be finemapped
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
#' @param ncore The number of cores used to parallelize susie over regions
#'
#' @param logfile the log file, if NULL will print log info on screen
#'
#' @param save_LD_R TRUE/FALSE. If TRUE, save correlation (R) matrices
#'
#' @param outputdir a string, the directory to store output
#'
#' @param outname a string, the output name
#'
#' @importFrom logging addHandler loginfo writeToFile
#'
#' @return finemapping results.
#'
#' @export
#'
finemap_region <- function(z_snp,
                           z_gene,
                           regionlist,
                           region_info,
                           gene_info,
                           region_tag,
                           L = 5,
                           group_prior = NULL,
                           group_prior_var = NULL,
                           use_null_weight = TRUE,
                           coverage = 0.95,
                           min_abs_corr = 0.5,
                           logfile = NULL,
                           save_LD_R = FALSE){

  if (!is.null(logfile)){
    addHandler(writeToFile, file=logfile, level='DEBUG')
  }

  loginfo("Run finemapping with L = %d for region %s", L, region_tag)

  # combine z-scores of different types
  zdf <- combine_z(z_snp, z_gene)

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
  b <- as.numeric(unlist(strsplit(region_tag, split = ":"))[1])
  rn <- unlist(strsplit(region_tag, split = ":"))[2]

  gidx <- regionlist[[b]][[rn]][["gidx"]]
  sidx <- regionlist[[b]][[rn]][["sidx"]]
  gid <- regionlist[[b]][[rn]][["gid"]]
  sid <- regionlist[[b]][[rn]][["sid"]]
  g_type <- zdf$type[match(gid, zdf$id)]
  s_type <- zdf$type[match(sid, zdf$id)]
  gs_type <- c(g_type, s_type)
  g_QTLtype <- zdf$QTLtype[match(gid, zdf$id)]

  p <- length(gidx) + length(sidx)

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

  # prepare R matrix
  if (!("regRDS" %in% names(regionlist[[b]][[rn]]))){
    stop("R matrix info not available for region", region_tag)
  }

  regRDS <- regionlist[[b]][[rn]][["regRDS"]]

  if (file.exists(regRDS)){
    # load precomputed correlation matrix
    R <- read_LD(regRDS)
  } else {
    # compute correlation matrix
    if (L >1) {
      R <- compute_cor()
      if (isTRUE(save_LD_R)) {
        LD_R_file <- file.path(outputdir, paste0(outname, "_LDR"), paste0("chr", b, "_reg", rn, ".R.RDS"))
        saveRDS(R, LD_R_file)
        regionlist[[b]][[rn]][["regRDS"]] <- LD_R_file
      }
    }else{
      # R does not matter for susie when L = 1
      R <- diag(length(z))
    }
  }

  # in susie, prior_variance is under standardized scale (if performed)
  susie_res <- ctwas_susie_rss(z = z,
                               R = R,
                               prior_weights = prior,
                               prior_variance = V,
                               L = L,
                               null_weight = null_weight,
                               coverage = coverage,
                               min_abs_corr = min_abs_corr)

  # annotate susie results with SNP, gene information
  # geneinfo_chr <- read_exprvar(ld_exprvarfs[b])
  gene_info_chr <- gene_info[gene_info$chrom == b, ]
  if (length(gene_info_chr) !=0){
    gene_anno <- data.frame(gene_info_chr[gidx,  c("chrom", "id", "p0")], type = type, QTLtype = QTLtype)
    colnames(gene_anno) <-  c("chrom", "id", "pos", "type", "QTLtype")
  }else{
    gene_anno <- NULL
  }

  snp_info_region <- do.call(rbind, lapply(regRDS, read_ld_Rvar_RDS))
  snp_anno <- data.frame(snp_info_region[sidx, c("chrom", "id", "pos")], type = "SNP", QTLtype = "SNP")
  colnames(snp_anno) <-  c("chrom", "id", "pos", "type", "QTLtype")

  susie_res <- anno_susie(susie_res,
                          gene_anno = gene_anno,
                          snp_anno = snp_anno,
                          region_tag = region_tag)

  return(susie_res)

}


# Run susie_rss for a single region
ctwas_susie_rss <- function(z,
                            R,
                            prior_weights = NULL,
                            prior_variance = NULL,
                            L = 5,
                            z_ld_weight = 0,
                            null_weight = NULL,
                            coverage = 0.95,
                            min_abs_corr = 0.5){

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
                         min_abs_corr = min_abs_corr)

  return(susie_res)
}

# annotate susie results
anno_susie <- function(susie_res,
                       gene_anno,
                       snp_anno,
                       region_tag) {

  res <- as.data.frame(rbind(gene_anno, snp_anno))

  res$region_tag <- region_tag

  res$susie_pip <- susie_res$pip

  p <- length(gidx) + length(sidx)
  res$mu2 <- colSums(susie_res$mu2[, seq(1, p)[1:p!=susie_res$null_index], drop = F]) #WARN: not sure for L>1

  res$cs_index <- 0
  if (!is.null(susie_res$sets$cs)){
    for (cs_i in susie_res$sets$cs_index){
      X.idx <- susie_res$sets$cs[[paste0("L", cs_i)]]
      X.idx <- X.idx[X.idx != susie_res$null_index] # susie_rss' bug
      res$cs_index[X.idx] <- cs_i
      #TODO: note this ignores the fact that some variants can belong to multiple CS
    }
  }

  return(res)
}

