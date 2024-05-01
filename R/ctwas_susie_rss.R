
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

#' annotate susie results with SNP and gene information
anno_susie <- function(susie_res,
                       gid,
                       sid,
                       region_id,
                       z = NULL,
                       g_type = "gene",
                       g_context = "gene",
                       g_group = "gene",
                       geneinfo = NULL,
                       snpinfo = NULL,
                       include_cs_index = TRUE) {

  if (!is.null(geneinfo)) {
    gidx <- match(gid, geneinfo$id)
    gene_anno <- data.frame(geneinfo[gidx,  c("chrom", "p0", "id")],
                            type = g_type, context = g_context, group = g_group)
    colnames(gene_anno) <-  c("chrom", "pos", "id", "type", "context", "group")
  } else {
    gene_anno <- data.frame(id = gid, type = g_type, context = g_context, group = g_group)
  }

  if (!is.null(snpinfo)) {
    sidx <- match(sid, snpinfo$id)
    snp_anno <- data.frame(snpinfo[sidx, c("chrom", "pos", "id")],
                           type = "SNP", context = "SNP", group = "SNP")
    colnames(snp_anno) <-  c("chrom", "pos", "id", "type", "context", "group")
  } else {
    snp_anno <- data.frame(id = sid, type = "SNP", context = "SNP", group = "SNP")
  }

  susie_res_df <- as.data.frame(rbind(gene_anno, snp_anno))

  if (!is.null(z)) {
    susie_res_df$z <- z
  }

  susie_res_df$region_id <- region_id

  susie_res_df$susie_pip <- susie_res$pip

  p <- length(gid) + length(sid)
  susie_res_df$mu2 <- colSums(susie_res$mu2[, seq(1, p)[1:p!=susie_res$null_index], drop = F]) #WARN: not sure for L>1

  if (isTRUE(include_cs_index)) {
    susie_res_df$cs_index <- 0
    if (!is.null(susie_res$sets$cs)){
      for (cs_i in susie_res$sets$cs_index){
        X.idx <- susie_res$sets$cs[[paste0("L", cs_i)]]
        X.idx <- X.idx[X.idx != susie_res$null_index] # susie_rss' bug
        susie_res_df$cs_index[X.idx] <- cs_i
        #TODO: note this ignores the fact that some variants can belong to multiple CS
      }
    }
  }

  return(susie_res_df)
}
