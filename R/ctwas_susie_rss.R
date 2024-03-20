
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
