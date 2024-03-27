
#' Compute correlation matrices for a region
compute_region_cor <- function(sids, gids, R_snp, weights, snpinfo) {

  # subset to genes in this region
  weights <- weights[gids]
  # extract wgt from weights
  wgtlist <- lapply(weights, "[[", "wgt")
  names(wgtlist) <- names(weights)

  # compute correlation matrices
  R_snp_gene <- matrix(NA, nrow(R_snp), length(gids))
  R_gene <- diag(length(gids))

  if (length(gids) > 0) {
    ldr <- list()
    # compute SNP-gene correlation matrix
    # loginfo("compute SNP-gene correlation matrix")
    for (i in 1:length(gids)){
      gid <- gids[i]
      wgt <- wgtlist[[gid]]
      snpnames <- intersect(rownames(wgt), snpinfo$id)
      wgt <- wgt[snpnames,,drop=F]
      ld.idx <- match(snpnames, snpinfo$id)
      ldr[[gid]] <- ld.idx
      R.s <- R_snp[ld.idx, ld.idx]
      R_snp_gene[,i] <- sapply(1:nrow(R_snp),
                               function(x){t(wgt)%*%R_snp[ld.idx,x]/sqrt(t(wgt)%*%R.s%*%wgt*R_snp[x,x])})
    }

    # compute gene-gene correlation matrix
    # loginfo("compute gene-gene correlation matrix")
    if (length(gids) > 1){
      gene_pairs <- combn(length(gids), 2)
      wgtr <- wgtlist[gids]
      gene_corrs <- apply(gene_pairs, 2, function(x){t(wgtr[[x[1]]])%*%R_snp[ldr[[x[1]]], ldr[[x[2]]]]%*%wgtr[[x[2]]]/(
        sqrt(t(wgtr[[x[1]]])%*%R_snp[ldr[[x[1]]], ldr[[x[1]]]]%*%wgtr[[x[1]]]) *
          sqrt(t(wgtr[[x[2]]])%*%R_snp[ldr[[x[2]]], ldr[[x[2]]]]%*%wgtr[[x[2]]]))})
      R_gene[t(gene_pairs)] <- gene_corrs
      R_gene[t(gene_pairs[c(2,1),])] <- gene_corrs
    }
  }

  # subset R_snp and R_snp_gene by sidx
  sidx <- match(sids, snpinfo$id)
  R_snp <- R_snp[sidx, sidx, drop = F]
  R_snp_gene <- R_snp_gene[sidx, , drop = F]

  if (anyNA(R_snp))
    stop("R_snp matrix contains missing values!\n")

  if (anyNA(R_snp_gene))
    stop("R_snp_gene matrix contains missing values!\n")

  if (anyNA(R_gene))
    stop("R_gene matrix contains missing values!\n")

  # loginfo("%d genes, %d SNPs in the correlation matrices", ncol(R_snp_gene), nrow(R_snp_gene))

  return(list("R_snp" = R_snp,
              "R_snp_gene" = R_snp_gene,
              "R_gene" = R_gene))
}
