
#' Compute correlation matrices for a single region
#'
#' @param sids SNP IDs
#'
#' @param gids gene IDs
#'
#' @param R_snp LD (R) matrix
#'
#' @param LD_sids SNP IDs for the rows and columns of the LD matrix
#'
#' @param weights a list of weights for all the genes
#'
#' @return a list of correlation matrices (R_snp, R_snp_gene and R_gene)
#'
#' @export
#'
compute_region_cor <- function(sids, gids, R_snp, LD_sids, weights) {

  # check input data
  if (!is.list(weights)){
    stop("'weights' should be a list.")
  }

  # subset weights to genes in this region
  weights <- weights[gids]

  # extract wgtlist and filter by SNPs in LD
  wgtlist <- lapply(weights, function(x){
    wgt <- x$wgt
    wgt[rownames(wgt) %in% LD_sids, , drop=FALSE]
    })
  names(wgtlist) <- names(weights)

  # compute correlation matrices
  R_snp_gene <- matrix(NA, nrow = nrow(R_snp), ncol = length(gids))
  R_gene <- diag(length(gids))

  if (length(gids) > 0) {
    ldr <- list()
    # compute SNP-gene correlation matrix
    for (i in 1:length(gids)){
      gid <- gids[i]
      wgt <- wgtlist[[gid]]
      snpnames <- rownames(wgt)
      ld.idx <- match(snpnames, LD_sids)
      ldr[[gid]] <- ld.idx
      R.s <- R_snp[ld.idx, ld.idx]
      R_snp_gene[,i] <- sapply(1:nrow(R_snp),
                               function(x){t(wgt)%*%R_snp[ld.idx,x]/sqrt(t(wgt)%*%R.s%*%wgt*R_snp[x,x])})
    }

    # compute gene-gene correlation matrix
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
  sidx <- match(sids, LD_sids)
  R_snp <- R_snp[sidx, sidx, drop = F]
  R_snp_gene <- R_snp_gene[sidx, , drop = F]

  # add rownames and colnames
  rownames(R_snp) <- sids
  colnames(R_snp) <- sids
  rownames(R_snp_gene) <- sids
  colnames(R_snp_gene) <- gids
  rownames(R_gene) <- gids
  colnames(R_gene) <- gids

  if (anyNA(R_snp))
    stop("R_snp matrix contains missing values!\n")

  if (anyNA(R_snp_gene))
    stop("R_snp_gene matrix contains missing values!\n")

  if (anyNA(R_gene))
    stop("R_gene matrix contains missing values!\n")

  return(list("R_snp" = R_snp,
              "R_snp_gene" = R_snp_gene,
              "R_gene" = R_gene))
}
