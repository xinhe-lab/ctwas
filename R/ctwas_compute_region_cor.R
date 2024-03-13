
#' Compute correlation matrices for a region
#'
#' @param regionlist regionlist
#' @param region_tag region tag
#' @param weight_list weight list
#'
#' @return a list of correlation matrices
#' @export
compute_region_cor <- function(regionlist,
                               region_tag,
                               weight_list) {

  loginfo("Compute correlation matrices for region %s ...", region_tag)

  R_snp <- lapply(regionlist[[region_tag]][["LD_matrix"]], read_LD)
  if (length(R_snp)==1){
    R_snp <- unname(R_snp[[1]])
  } else {
    R_snp <- suppressWarnings(as.matrix(Matrix::bdiag(R_snp)))
  }

  ld_snpinfo <- lapply(regionlist[[region_tag]][["SNP_info"]],read_LD_SNP_file)
  if (length(ld_snpinfo) > 1){
    ld_snpinfo <- do.call(rbind,ld_snpinfo)
  }

  sidx <- match(regionlist[[region_tag]][["sid"]], ld_snpinfo$id)
  gnames <- regionlist[[region_tag]][["gid"]]

  # subset weight_list for genes in this region
  weight_list <- weight_list[gnames]

  loginfo("%d SNPs, %d genes in the region.", length(sidx), length(gnames))

  R_snp_gene <- matrix(NA, nrow(R_snp), length(gnames))
  R_gene <- diag(length(gnames))

  if (length(gnames) > 0) {
    ldr <- list()

    # compute SNP-gene correlation matrix
    loginfo("Compute SNP-gene correlation matrix")

    for (i in 1:length(gnames)){
      gname <- gnames[i]
      wgt <- weight_list[[gname]]
      snpnames <- rownames(wgt)
      ld.idx <- match(snpnames, ld_snpinfo$id)
      ldr[[gname]] <- ld.idx
      R.s <- R_snp[ld.idx, ld.idx]
      R_snp_gene[,i] <- sapply(1:nrow(R_snp),
                               function(x){t(wgt)%*%R_snp[ld.idx,x]/sqrt(t(wgt)%*%R.s%*%wgt*R_snp[x,x])})
    }

    # compute gene-gene correlation matrix
    if (length(gnames) > 1){
      loginfo("Compute gene-gene correlation matrix")

      gene_pairs <- combn(length(gnames), 2)
      wgtr <- weight_list[gnames]
      gene_corrs <- apply(gene_pairs, 2, function(x){t(wgtr[[x[1]]])%*%R_snp[ldr[[x[1]]], ldr[[x[2]]]]%*%wgtr[[x[2]]]/(
        sqrt(t(wgtr[[x[1]]])%*%R_snp[ldr[[x[1]]], ldr[[x[1]]]]%*%wgtr[[x[1]]]) *
          sqrt(t(wgtr[[x[2]]])%*%R_snp[ldr[[x[2]]], ldr[[x[2]]]]%*%wgtr[[x[2]]]))})
      R_gene[t(gene_pairs)] <- gene_corrs
      R_gene[t(gene_pairs[c(2,1),])] <- gene_corrs
    }
  }

  # subset R_snp by sidx
  R_snp <- R_snp[sidx, sidx, drop = F]

  # subset R_snp_gene by sidx
  R_snp_gene <- R_snp_gene[sidx, , drop = F]

  return(list("R_snp" = R_snp,
              "R_snp_gene" = R_snp_gene,
              "R_gene" = R_gene))
}
