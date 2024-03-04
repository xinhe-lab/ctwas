
#' Compute correlation matrices for a region
#'
#' @param R_snp Reference LD matrix
#' @param ld_snpinfo SNP information of the SNPs in \code{R_snp} matrix
#' @param region_idx region indices of the region from regionlist
#' @param wgtlist a list of weights for genes in the region
#' @param save TRUE/FALSE, if TRUE, save correlation matrices in \code{outputdir}
#' @param outputdir directory of correlation matrix result
#' @param outname filename of correlation result
#' @param verbose TRUE/FALSE, if TRUE,
#'
#' @return a list of correlation matrices (R_snp, R_snp_gene, R_gene), SNP indices and gene id
#' @export
compute_region_cor <- function(R_snp,
                               ld_snpinfo,
                               region_idx,
                               wgtlist,
                               save = FALSE,
                               outputdir = NULL,
                               outname = NULL,
                               verbose = FALSE) {

  sidx <- match(region_idx[["sid"]], ld_snpinfo$id)
  gnames <- region_idx[["gid"]]

  # subset wgtlist for genes in this region
  wgtlist <- wgtlist[gnames]

  if (isTRUE(verbose)){
    loginfo("R_snp: %d rows, %d columns \n", nrow(R_snp), ncol(R_snp))
    loginfo("%d SNPs, %d genes in the region.", length(sidx), length(gnames))
  }

  R_snp_gene <- matrix(NA, nrow(R_snp), length(gnames))
  R_gene <- diag(length(gnames))

  if (length(gnames) > 0) {
    ldr <- list()

    # compute SNP-gene correlation matrix
    if (isTRUE(verbose))
      loginfo("Compute SNP-gene correlation matrix")

    for (i in 1:length(gnames)){
      gname <- gnames[i]
      wgt <- wgtlist[[gname]]
      snpnames <- rownames(wgt)
      ld.idx <- match(snpnames, ld_snpinfo$id)
      ldr[[gname]] <- ld.idx
      R.s <- R_snp[ld.idx, ld.idx]
      R_snp_gene[,i] <- sapply(1:nrow(R_snp),
                               function(x){t(wgt)%*%R_snp[ld.idx,x]/sqrt(t(wgt)%*%R.s%*%wgt*R_snp[x,x])})
    }

    # compute gene-gene correlation matrix
    if (length(gnames) > 1){
      if (isTRUE(verbose))
        loginfo("Compute gene-gene correlation matrix")

      gene_pairs <- combn(length(gnames), 2)
      wgtr <- wgtlist[gnames]
      gene_corrs <- apply(gene_pairs, 2, function(x){t(wgtr[[x[1]]])%*%R_snp[ldr[[x[1]]], ldr[[x[2]]]]%*%wgtr[[x[2]]]/(
        sqrt(t(wgtr[[x[1]]])%*%R_snp[ldr[[x[1]]], ldr[[x[1]]]]%*%wgtr[[x[1]]]) *
          sqrt(t(wgtr[[x[2]]])%*%R_snp[ldr[[x[2]]], ldr[[x[2]]]]%*%wgtr[[x[2]]]))})
      R_gene[t(gene_pairs)] <- gene_corrs
      R_gene[t(gene_pairs[c(2,1),])] <- gene_corrs
    }
  }

  # extract R_snp subset by sidx
  R_snp <- R_snp[sidx, sidx, drop = F]

  # R_snp_gene <- R_snp_gene[sidx, , drop = F]

  res <- list("R_snp" = R_snp,
              "R_snp_gene" = R_snp_gene,
              "R_gene" = R_gene,
              "snp_info" = ld_snpinfo[sidx,],
              "gene_id" = region_idx$gid)

  # save correlation matrices
  if (isTRUE(save)) {
    cor_dir <- outputdir
    cor_file <- paste0(outname, ".RDS")
    if (isTRUE(verbose)){
      loginfo("save correlation matrices to %s", file.path(cor_dir, cor_file))

      loginfo("R_snp: %d rows, %d columns, size: %s \n",
              nrow(R_snp), ncol(R_snp), format(object.size(R_snp), units = "Mb"))

      loginfo("R_snp_gene: %d rows, %d columns, size: %s \n",
              nrow(R_snp_gene), ncol(R_snp_gene), format(object.size(R_snp_gene), units = "Mb"))

      loginfo("R_gene: %d rows, %d columns, size: %s",
              nrow(R_gene), ncol(R_gene), format(object.size(R_gene), units = "Mb"))
    }

    if (!dir.exists(cor_dir))
      dir.create(cor_dir, recursive = TRUE)

    saveRDS(res, file = file.path(cor_dir, cor_file))
    res$cor_dir <- cor_dir
    res$cor_file <- cor_file
  }

  return(res)

}
