
compute_region_cor <- function(regionlist,
                               wgtlist,
                               b,
                               rn,
                               save = FALSE,
                               outputdir = NULL,
                               outname = NULL) {

  loginfo("Compute correlation for region %s_%s", b, rn)

  region_idx <- regionlist[[b]][[rn]]
  gnames <- region_idx[["gid"]]

  # subset wgtlist for genes in this region
  wgtlist <- wgtlist[gnames]

  # load precomputed LD matrix
  R_snp <- lapply(region_idx[["LD_matrix"]], ctwas:::read_LD)

  if (length(R_snp)==1){
    R_snp <- unname(R_snp[[1]])
  } else {
    R_snp <- suppressWarnings(as.matrix(Matrix::bdiag(R_snp)))
  }

  loginfo("R_snp (reference LD matrix): %d rows, %d columns, size: %s",
          nrow(R_snp), ncol(R_snp), format(object.size(R_snp), units = "Mb"))

  ld_snpinfo <- ctwas:::read_LD_SNP_file(region_idx[["SNP_info"]])
  sidx <-  match(region_idx[["sid"]], ld_snpinfo$id)
  loginfo("%d SNPs, %d genes in the region.", length(sidx), length(gnames))

  R_snp_gene <- matrix(NA, nrow(R_snp), length(gnames))
  R_gene <- diag(length(gnames))

  if (length(gnames) > 0) {
    ldr <- list()

    # compute SNP-gene correlation matrix
    loginfo("Compute SNP-gene correlation matrix")
    for (i in 1:length(gnames)){
      gname <- gnames[i]
      wgt <- wgtlist[[gname]]
      snpnames <- rownames(wgt)
      ld.idx <- match(snpnames, ld_snpinfo$id)
      ldr[[gname]] <- ld.idx
      R.s <- R_snp[ld.idx, ld.idx]
      R_snp_gene[,i] <- sapply(1:nrow(R_snp), function(x){t(wgt)%*%R_snp[ld.idx,x]/sqrt(t(wgt)%*%R.s%*%wgt*R_snp[x,x])})
    }

    # compute gene-gene correlation matrix
    if (length(gnames) > 1){
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

  res <- list("R_snp" = R_snp,
              "R_snp_gene" = R_snp_gene,
              "R_gene" = R_gene,
              "sidx" = sidx,
              "gid" = gnames,
              "LD_matrix" = region_idx[["LD_matrix"]],
              "SNP_info" = region_idx[["SNP_info"]],
              "wgtlist" = wgtlist)

  # summary
  loginfo("R_snp (subset by sidx): %d rows, %d columns, size: %s \n",
              nrow(R_snp), ncol(R_snp), format(object.size(R_snp), units = "Mb"))

  loginfo("R_snp_gene: %d rows, %d columns, size: %s \n",
              nrow(R_snp_gene), ncol(R_snp_gene), format(object.size(R_snp_gene), units = "Mb"))

  loginfo("R_gene: %d rows, %d columns, size: %s",
              nrow(R_gene), ncol(R_gene), format(object.size(R_gene), units = "Mb"))

  # save correlation matrices
  if (isTRUE(save)) {
    cor_dir <- outputdir
    cor_file <- paste0(outname, ".chr", b, "_region", rn, ".R.matrices.RDS")
    res$cor_dir <- cor_dir
    res$cor_file <- cor_file
    if (!dir.exists(cor_dir)) {
      dir.create(cor_dir, recursive = TRUE)
    }
    loginfo("save correlation matrices to %s", file.path(cor_dir, cor_file))
    save(res, file = file.path(cor_dir, cor_file))
  }

  return(res)

}
