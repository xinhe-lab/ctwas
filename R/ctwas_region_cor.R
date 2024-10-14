
#' @title Gets correlation matrices for a single region.
#'
#' @description Loads precomputed correlation matrices if available in \code{cor_dir},
#' otherwise, computes correlation matrices if no precomputed correlation matrices,
#' or force_compute_cor = TRUE.
#'
#' @param region_id a character string of region id.
#'
#' @param sids SNP IDs in the region
#'
#' @param gids gene IDs in the region
#'
#' @param LD_map a data frame with filenames of LD matrices and SNP information for the regions.
#'
#' @param weights a list of preprocessed weights
#'
#' @param force_compute_cor If TRUE, force computing correlation (R) matrices
#'
#' @param save_cor If TRUE, save correlation (R) matrices to \code{cor_dir}
#'
#' @param cor_dir a string, the directory to store correlation (R) matrices
#'
#' @param LD_format file format for LD matrix. If "custom", use a user defined
#' \code{LD_loader_fun()} function to load LD matrix.
#'
#' @param LD_loader_fun a user defined function to load LD matrix when \code{LD_format = "custom"}.
#'
#' @param snpinfo_loader_fun a user defined function to load SNP information file,
#' if SNP information files are not in standard cTWAS reference format.
#'
#' @param verbose If TRUE, print detail messages
#'
#' @return correlation matrices (R_snp, R_snp_gene and R_gene)
#'
#' @importFrom logging loginfo
#' @importFrom Matrix bdiag
#'
#' @export
get_region_cor <- function(region_id,
                           sids,
                           gids,
                           LD_map,
                           weights,
                           force_compute_cor = FALSE,
                           save_cor = FALSE,
                           cor_dir = NULL,
                           LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
                           LD_loader_fun = NULL,
                           snpinfo_loader_fun = NULL,
                           verbose = FALSE) {

  LD_format <- match.arg(LD_format)

  if (save_cor) {
    if (is.null(cor_dir))
      stop("cor_dir is required for saving correlation matrices")
  }

  if (!is.null(cor_dir)) {
    if (!dir.exists(cor_dir))
      dir.create(cor_dir, showWarnings = FALSE, recursive = TRUE)

    R_sg_file <- file.path(cor_dir, paste0("region.", region_id, ".R_snp_gene.RDS"))
    R_g_file <- file.path(cor_dir, paste0("region.", region_id, ".R_gene.RDS"))
    R_s_file <- file.path(cor_dir, paste0("region.", region_id, ".R_snp.RDS"))
  }

  cor_files_exist <- isTRUE(!is.null(cor_dir) && file.exists(R_sg_file) && file.exists(R_g_file) && file.exists(R_s_file))

  if (cor_files_exist && !force_compute_cor) {
    if (verbose)
      loginfo("Load precomputed correlation matrices for region %s", region_id)

    # load precomputed correlation matrices
    R_snp_gene <- load_LD(R_sg_file)
    R_gene <- load_LD(R_g_file)
    R_snp <- load_LD(R_s_file)
  } else {
    # if no precomputed correlation matrices, or force_compute_cor = TRUE,
    # compute correlation matrices
    if (missing(sids))
      stop("'sids' is required for computing correlation matrices")

    if (missing(gids))
      stop("'gids' is required for computing correlation matrices")

    if (missing(weights))
      stop("'weights' is required for computing correlation matrices!")

    if (!inherits(weights,"list"))
      stop("'weights' should be a list!")

    if (missing(LD_map))
      stop("LD_map is required for computing correlation matrices")

    # load LD matrix of the region
    LD_matrix_files <- unlist(strsplit(LD_map$LD_file[LD_map$region_id == region_id], split = ","))
    stopifnot(all(file.exists(LD_matrix_files)))
    if (verbose)
      loginfo("Load LD matrix for region %s", region_id)
    if (length(LD_matrix_files) > 1) {
      R_snp <- lapply(LD_matrix_files, load_LD, format = LD_format, LD_loader_fun = LD_loader_fun)
      R_snp <- suppressWarnings(as.matrix(bdiag(R_snp)))
    } else {
      R_snp <- load_LD(LD_matrix_files, format = LD_format, LD_loader_fun = LD_loader_fun)
    }

    # load SNP info of the region
    if (verbose)
      loginfo("Load SNP info for region %s", region_id)
    SNP_info_files <- unlist(strsplit(LD_map$SNP_file[LD_map$region_id == region_id], split = ","))
    stopifnot(all(file.exists(SNP_info_files)))
    snpinfo <- read_snp_info_files(SNP_info_files, snpinfo_loader_fun = snpinfo_loader_fun)

    # compute correlation matrices
    if (verbose)
      loginfo("Compute correlation matrices")
    res <- compute_region_cor(sids, gids, R_snp, snpinfo$id, weights)
    R_snp <- res$R_snp
    R_snp_gene <- res$R_snp_gene
    R_gene <- res$R_gene
    rm(res)
    # save correlation matrices
    if (isTRUE(save_cor && !is.null(cor_dir))) {
      if (verbose)
        loginfo("Save correlation matrices")
      saveRDS(R_snp_gene, file=R_sg_file)
      saveRDS(R_gene, file=R_g_file)
      saveRDS(R_snp, file=R_s_file)
    }
  }
  return(list("R_snp" = R_snp,
              "R_snp_gene" = R_snp_gene,
              "R_gene" = R_gene))
}


#' @title Loads precomputed correlation matrices for a single region.
#'
#' @description It loads precomputed correlation matrices for a single region.
#' It could load correlation matrices by \code{region_id} and
#' directory of correlation matrices \code{cor_dir}. Otherwise,
#' it loads correlation matrices by the
#' filenames (\code{R_sg_file}, \code{R_sg_file}, \code{R_s_file})
#' if they are provided.
#'
#' @param region_id a character string of region id.
#'
#' @param cor_dir a string, the directory to store correlation matrices.
#'
#' @param R_sg_file filename of SNP-gene correlations.
#'
#' @param R_g_file filename of gene-gene correlations.
#'
#' @param R_s_file filename of SNP-SNP correlations.
#'
#' @return correlation matrices (R_snp, R_snp_gene and R_gene)
#'
#' @export
load_region_cor <- function(region_id,
                            cor_dir,
                            R_sg_file,
                            R_g_file,
                            R_s_file) {

  if (missing(R_sg_file)){
    R_sg_file <- file.path(cor_dir, paste0("region.", region_id, ".R_snp_gene.RDS"))
  }

  if (missing(R_g_file)){
    R_g_file <- file.path(cor_dir, paste0("region.", region_id, ".R_gene.RDS"))
  }

  if (missing(R_s_file)){
    R_s_file <- file.path(cor_dir, paste0("region.", region_id, ".R_snp.RDS"))
  }

  # load precomputed correlation matrices
  R_snp_gene <- load_LD(R_sg_file)
  R_gene <- load_LD(R_g_file)
  R_snp <- load_LD(R_s_file)

  return(list("R_snp" = R_snp,
              "R_snp_gene" = R_snp_gene,
              "R_gene" = R_gene))
}

#' @title Computes correlation matrices for a single region
#'
#' @param sids SNP IDs in the region
#'
#' @param gids gene IDs in the region
#'
#' @param R_snp LD (R) matrix
#'
#' @param LD_sids SNP IDs for the rows and columns of the LD matrix
#'
#' @param weights a list of weights for all the genes
#'
#' @return a list of correlation matrices (R_snp, R_snp_gene and R_gene)
#'
#' @importFrom utils combn
compute_region_cor <- function(sids, gids, R_snp, LD_sids, weights) {

  # check input data
  if (!is.list(weights)){
    stop("'weights' should be a list.")
  }

  R_snp <- as.matrix(R_snp)

  # subset weights to genes in this region
  weights <- weights[gids]

  # extract wgtlist and filter by SNPs in LD
  wgtlist <- lapply(weights, function(x){
    wgt <- x$wgt
    wgt[rownames(wgt) %in% LD_sids, , drop=FALSE]
  })
  names(wgtlist) <- names(weights)

  # subset LD matrix to SNPs (sids) in the region and SNPs in weights
  all_wgt_snps <- unique(as.character(unlist(lapply(wgtlist, function(x){rownames(x)}))))
  sidx <- which(LD_sids %in% c(sids, all_wgt_snps))
  LD_sids <- LD_sids[sidx]
  R_snp <- R_snp[sidx, sidx]
  rm(all_wgt_snps)
  rm(sidx)

  # compute correlation matrices
  R_snp_gene <- matrix(NA, nrow = nrow(R_snp), ncol = length(gids))
  R_gene <- diag(length(gids))

  if (length(gids) > 0) {
    wgt_sidx_list <- list()
    # compute SNP-gene correlation matrix
    for (i in 1:length(gids)){
      gid <- gids[i]
      wgt <- wgtlist[[gid]]
      wgt_snps <- rownames(wgt)
      wgt_sidx <- match(wgt_snps, LD_sids)
      wgt_sidx_list[[gid]] <- wgt_sidx
      R.s <- R_snp[wgt_sidx, wgt_sidx]
      R_snp_gene[,i] <- sapply(1:nrow(R_snp),
                               function(x){t(wgt)%*%R_snp[wgt_sidx,x]/sqrt(t(wgt)%*%R.s%*%wgt*R_snp[x,x])})
    }

    # compute gene-gene correlation matrix
    if (length(gids) > 1){
      gene_pairs <- combn(length(gids), 2)
      gene_corrs <- apply(gene_pairs, 2, function(x){
        wgt1 <- wgtlist[[x[1]]]
        wgt2 <- wgtlist[[x[2]]]
        wgt1_sidx <- wgt_sidx_list[[x[1]]]
        wgt2_sidx <- wgt_sidx_list[[x[2]]]
        t(wgt1)%*%R_snp[wgt1_sidx, wgt2_sidx]%*%wgt2/(
          sqrt(t(wgt1)%*%R_snp[wgt1_sidx, wgt1_sidx]%*%wgt1) *
            sqrt(t(wgt2)%*%R_snp[wgt2_sidx, wgt2_sidx]%*%wgt2))
      })
      R_gene[t(gene_pairs)] <- gene_corrs
      R_gene[t(gene_pairs[c(2,1),])] <- gene_corrs
    }
  }

  # subset R_snp and R_snp_gene by sidx
  sidx <- match(sids, LD_sids)
  R_snp <- R_snp[sidx, sidx, drop = FALSE]
  R_snp_gene <- R_snp_gene[sidx, , drop = FALSE]

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
