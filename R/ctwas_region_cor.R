
#' @title Gets correlation matrices for a single region.
#'
#' @description Loads precomputed correlation matrices if available in \code{cor_dir},
#' otherwise, computes correlation matrices if no precomputed correlation matrices,
#' or force_compute_cor = TRUE.
#'
#' @param region_id a character string of region id to be finemapped
#'
#' @param region_data a list object with data for the regions
#'
#' @param LD_map a data frame with filenames of LD matrices for each of the regions. Required when \code{use_LD = TRUE}.
#'
#' @param snp_map  list of data frames with SNP-to-region map for the reference. Required when \code{use_LD = TRUE}.
#'
#' @param weights a list of weights
#'
#' @param force_compute_cor TRUE/FALSE. If TRUE, force computing correlation (R) matrices
#'
#' @param save_cor TRUE/FALSE. If TRUE, save correlation (R) matrices to \code{cor_dir}
#'
#' @param cor_dir a string, the directory to store correlation (R) matrices
#'
#' @param LD_format file format for LD matrix. If "custom", use a user defined
#' \code{LD_loader_fun()} function to load LD matrix.
#'
#' @param LD_loader_fun a user defined function to load LD matrix when \code{LD_format = "custom"}.
#'
#' @param verbose TRUE/FALSE. If TRUE, print detail messages
#'
#' @return a list of correlation matrices (R_snp, R_snp_gene and R_gene)
#'
#' @importFrom logging loginfo
#' @importFrom Matrix bdiag
#'
#' @export
get_region_cor <- function(region_id,
                           region_data = NULL,
                           LD_map = NULL,
                           snp_map = NULL,
                           weights = NULL,
                           force_compute_cor = FALSE,
                           save_cor = FALSE,
                           cor_dir = NULL,
                           LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
                           LD_loader_fun,
                           verbose = FALSE) {

  LD_format <- match.arg(LD_format)

  if (!is.null(cor_dir)) {
    if (!dir.exists(cor_dir))
      dir.create(cor_dir, recursive = TRUE)
    R_sg_file <- file.path(cor_dir, paste0("region.", region_id, ".R_snp_gene.RDS"))
    R_g_file <- file.path(cor_dir, paste0("region.", region_id, ".R_gene.RDS"))
    R_s_file <- file.path(cor_dir, paste0("region.", region_id, ".R_snp.RDS"))
  }

  cor_files_exist <- isTRUE(!is.null(cor_dir) && file.exists(R_sg_file) && file.exists(R_g_file) && file.exists(R_s_file))

  if (cor_files_exist && !force_compute_cor) {
    if (verbose){
      loginfo("Load correlation matrices for region %s", region_id)
    }
    # load precomputed correlation matrices
    R_snp_gene <- load_LD(R_sg_file)
    R_gene <- load_LD(R_g_file)
    R_snp <- load_LD(R_s_file)
  } else {
    # if no precomputed correlation matrices, or force_compute_cor = TRUE,
    # compute correlation matrices
    if (verbose){
      loginfo("Compute correlation matrices for region %s", region_id)
    }
    if (is.null(region_data)) {
      stop("region_data is needed for computing correlation matrices")
    }
    # load LD matrix of the region
    if (is.null(LD_map) || is.null(snp_map)) {
      stop("LD_map and snp_map are required for computing correlation matrices")
    }
    LD_matrix_files <- unlist(strsplit(LD_map$LD_file[LD_map$region_id == region_id], split = ";"))
    stopifnot(all(file.exists(LD_matrix_files)))

    if (length(LD_matrix_files) > 1) {
      R_snp <- lapply(LD_matrix_files, load_LD, format = LD_format, LD_loader_fun = LD_loader_fun)
      R_snp <- suppressWarnings(as.matrix(bdiag(R_snp)))
    } else {
      R_snp <- load_LD(LD_matrix_files, format = LD_format, LD_loader_fun = LD_loader_fun)
    }

    # load SNP info of the region
    snpinfo <- snp_map[[region_id]]

    # compute correlation matrices
    sids <- region_data[[region_id]]$sid
    gids <- region_data[[region_id]]$gid
    res <- compute_region_cor(sids, gids, R_snp, snpinfo$id, weights)
    R_snp <- res$R_snp
    R_snp_gene <- res$R_snp_gene
    R_gene <- res$R_gene
    rm(res)
    # save correlation matrices
    if (isTRUE(save_cor && !is.null(cor_dir))) {
      saveRDS(R_snp_gene, file=R_sg_file)
      saveRDS(R_gene, file=R_g_file)
      saveRDS(R_snp, file=R_s_file)
    }
  }
  return(list("R_snp" = R_snp,
              "R_snp_gene" = R_snp_gene,
              "R_gene" = R_gene))
}


# @title Computes correlation matrices for a single region
#
# @param sids SNP IDs
#
# @param gids gene IDs
#
# @param R_snp LD (R) matrix
#
# @param LD_sids SNP IDs for the rows and columns of the LD matrix
#
# @param weights a list of weights for all the genes
#
# @return a list of correlation matrices (R_snp, R_snp_gene and R_gene)
#
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
      gene_corrs <- apply(gene_pairs, 2, function(x){
        t(wgtr[[x[1]]])%*%R_snp[ldr[[x[1]]], ldr[[x[2]]]]%*%wgtr[[x[2]]]/(
          sqrt(t(wgtr[[x[1]]])%*%R_snp[ldr[[x[1]]], ldr[[x[1]]]]%*%wgtr[[x[1]]]) *
            sqrt(t(wgtr[[x[2]]])%*%R_snp[ldr[[x[2]]], ldr[[x[2]]]]%*%wgtr[[x[2]]]))})
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
