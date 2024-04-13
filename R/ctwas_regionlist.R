#' Assemble data for all the regions
#'
#' @param region_info a data frame of region definition and associated file names
#'
#' @param z_snp A data frame with columns: "id", "z", giving the z-scores for SNPs.
#'
#' @param z_gene A data frame with columns: "id", "z", giving the z-scores for genes.
#'
#' @param weights a list of weights
#'
#' @param gene_info a data frame of gene information obtained from \code{get_gene_info}
#'
#' @param thin The proportion of SNPs to be used for the parameter estimation and
#' initial screening region steps.
#' Smaller \code{thin} parameters reduce runtime at the expense of accuracy.
#' The fine mapping step is rerun using full SNPs for regions with strong gene signals.
#'
#' @param maxSNP Inf or integer. Maximum number of SNPs in a region. Default is
#' Inf, no limit. This can be useful if there are many SNPs in a region and you don't
#' have enough memory to run the program. This applies to the last rerun step
#' (using full SNPs and rerun susie for regions with strong gene signals) only.
#'
#' @param trim_by remove SNPs if the total number of SNPs exceeds limit, options: "random",
#' or "z" (trim SNPs with lower |z|) See parameter `maxSNP` for more information.
#'
#' @param minvar minimum number of SNPs (and genes) in a region
#'
#' @param mingene minimum number of genes in a region
#'
#' @param adjust_boundary_genes identify cross-boundary genes, adjust region_data and update weighs
#'
#' @param ncore The number of cores used to parallelize susie over regions
#'
#' @param seed seed for random sampling
#'
#' @importFrom logging loginfo
#'
#' @return a list with region_data, updated weights, and cross-bounary genes
#'
#' @export
#'
assemble_region_data <- function(region_info,
                                 z_snp,
                                 z_gene,
                                 weights,
                                 thin = 1,
                                 maxSNP = Inf,
                                 trim_by = c("random", "z"),
                                 minvar = 1,
                                 mingene = 0,
                                 adjust_boundary_genes = TRUE,
                                 ncore = 1,
                                 seed = 99) {

  trim_by <- match.arg(trim_by)

  loginfo("No. regions in total: %d", nrow(region_info))

  if (thin > 1 | thin <= 0){
    stop("thin needs to be in (0,1]")
  }

  region_info <- region_info[order(region_info$chrom, region_info$start),]

  # get gene info from weights
  loginfo("Get gene info from weights")
  gene_info <- get_gene_info(weights)

  region_data <- list()
  # get region_data for each chromosome
  for (b in unique(region_info$chrom)){

    # select regions in the chromosome
    regioninfo <- region_info[region_info$chrom == b, ]

    # select genes in the chromosome
    geneinfo <- gene_info[gene_info$chrom == b, ]

    # get SNP info in LD in the chromosome
    snpinfo <- read_LD_SNP_files(regioninfo$SNP_info)

    # remove SNPs not in z_snp
    snpinfo <- snpinfo[snpinfo$id %in% z_snp$id, , drop=FALSE]

    # remove genes not in z_gene
    geneinfo <- geneinfo[geneinfo$id %in% z_gene$id, , drop=FALSE]

    # get region_data for the chromosome
    region_data_chr <- assign_region_ids(regioninfo, geneinfo, snpinfo,
                                         thin = thin, seed = seed,
                                         minvar = minvar, mingene = mingene)
    loginfo("No. regions in chr%s: %d", b, length(region_data_chr))
    region_data <- c(region_data, region_data_chr)
  }

  # adjust region_data for boundary genes
  if (isTRUE(adjust_boundary_genes)){
    loginfo("Adjust region_data for across boundary genes...")
    gene_info <- get_gene_regions(gene_info, region_info)
    boundary_genes <- gene_info[gene_info$n_regions > 1, ]
    boundary_genes <- boundary_genes[with(boundary_genes, order(chrom, p0)), ]
    rownames(boundary_genes) <- NULL
    loginfo("No. boundary genes: %d", nrow(boundary_genes))
    if (nrow(boundary_genes) > 0) {
      region_data <- adjust_boundary_genes(boundary_genes, region_info, weights, region_data)
    }
  }

  # trim regions with SNPs more than maxSNP
  region_data <- trim_region_data(region_data, z_snp, trim_by = trim_by, maxSNP = maxSNP, seed = seed)

  # add z-scores to region_data
  loginfo("Add z-scores to region_data...")
  region_data <- add_z_to_region_data(region_data, z_snp, z_gene, ncore = ncore)

  return(list(region_data=region_data, boundary_genes=boundary_genes))
}

#' assign gene and SNP IDs for regions in the regioninfo
assign_region_ids <- function(regioninfo,
                              geneinfo,
                              snpinfo,
                              thin = 1,
                              minvar = 1,
                              mingene = 0,
                              seed = 99) {

  # downsampling for SNPs
  if (thin < 1) {
    set.seed(seed)
    n_kept <- round(nrow(snpinfo) * thin)
    idx_kept <- sample(1:nrow(snpinfo), n_kept)
    snpinfo <- snpinfo[idx_kept, , drop = FALSE]
  }

  region_data <- list()

  # assign genes and SNPs to the region
  for (i in 1:nrow(regioninfo)){

    region_id <- regioninfo$region_id[i]
    region_chrom <- regioninfo$chrom[i]
    region_start <- regioninfo$start[i]
    region_stop <- regioninfo$stop[i]

    # assign genes to regions based gene p0 positions
    # for genes across region boundaries, assign to the first region, and adjust later
    gidx <- which(geneinfo$chrom == region_chrom & geneinfo$p0 >= region_start & geneinfo$p0 < region_stop)

    sidx <- which(snpinfo$chrom == region_chrom & snpinfo$pos >= region_start & snpinfo$pos < region_stop)

    if (length(gidx) + length(sidx) < minvar) {next}

    if (length(gidx) < mingene) {next}

    gid <- geneinfo$id[gidx]
    sid <- snpinfo$id[sidx]

    minpos <- min(c(geneinfo$p0[gidx], snpinfo$pos[sidx]))
    maxpos <- max(c(geneinfo$p1[gidx], snpinfo$pos[sidx]))

    region_data[[region_id]] <- list("region_id" = region_id,
                                     "chrom" = region_chrom,
                                     "start" = region_start,
                                     "stop" = region_stop,
                                     "minpos" = minpos,
                                     "maxpos" = maxpos,
                                     "gid" = gid,
                                     "sid" = sid,
                                     "thin" = thin)
  }

  return(region_data)

}


#' Remove SNPs from region_data if the total number of SNPs exceeds limit
trim_region_data <- function(region_data, z_snp, trim_by = c("random", "z"), maxSNP = Inf, seed = 99){

  trim_by <- match.arg(trim_by)

  if (maxSNP < Inf){
    if (trim_by == "z") {
      # trim SNPs with lower |z|
      for (region_id in names(region_data)){
        if (length(region_data[[region_id]][["sid"]]) > maxSNP){
          loginfo("Trim region %s with SNPs more than %s", region_id, maxSNP)
          idx <- match(region_data[[region_id]][["sid"]], z_snp$id)
          z.abs <- abs(z_snp[idx, "z"])
          ifkeep <- rank(-z.abs) <= maxSNP
          region_data[[region_id]][["sid"]] <- region_data[[region_id]][["sid"]][ifkeep]
        }
      }
    } else {
      # randomly trim snps
      for (region_id in names(region_data)){
        if (length(region_data[[region_id]][["sid"]]) > maxSNP){
          loginfo("Trim region %s with SNPs more than %s", region_id, maxSNP)
          n.snps <- length(region_data[[region_id]][["sid"]])
          ifkeep <- rep(FALSE, n.snps)
          set.seed(seed)
          ifkeep[sample.int(n.snps, size = maxSNP)] <- TRUE
          region_data[[region_id]][["sid"]] <-  region_data[[region_id]][["sid"]][ifkeep]
        }
      }
    }
  }

  return(region_data)
}

#' add z-scores from z_snp and z_gene to region_data
add_z_to_region_data <- function(region_data,
                                 z_snp,
                                 z_gene,
                                 ncore = 1){

  # Combine z-scores from z_snp and z_gene
  zdf <- combine_z(z_snp, z_gene)

  cl <- parallel::makeCluster(ncore, outfile = "")
  doParallel::registerDoParallel(cl)
  corelist <- region2core(region_data, ncore)

  region_ids <- names(region_data)
  region_data2 <- foreach (core = 1:length(corelist), .combine = "c") %dopar% {
    region_data2.core <- list()
    region_ids.core <- corelist[[core]]
    for (region_id in region_ids.core) {
      # add z-scores and types of the region to the region_data
      region_data2.core[[region_id]] <- region_data[[region_id]]
      sid <- region_data[[region_id]][["sid"]]
      gid <- region_data[[region_id]][["gid"]]
      g_idx <- match(gid, zdf$id)
      s_idx <- match(sid, zdf$id)
      gs_idx <- c(g_idx, s_idx)
      # region_data2.core[[region_id]][["zdf"]] <- zdf[gs_idx,]
      region_data2.core[[region_id]][["z"]] <- zdf$z[gs_idx]
      region_data2.core[[region_id]][["gs_type"]] <- zdf$type[gs_idx]
      region_data2.core[[region_id]][["gs_context"]] <- zdf$context[gs_idx]
      region_data2.core[[region_id]][["gs_group"]] <- zdf$group[gs_idx]
      region_data2.core[[region_id]][["g_type"]] <- zdf$type[g_idx]
      region_data2.core[[region_id]][["g_context"]] <- zdf$context[g_idx]
      region_data2.core[[region_id]][["g_group"]] <- zdf$group[g_idx]
    }
    region_data2.core
  }
  parallel::stopCluster(cl)

  return(region_data2)
}


#' adjust region_data for boundary genes
adjust_boundary_genes <- function(boundary_genes, region_info, weights, region_data){

  for (i in 1:nrow(boundary_genes)){
    gname <- boundary_genes[i, "id"]
    region_ids <- unlist(strsplit(boundary_genes[i, "region_id"], split = ";"))
    wgt <- weights[[gname]][["wgt"]]

    region_r2 <- sapply(region_ids, function(x){
      ld_snpinfo <- read_LD_SNP_files(region_info[region_info$region_id == x, "SNP_info"])
      sum(wgt[which(rownames(wgt) %in% ld_snpinfo$id)]^2)})

    # assign boundary gene to the region with max r2, and remove it from other regions
    selected_region_id <- region_ids[which.max(region_r2)]
    unselected_region_ids <- setdiff(region_ids, selected_region_id)
    region_data[[selected_region_id]][["gid"]] <- unique(c(region_data[[selected_region_id]][["gid"]],gname))
    for(unselected_region_id in unselected_region_ids){
      region_data[[unselected_region_id]][["gid"]] <- region_data[[unselected_region_id]][["gid"]][region_data[[unselected_region_id]][["gid"]]!=gname]
    }
  }

  return(region_data)
}


#' Expand region_data with full SNPs
#'
#' @param region_data a list of region gene IDs and SNP IDs and associated file names
#'
#' @param region_info a data frame of region definition and associated LD file names
#'
#' @param z_snp A data frame with columns: "id", "z", giving the z-scores for SNPs.
#'
#' @param z_gene A data frame with columns: "id", "z", giving the z-scores for genes.
#'
#' @param trim_by remove SNPs if the total number of SNPs exceeds limit, options: "random",
#' or "z" (trim SNPs with lower |z|) See parameter `maxSNP` for more information.
#'
#' @param maxSNP Inf or integer. Maximum number of SNPs in a region. Default is
#' Inf, no limit. This can be useful if there are many SNPs in a region and you don't
#' have enough memory to run the program.
#'
#' @param ncore The number of cores used to parallelize susie over regions
#'
#' @param seed seed for random sampling
#'
#' @importFrom logging loginfo
#'
#' @return updated region_data with full SNPs
#'
#' @export
#'
expand_region_data <- function(region_data,
                               region_info,
                               z_snp,
                               z_gene,
                               trim_by = c("z", "random"),
                               maxSNP = Inf,
                               ncore = 1,
                               seed = 99) {

  trim_by <- match.arg(trim_by)

  region_ids <- names(region_data)
  loginfo("Expand region_data for %d regions with full SNPs", length(region_ids))

  pb <- txtProgressBar(min = 0, max = length(region_ids), initial = 0, style = 3)

  # update SNP IDs for each region
  for (i in 1:length(region_ids)){
    region_id <- region_ids[i]

    if (region_data[[region_id]][["thin"]] == 1){
      next
    }

    # load all SNPs in the region
    snpinfo <- read_LD_SNP_files(region_info[region_info$region_id == region_id, "SNP_info"])
    # remove SNPs not in z_snp
    snpinfo <- snpinfo[snpinfo$id %in% z_snp$id, , drop=FALSE]
    region_data[[region_id]][["sid"]] <- snpinfo$id

    # update minpos and maxpos in the region
    region_data[[region_id]][["minpos"]] <- min(c(region_data[[region_id]][["minpos"]], snpinfo$pos))
    region_data[[region_id]][["maxpos"]] <- max(c(region_data[[region_id]][["maxpos"]], snpinfo$pos))

    setTxtProgressBar(pb, i)
  }

  close(pb)

  # Trim regions with SNPs more than maxSNP
  region_data <- trim_region_data(region_data, z_snp, trim_by = trim_by, maxSNP = maxSNP, seed = seed)

  loginfo("Add z-scores to region_data ...")
  region_data <- add_z_to_region_data(region_data, z_snp, z_gene, ncore = ncore)

  return(region_data)

}

