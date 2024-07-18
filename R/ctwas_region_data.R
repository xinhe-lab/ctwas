#' @title Assembles data for all the regions
#'
#' @param region_info a data frame of region definitions.
#'
#' @param z_snp A data frame with columns: "id", "z", giving the z-scores for SNPs.
#'
#' @param z_gene A data frame with columns: "id", "z", giving the z-scores for genes.
#'
#' @param weights a list of preprocessed weights.
#'
#' @param snp_map a list of data frames with SNP-to-region map for the reference.
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
#' @param trim_by remove SNPs if the total number of SNPs exceeds limit,
#' options: "random", or "z" (trim SNPs with lower |z|).
#' See parameter `maxSNP` for more information.
#'
#' @param adjust_boundary_genes identify cross-boundary genes, adjust region_data
#'
#' @param thin_gwas_snps TRUE/FALSE, if TRUE, only apply thin to GWAS SNPs.
#' Otherwise, apply thins to all SNPs.
#'
#' @param add_z TRUE/FALSE, if TRUE, add z-scores to the region_data
#'
#' @param ncore The number of cores used to parallelize susie over regions
#'
#' @param seed seed for random sampling
#'
#' @param logfile path to the log file, if NULL will print log info on screen.
#'
#' @return a list with region_data, updated weights, and cross-bounary genes
#'
#' @importFrom logging loginfo
#' @importFrom data.table rbindlist
#'
#' @export
#'
assemble_region_data <- function(region_info,
                                 z_snp,
                                 z_gene,
                                 weights,
                                 snp_map,
                                 thin = 0.1,
                                 maxSNP = Inf,
                                 trim_by = c("random", "z"),
                                 adjust_boundary_genes = TRUE,
                                 thin_gwas_snps = TRUE,
                                 add_z = TRUE,
                                 ncore = 1,
                                 seed = 99,
                                 logfile = NULL) {

  # check inputs
  trim_by <- match.arg(trim_by)

  if (anyNA(z_snp)){
    stop("z_snp contains missing values!")
  }

  if (anyNA(z_gene)){
    stop("z_gene contains missing values!")
  }

  if (!inherits(weights,"list")){
    stop("'weights' should be a list object.")
  }

  if (any(sapply(weights, is.null))) {
    stop("weights contain NULL, remove empty weights!")
  }

  if (thin > 1 | thin <= 0){
    stop("thin needs to be in (0,1]")
  }

  snp_info <- as.data.frame(rbindlist(snp_map, idcol = "region_id"))

  # begin assembling region_data
  region_ids <- region_info$region_id
  loginfo("Assembling region_data ...")
  loginfo("Number of regions in total: %d", length(region_ids))
  loginfo("thin = %s", thin)

  region_info <- region_info[order(region_info$chrom, region_info$start),]

  # get gene info from weights
  gene_info <- get_gene_info(weights)

  region_data <- list()
  # get region_data for each chromosome
  for (b in unique(region_info$chrom)){

    # select regions in the chromosome
    regioninfo <- region_info[region_info$chrom == b, ]

    # select genes in the chromosome
    geneinfo <- gene_info[gene_info$chrom == b, ]

    # read SNP info in the chromosome
    snpinfo <- snp_info[snp_info$chrom == b, ]

    # select SNPs
    snpinfo$keep <- rep(1, nrow(snpinfo))
    # remove SNPs not in z_snp
    snpinfo$keep[!(snpinfo$id %in% z_snp$id)] <- 0

    if (nrow(geneinfo)!=0){
      # select genes
      geneinfo$keep <- 1
      # remove genes not in z_gene
      geneinfo[!(geneinfo$id %in% z_gene$id), "keep"] <- 0
    }

    # get region_data for the chromosome
    region_data_chr <- assign_region_ids(regioninfo,
                                         geneinfo,
                                         snpinfo,
                                         thin = thin,
                                         thin_gwas_snps = thin_gwas_snps,
                                         seed = seed)
    loginfo("Number of regions in chr%s: %d", b, length(region_data_chr))
    region_data <- c(region_data, region_data_chr)
  }

  # adjust region_data for boundary genes
  if (isTRUE(adjust_boundary_genes) && nrow(region_info) > 1){
    gene_info <- get_gene_regions(gene_info, region_info)
    boundary_genes <- gene_info[gene_info$n_regions > 1, ]
    boundary_genes <- boundary_genes[with(boundary_genes, order(chrom, p0)), ]
    rownames(boundary_genes) <- NULL
    loginfo("Number of boundary genes: %d", nrow(boundary_genes))
    if (nrow(boundary_genes) > 0) {
      region_data <- adjust_boundary_genes(boundary_genes, weights, region_data, snp_map)
    }
  }else{
    boundary_genes <- NULL
  }

  # trim regions with SNPs more than maxSNP
  region_data <- trim_region_data(region_data, z_snp, trim_by = trim_by, maxSNP = maxSNP, seed = seed)

  # add z-scores to region_data
  if (add_z) {
    region_data <- add_z_to_region_data(region_data, z_snp, z_gene, ncore = ncore)
  }

  return(list(region_data=region_data, boundary_genes=boundary_genes))
}

# Assign gene and SNP IDs for regions in the regioninfo
assign_region_ids <- function(regioninfo,
                              geneinfo,
                              snpinfo,
                              thin = 0.1,
                              thin_gwas_snps = TRUE,
                              seed = 99) {

  # downsampling for SNPs if thin < 1
  if (thin < 1) {
    set.seed(seed)
    if (isTRUE(thin_gwas_snps)){
      # only thin GWAS snps with keep label = 1
      thin_idx <- which(snpinfo$keep == 1)
    } else {
      thin_idx <- 1:nrow(snpinfo)
    }
    snpinfo$thin_tag <- rep(0, nrow(snpinfo))
    nkept <- round(length(thin_idx) * thin)
    snpinfo$thin_tag[sample(thin_idx, nkept)] <- 1
  } else {
    snpinfo$thin_tag <- 1
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
    gidx <- which(geneinfo$chrom == region_chrom & geneinfo$p0 >= region_start & geneinfo$p0 < region_stop
                  & geneinfo$keep == 1)

    sidx <- which(snpinfo$chrom == region_chrom & snpinfo$pos >= region_start & snpinfo$pos < region_stop
                  & snpinfo$keep == 1 & snpinfo$thin_tag == 1)

    if (length(gidx) + length(sidx) < 1) {
      gid <- NULL
      sid <- NULL
      minpos <- NULL
      maxpos <- NULL
    }else{
      gid <- geneinfo$id[gidx]
      sid <- snpinfo$id[sidx]
      minpos <- min(c(geneinfo$p0[gidx], snpinfo$pos[sidx]))
      maxpos <- max(c(geneinfo$p1[gidx], snpinfo$pos[sidx]))
    }

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


# Trim SNPs from region_data if the total number of SNPs exceeds limit
#' @importFrom logging loginfo
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
      # randomly trim SNPs
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

# add z-scores from z_snp and z_gene to region_data
#' @importFrom logging loginfo
#' @importFrom parallel mclapply
add_z_to_region_data <- function(region_data,
                                 z_snp,
                                 z_gene,
                                 ncore = 1){
  loginfo("Adding z-scores to region_data ...")

  # Combine z-scores from z_snp and z_gene
  zdf <- combine_z(z_snp, z_gene)

  region_ids <- names(region_data)
  region_data2 <- mclapply_check(region_ids, function(region_id){
    # add z-scores and types of the region to the region_data
    regiondata <- region_data[[region_id]]
    gid <- regiondata[["gid"]]
    sid <- regiondata[["sid"]]
    region_zdf <- zdf[zdf$id %in% c(gid, sid), ]
    region_zdf <- region_zdf[match(c(gid, sid), region_zdf$id),]
    regiondata[["z"]] <- region_zdf$z
    regiondata[["gs_group"]] <- region_zdf$group
    regiondata
  }, mc.cores = ncore)

  names(region_data2) <- region_ids

  return(region_data2)
}

# Adjust region_data for boundary genes
#' @importFrom logging loginfo
#' @importFrom data.table rbindlist
adjust_boundary_genes <- function(boundary_genes,
                                  weights,
                                  region_data,
                                  snp_map){
  loginfo("Adjusting for boundary genes ...")

  if (!inherits(weights,"list")){
    stop("'weights' should be a list.")
  }

  for (i in 1:nrow(boundary_genes)){
    gname <- boundary_genes[i, "id"]
    region_ids <- unlist(strsplit(boundary_genes[i, "region_id"], split = ";"))
    wgt <- weights[[gname]][["wgt"]]

    region_r2 <- sapply(region_ids, function(region_id){
      snpinfo <- snp_map[[region_id]]
      sum(wgt[which(rownames(wgt) %in% snpinfo$id)]^2)
    })

    # assign boundary gene to the region with max r2, and remove it from other regions
    selected_region_id <- region_ids[which.max(region_r2)]
    unselected_region_ids <- setdiff(region_ids, selected_region_id)
    region_data[[selected_region_id]][["gid"]] <- union(region_data[[selected_region_id]][["gid"]], gname)
    for(unselected_region_id in unselected_region_ids){
      region_data[[unselected_region_id]][["gid"]] <- setdiff(region_data[[unselected_region_id]][["gid"]], gname)
    }
  }

  return(region_data)
}

#' @title Expands region_data with full SNPs
#'
#' @param region_data a list of region gene IDs and SNP IDs and associated file names
#'
#' @param snp_map a list of data frames with SNP-to-region map for the reference.
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
#' @param add_z TRUE/FALSE, if TRUE, add z-scores to the region_data
#'
#' @param ncore The number of cores used to parallelize susie over regions
#'
#' @param seed seed for random sampling
#'
#' @return updated region_data with full SNPs
#'
#' @importFrom logging loginfo
#' @importFrom data.table rbindlist
#' @importFrom parallel mclapply
#'
#' @export
#'
expand_region_data <- function(region_data,
                               snp_map,
                               z_snp,
                               z_gene,
                               trim_by = c("z", "random"),
                               maxSNP = Inf,
                               add_z = TRUE,
                               ncore = 1,
                               seed = 99) {
  # check arguments
  trim_by <- match.arg(trim_by)

  if (anyNA(z_snp)){
    stop("z_snp contains missing values!")
  }

  if (anyNA(z_gene)){
    stop("z_gene contains missing values!")
  }

  # update SNP IDs for each region
  thin <- sapply(region_data, "[[", "thin")
  loginfo("Expanding %d regions with full SNPs ...", length(which(thin < 1)))

  region_ids <- names(region_data)
  region_data <- mclapply_check(region_ids, function(region_id){
    # add z-scores and types of the region to the region_data
    regiondata <- region_data[[region_id]]
    if (regiondata[["thin"]] < 1){

      # load all SNPs in the region
      snpinfo <- snp_map[[region_id]]

      # update sid in the region
      snpinfo$keep <- rep(1, nrow(snpinfo))
      # remove SNPs not in z_snp
      snpinfo$keep[!(snpinfo$id %in% z_snp$id)] <- 0
      sid <- snpinfo$id[snpinfo$keep == 1]
      sidx <- match(sid, snpinfo$id)
      regiondata[["sid"]] <- sid
      # update minpos and maxpos in the region
      regiondata[["minpos"]] <- min(c(regiondata[["minpos"]], snpinfo$pos[sidx]))
      regiondata[["maxpos"]] <- max(c(regiondata[["maxpos"]], snpinfo$pos[sidx]))
      # set thin to 1 after expanding SNPs
      regiondata[["thin"]] <- 1
    }
    regiondata
  }, mc.cores = ncore)

  names(region_data) <- region_ids

  # Trim regions with SNPs more than maxSNP
  region_data <- trim_region_data(region_data, z_snp, trim_by = trim_by, maxSNP = maxSNP, seed = seed)

  # add z-scores to region_data
  if (add_z) {
    region_data <- add_z_to_region_data(region_data, z_snp, z_gene, ncore = ncore)
  }

  return(region_data)

}

