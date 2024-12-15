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
#' @param min_group_size Minimum number of genes for a group to be included.
#'
#' @param trim_by remove SNPs if the total number of SNPs exceeds limit,
#' options: "random", or "z" (trim SNPs with lower |z|).
#' See parameter `maxSNP` for more information.
#'
#' @param thin_by options for thinning SNPs,
#' "reference": thin reference SNPs,
#' "gwas": thin GWAS SNPs.
#'
#' @param adjust_boundary_genes If TRUE, identify cross-boundary genes, and
#' adjust region_data.
#'
#' @param ncore The number of cores used to parallelize susie over regions
#'
#' @param seed seed for random sampling
#'
#' @param logfile path to the log file, if NULL will print log info on screen.
#'
#' @return a list with region_data, updated weights, and cross-bounary genes
#'
#' @importFrom logging addHandler loginfo logwarn writeToFile
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
                                 min_group_size = 100,
                                 trim_by = c("random", "z"),
                                 thin_by = c("ref", "gwas"),
                                 adjust_boundary_genes = TRUE,
                                 ncore = 1,
                                 seed = 99,
                                 logfile = NULL) {

  if (!is.null(logfile)) {
    addHandler(writeToFile, file=logfile, level='DEBUG')
  }

  # check inputs
  trim_by <- match.arg(trim_by)
  thin_by <- match.arg(thin_by)

  if (anyNA(z_snp))
    stop("z_snp contains missing values!")

  if (anyNA(z_gene))
    stop("z_gene contains missing values!")

  if (!inherits(weights,"list"))
    stop("'weights' should be a list object.")

  if (any(sapply(weights, is.null)))
    stop("weights contain NULL, remove empty weights!")

  if (thin > 1 | thin <= 0)
    stop("thin needs to be in (0,1]")

  if (!inherits(snp_map,"list"))
    stop("'snp_map' should be a list object.")

  snp_info <- as.data.frame(rbindlist(snp_map, idcol = "region_id"))

  # begin assembling region_data
  region_ids <- region_info$region_id
  loginfo("Assembling region_data ...")
  loginfo("Number of regions in total: %d", length(region_ids))
  loginfo("thin = %s", thin)

  region_info <- region_info[order(region_info$chrom, region_info$start),]

  # get gene info from weights
  gene_info <- get_gene_info(weights)

  # filter groups with too few genes
  z_gene <- filter_z_gene_by_group_size(z_gene, min_group_size)

  # remove genes not in z_gene
  gene_info <- gene_info[gene_info$id %in% z_gene$id, , drop=FALSE]

  # flag SNPs not in z_snp
  snp_info$keep <- rep(1, nrow(snp_info))
  snp_info$keep[!(snp_info$id %in% z_snp$id)] <- 0

  region_data <- list()
  # get region_data for each chromosome
  for (b in unique(region_info$chrom)){

    # select regions in the chromosome
    region_info_chr <- region_info[region_info$chrom == b, ]

    # read SNPs in the chromosome
    snp_info_chr <- snp_info[snp_info$chrom == b, ]

    if (nrow(snp_info_chr) == 0) {
      stop(paste("No reference SNP info in chrom", b, "!"))
    }

    # select genes in the chromosome
    gene_info_chr <- gene_info[gene_info$chrom == b, ]

    # get region_data for the chromosome
    region_data_chr <- assign_region_data(region_info_chr,
                                          snp_info_chr,
                                          gene_info_chr,
                                          thin = thin,
                                          thin_by = thin_by,
                                          seed = seed,
                                          ncore = ncore)

    loginfo("Number of regions in chr%s: %d", b, length(region_data_chr))
    region_data <- c(region_data, region_data_chr)
  }

  # adjust region_data for boundary genes
  if (adjust_boundary_genes && nrow(region_info) > 1){
    gene_info <- get_gene_regions(gene_info, region_info)
    boundary_genes <- gene_info[gene_info$n_regions > 1, ]
    boundary_genes <- boundary_genes[with(boundary_genes, order(chrom, p0)), ]
    rownames(boundary_genes) <- NULL
    loginfo("Number of boundary genes: %d", nrow(boundary_genes))
    if (nrow(boundary_genes) > 0) {
      region_data <- adjust_boundary_genes(boundary_genes, weights, region_data, snp_map)
    }
  } else {
    boundary_genes <- NULL
  }

  # trim regions with SNPs more than maxSNP
  region_data <- trim_region_data(region_data, z_snp, trim_by = trim_by,
                                  maxSNP = maxSNP, seed = seed)

  # add z-scores to region_data
  loginfo("Adding region z-scores ...")
  region_data <- update_region_z(region_data, z_snp, z_gene, ncore = ncore)

  return(list("region_data" = region_data,
              "boundary_genes" = boundary_genes))
}

# Assign genes and SNPs in regions
assign_region_data <- function(region_info,
                               snp_info,
                               gene_info,
                               thin = 0.1,
                               thin_by = c("ref", "gwas"),
                               seed = 99,
                               ncore = 1) {

  thin_by <- match.arg(thin_by)

  # down-sampling SNPs
  if (thin < 1) {
    set.seed(seed)
    if (thin_by == "ref"){
      # thin reference SNPs
      thin_idx <- 1:nrow(snp_info)
    } else if (thin_by == "gwas"){
      # thin GWAS SNPs
      thin_idx <- which(snp_info$keep == 1)
    }
    snp_info$thin_tag <- 0
    nkept <- round(length(thin_idx) * thin)
    snp_info$thin_tag[sample(thin_idx, nkept)] <- 1
  } else {
    snp_info$thin_tag <- 1
  }

  # assign genes and SNPs to the region
  region_data <- mclapply_check(1:nrow(region_info), function(i){
    region_id <- region_info$region_id[i]
    region_chrom <- region_info$chrom[i]
    region_start <- region_info$start[i]
    region_stop <- region_info$stop[i]

    # assign genes to regions based gene p0 positions
    # for genes across region boundaries, assign to the first region, and adjust later
    gidx <- which(gene_info$chrom == region_chrom & gene_info$p0 >= region_start & gene_info$p0 < region_stop)

    sidx <- which(snp_info$chrom == region_chrom & snp_info$pos >= region_start & snp_info$pos < region_stop
                  & snp_info$keep == 1 & snp_info$thin_tag == 1)

    if (length(gidx) + length(sidx) < 1) {
      gid <- NULL
      sid <- NULL
      minpos <- NULL
      maxpos <- NULL
    } else {
      gid <- gene_info$id[gidx]
      sid <- snp_info$id[sidx]
      minpos <- min(c(gene_info$p0[gidx], snp_info$pos[sidx]))
      maxpos <- max(c(gene_info$p1[gidx], snp_info$pos[sidx]))
    }

    list("region_id" = region_id,
         "chrom" = region_chrom,
         "start" = region_start,
         "stop" = region_stop,
         "minpos" = minpos,
         "maxpos" = maxpos,
         "thin" = thin,
         "gid" = gid,
         "sid" = sid)
  }, mc.cores = ncore)

  names(region_data) <- region_info$region_id

  return(region_data)
}


# Trim SNPs from region_data if the total number of SNPs exceeds limit
#' @importFrom logging loginfo
trim_region_data <- function(region_data,
                             z_snp,
                             trim_by = c("random", "z"),
                             maxSNP = Inf,
                             seed = 99){

  trim_by <- match.arg(trim_by)

  if (!inherits(region_data,"list")){
    stop("'region_data' should be a list.")
  }

  if (maxSNP < Inf){
    if (trim_by == "z") {
      # trim SNPs with lower |z|
      for (region_id in names(region_data)){
        if (length(region_data[[region_id]][["sid"]]) > maxSNP){
          loginfo("Trim region %s with SNPs more than %s", region_id, maxSNP)
          idx <- match(region_data[[region_id]][["sid"]], z_snp$id)
          z.abs <- abs(z_snp[idx, "z"])
          keep_idx <- which(rank(-z.abs) <= maxSNP)
          region_data[[region_id]][["sid"]] <- region_data[[region_id]][["sid"]][keep_idx]
        }
      }
    } else {
      # randomly trim SNPs
      set.seed(seed)
      for (region_id in names(region_data)){
        if (length(region_data[[region_id]][["sid"]]) > maxSNP){
          loginfo("Trim region %s with SNPs more than %s", region_id, maxSNP)
          n.snps <- length(region_data[[region_id]][["sid"]])
          keep_idx <- sample.int(n.snps, size = maxSNP)
          region_data[[region_id]][["sid"]] <-  region_data[[region_id]][["sid"]][keep_idx]
        }
      }
    }
  }

  return(region_data)
}

#' @title Adds or updates z-scores in region_data based on z_snp and z_gene.
#' this will also update sid and gid based on z_snp and z_gene.
#
#' @param region_data a list of region gene IDs and SNP IDs and associated file names
#'
#' @param z_snp A data frame with columns: "id", "z", giving the z-scores for SNPs.
#'
#' @param z_gene A data frame with columns: "id", "z", giving the z-scores for genes.
#'
#' @param update options to update z-scores in region data.
#' "all": update all data (default),
#' "snps": updates SNP data only,
#' "genes": updates gene data only.
#'
#' @param ncore The number of cores used to parallelize susie over regions
#'
#' @importFrom parallel mclapply
#'
#' @export
update_region_z <- function(region_data,
                            z_snp,
                            z_gene,
                            update = c("all", "snps", "genes"),
                            ncore = 1){

  update <- match.arg(update)

  if (!inherits(region_data,"list")){
    stop("'region_data' should be a list.")
  }

  # combine z_snp, z_gene and group information
  if (update == "snps") {
    z_snp$type <- "SNP"
    z_snp$context <- "SNP"
    z_snp$group <- "SNP"
    z_df <- z_snp[, c("id", "z", "type", "context", "group")]
  } else if (update == "genes") {
    z_df <- z_gene[, c("id", "z", "type", "context", "group")]
  } else {
    z_df <- combine_z(z_snp, z_gene)
  }

  region_ids <- names(region_data)
  region_data2 <- mclapply_check(region_ids, function(region_id){
    regiondata <- region_data[[region_id]]
    gid <- regiondata[["gid"]]
    sid <- regiondata[["sid"]]
    region_z_df <- z_df[z_df$id %in% c(gid, sid), , drop=FALSE]

    # update z-scores for genes
    if (update == "genes" || update == "all") {
      gid <- intersect(gid, region_z_df$id)
      regiondata[["gid"]] <- gid
      region_z_gene <- region_z_df[match(gid, region_z_df$id), ]
      region_z_gene <- region_z_gene[complete.cases(region_z_gene),]
      regiondata[["z_gene"]] <- region_z_gene
      rownames(regiondata[["z_gene"]]) <- NULL
    }

    # update z-scores for snps
    if (update == "snps" || update == "all") {
      sid <- intersect(sid, region_z_df$id)
      regiondata[["sid"]] <- sid
      region_z_snp <- region_z_df[match(sid, region_z_df$id), ]
      region_z_snp <- region_z_snp[complete.cases(region_z_snp),]
      regiondata[["z_snp"]] <- region_z_snp
      rownames(regiondata[["z_snp"]]) <- NULL
    }

    regiondata[["types"]] <- c(unique(regiondata[["z_gene"]]$type), "SNP")
    regiondata[["contexts"]] <- c(unique(regiondata[["z_gene"]]$context), "SNP")
    regiondata[["groups"]] <- c(unique(regiondata[["z_gene"]]$group), "SNP")

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
  loginfo("Adjust region assignment for boundary genes")

  if (!inherits(weights,"list")){
    stop("'weights' should be a list.")
  }

  if (!inherits(region_data,"list")){
    stop("'region_data' should be a list.")
  }

  if (!inherits(snp_map,"list")){
    stop("'snp_map' should be a list.")
  }

  # assign boundary gene to the region with max weights
  for (i in 1:nrow(boundary_genes)){
    gname <- boundary_genes[i, "id"]
    region_ids <- unlist(strsplit(boundary_genes[i, "region_id"], split = ","))
    wgt <- weights[[gname]][["wgt"]]

    region_sum_wgt <- sapply(region_ids, function(region_id){
      snpinfo <- snp_map[[region_id]]
      sum(abs(wgt)[which(rownames(wgt) %in% snpinfo$id)])
    })
    selected_region_id <- region_ids[which.max(region_sum_wgt)]
    unselected_region_ids <- setdiff(region_ids, selected_region_id)
    region_data[[selected_region_id]][["gid"]] <- union(region_data[[selected_region_id]][["gid"]], gname)
    for(unselected_region_id in unselected_region_ids){
      region_data[[unselected_region_id]][["gid"]] <- setdiff(region_data[[unselected_region_id]][["gid"]], gname)
    }
  }

  return(region_data)
}

#' @title Expands region_data with all SNPs
#'
#' @param region_data a list of region gene IDs and SNP IDs and associated file names
#'
#' @param snp_map a list of data frames with SNP-to-region map for the reference.
#'
#' @param z_snp A data frame with columns: "id", "z", giving the z-scores for SNPs.
#'
#' @param maxSNP Inf or integer. Maximum number of SNPs in a region. Default is
#' Inf, no limit. This can be useful if there are many SNPs in a region and you don't
#' have enough memory to run the program.
#'
#' @param ncore The number of cores used to parallelize susie over regions
#'
#' @return updated region_data with all SNPs
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
                               maxSNP = Inf,
                               ncore = 1) {
  # check inputs
  if (!inherits(region_data,"list"))
    stop("'region_data' should be a list.")

  if (!inherits(snp_map,"list"))
    stop("'snp_map' should be a list object.")

  if (anyNA(z_snp))
    stop("z_snp contains missing values!")

  # update SNP IDs for each region
  thin <- sapply(region_data, "[[", "thin")

  if (length(which(thin < 1)) > 0) {
    loginfo("Expand %d regions with all SNPs", length(which(thin < 1)))
    region_ids <- names(region_data)
    region_data <- mclapply_check(region_ids, function(region_id){
      # add z-scores and types of the region to the region_data
      regiondata <- region_data[[region_id]]
      if (regiondata[["thin"]] < 1){
        # load all SNPs in the region
        snpinfo <- snp_map[[region_id]]
        # remove SNPs not in z_snp
        snpinfo <- snpinfo[snpinfo$id %in% z_snp$id, ,drop=FALSE]
        # update minpos and maxpos in the region
        regiondata[["minpos"]] <- min(c(regiondata[["minpos"]], snpinfo$pos))
        regiondata[["maxpos"]] <- max(c(regiondata[["maxpos"]], snpinfo$pos))
        # set thin to 1 after expanding SNPs
        regiondata[["thin"]] <- 1
        # update SNPs in the region
        regiondata[["sid"]] <- snpinfo$id
      }
      regiondata
    }, mc.cores = ncore)

    names(region_data) <- region_ids

    # Trim regions with SNPs more than maxSNP
    region_data <- trim_region_data(region_data, z_snp, trim_by = "z", maxSNP = maxSNP)

    # add z-scores to region_data
    loginfo("Updating region z-scores ...")
    region_data <- update_region_z(region_data, z_snp, update = "snps", ncore = ncore)
  }

  return(region_data)

}

# extract data for a region from region_data
extract_region_data <- function(region_data,
                                region_id,
                                groups,
                                snps_only = FALSE){

  if (!inherits(region_data,"list")){
    stop("'region_data' should be a list.")
  }

  regiondata <- region_data[[region_id]]

  stopifnot(all.equal(regiondata$sid, regiondata$z_snp$id))
  stopifnot(all.equal(regiondata$gid, regiondata$z_gene$id))

  if (snps_only) {
    regiondata$gid <- NULL
    regiondata$z_gene <- NULL
    regiondata$z <- regiondata$z_snp$z
    regiondata$gs_type <- regiondata$z_snp$type
    regiondata$gs_context <- regiondata$z_snp$context
    regiondata$gs_group <- regiondata$z_snp$group
    regiondata$g_type <- NULL
    regiondata$g_context <- NULL
    regiondata$g_group <- NULL
    regiondata$types <- "SNP"
    regiondata$contexts <- "SNP"
    regiondata$groups <- "SNP"
  } else {
    if (!missing(groups)) {
      regiondata$z_gene <- regiondata$z_gene[regiondata$z_gene$group %in% groups,,drop=FALSE]
      regiondata$z_snp <- regiondata$z_snp[regiondata$z_snp$group %in% groups,,drop=FALSE]
    }
    regiondata$gid <- regiondata$z_gene$id
    regiondata$sid <- regiondata$z_snp$id
    regiondata$z <- c(regiondata$z_gene$z, regiondata$z_snp$z)
    regiondata$gs_type <- c(regiondata$z_gene$type, regiondata$z_snp$type)
    regiondata$gs_context <- c(regiondata$z_gene$context, regiondata$z_snp$context)
    regiondata$gs_group <- c(regiondata$z_gene$group, regiondata$z_snp$group)
    regiondata$g_type <- regiondata$z_gene$type
    regiondata$g_context <- regiondata$z_gene$context
    regiondata$g_group <- regiondata$z_gene$group
    regiondata$types <- unique(regiondata$gs_type)
    regiondata$contexts <- unique(regiondata$gs_context)
    regiondata$groups <- unique(regiondata$gs_group)
  }

  return(regiondata)
}

# filter region_data in selected groups
filter_region_data <- function(region_data, groups) {
  region_ids <- names(region_data)
  region_data2 <- lapply(region_ids, function(region_id){
    regiondata <- region_data[[region_id]]
    regiondata$z_gene <- regiondata$z_gene[regiondata$z_gene$group %in% groups,,drop=FALSE]
    regiondata$z_snp <- regiondata$z_snp[regiondata$z_snp$group %in% groups,,drop=FALSE]
    regiondata$gid <- regiondata$z_gene$id
    regiondata$sid <- regiondata$z_snp$id
    regiondata$z <- c(regiondata$z_gene$z, regiondata$z_snp$z)
    regiondata$gs_type <- c(regiondata$z_gene$type, regiondata$z_snp$type)
    regiondata$gs_context <- c(regiondata$z_gene$context, regiondata$z_snp$context)
    regiondata$gs_group <- c(regiondata$z_gene$group, regiondata$z_snp$group)
    regiondata$g_type <- regiondata$z_gene$type
    regiondata$g_context <- regiondata$z_gene$context
    regiondata$g_group <- regiondata$z_gene$group
    regiondata$types <- unique(regiondata$gs_type)
    regiondata$contexts <- unique(regiondata$gs_context)
    regiondata$groups <- unique(regiondata$gs_group)
    regiondata
  })
  names(region_data2) <- region_ids

  return(region_data2)
}

# get group sizes from region data
get_group_size_from_region_data <- function(region_data) {
  z_gene_all <- lapply(region_data, "[[", "z_gene")
  z_gene_all <- as.data.frame(rbindlist(z_gene_all, idcol = "region_id"))

  z_snp_all <- lapply(region_data, "[[", "z_snp")
  z_snp_all <- as.data.frame(rbindlist(z_snp_all, idcol = "region_id"))

  z_df <- combine_z(z_snp_all, z_gene_all)
  group_size <- table(z_df$group)
  group_names <- names(group_size)
  group_size <- as.numeric(group_size)
  names(group_size) <- group_names
  group_names <- c(setdiff(group_names, "SNP"), "SNP")
  group_size <- group_size[group_names]

  return(group_size)
}
