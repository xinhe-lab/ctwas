#' Get regionlist for each region in region_info
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
#' @param adjust_boundary_genes identify cross-boundary genes, adjust regionlist and update weighs
#'
#' @param ncore The number of cores used to parallelize susie over regions
#'
#' @param seed seed for random sampling
#'
#' @importFrom logging loginfo
#'
#' @return a list with regionlist, updated weights, and cross-bounary genes
#'
#' @export
#'
get_regionlist <- function(region_info,
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
  gene_info <- get_gene_info(weights, region_info)

  regionlist <- list()
  # get regionlist for each chromosome
  for (b in unique(region_info$chrom)){
    # select regions in the chromosome
    regioninfo <- region_info[region_info$chrom == b, ]

    # select genes in the chromosome
    geneinfo <- gene_info[gene_info$chrom == b, ]

    # get SNP info in LD in the chromosome
    snpinfo <- read_LD_SNP_files(regioninfo$SNP_info)

    # select SNPs
    snpinfo$keep <- rep(1, nrow(snpinfo))
    # remove SNPs not in z_snp
    snpinfo$keep[!(snpinfo$id %in% z_snp$id)] <- 0

    # downsampling for SNPs
    snpinfo$thin_tag <- rep(0, nrow(snpinfo))
    nkept <- round(nrow(snpinfo) * thin)
    set.seed(seed)
    snpinfo$thin_tag[sample(1:nrow(snpinfo), nkept)] <- 1

    if (nrow(geneinfo)!=0){
      # select genes
      geneinfo$keep <- 1
      # remove genes not in z_gene
      geneinfo[!(geneinfo$id %in% z_gene$id), "keep"] <- 0
    }

    # get regionlist for the chromosome
    regionlist_chr <- assign_region_ids(regioninfo,geneinfo,snpinfo)
    loginfo("No. regions with at least one SNP/gene for chr%s: %d", b, length(regionlist_chr))
    regionlist <- c(regionlist, regionlist_chr)
  }

  # adjust regionlist for boundary genes
  if (isTRUE(adjust_boundary_genes)){
    loginfo("Adjust regionlist for across boundary genes...")
    boundary_genes <- gene_info[gene_info$n_regions > 1, ]
    boundary_genes <- boundary_genes[with(boundary_genes, order(chrom, p0)), ]
    rownames(boundary_genes) <- NULL
    if (nrow(boundary_genes) > 0) {
      regionlist <- adjust_boundary_genes(boundary_genes, region_info, weights, regionlist)
    }
  }

  # trim regions with SNPs more than maxSNP
  regionlist <- trim_regionlist(regionlist, z_snp, trim_by = trim_by, maxSNP = maxSNP, seed = seed)

  # add z-scores to regionlist
  loginfo("Add z-scores to regionlist...")
  regionlist <- add_z_to_regionlist(regionlist, z_snp, z_gene, ncore = ncore)

  return(list(regionlist=regionlist, boundary_genes=boundary_genes))
}

#' Remove SNPs from regionlist if the total number of SNPs exceeds limit
trim_regionlist <- function(regionlist, z_snp, trim_by = c("random", "z"), maxSNP = Inf, seed = 99){

  trim_by <- match.arg(trim_by)

  if (maxSNP < Inf){
    if (trim_by == "z") {
      # trim SNPs with lower |z|
      for (region_tag in names(regionlist)){
        if (length(regionlist[[region_tag]][["sid"]]) > maxSNP){
          loginfo("Trim region %s with SNPs more than %s", region_tag, maxSNP)
          idx <- match(regionlist[[region_tag]][["sid"]], z_snp$id)
          z.abs <- abs(z_snp[idx, "z"])
          ifkeep <- rank(-z.abs) <= maxSNP
          regionlist[[region_tag]][["sid"]] <- regionlist[[region_tag]][["sid"]][ifkeep]
        }
      }
    } else {
      # randomly trim snps
      for (region_tag in names(regionlist)){
        if (length(regionlist[[region_tag]][["sid"]]) > maxSNP){
          loginfo("Trim region %s with SNPs more than %s", region_tag, maxSNP)
          n.snps <- length(regionlist[[region_tag]][["sid"]])
          ifkeep <- rep(FALSE, n.snps)
          set.seed(seed)
          ifkeep[sample.int(n.snps, size = maxSNP)] <- TRUE
          regionlist[[region_tag]][["sid"]] <-  regionlist[[region_tag]][["sid"]][ifkeep]
        }
      }
    }
  }

  return(regionlist)
}
