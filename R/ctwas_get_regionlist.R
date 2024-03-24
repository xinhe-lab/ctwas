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
#' @param filter_z_ids If TRUE, only keep SNPs and genes in z_snp and z_gene.
#'
#' @param trim_by remove SNPs if the total number of SNPs exceeds limit, options: "random",
#' or "z" (trim SNPs with lower |z|) See parameter `maxSNP` for more information.
#'
#' @param seed seed for.random sampling
#'
#' @param minvar minimum number of SNPs (and genes) in a region
#'
#' @param mingene minimum number of genes in a region
#'
#' @param adjust_boundary identify cross-boundary genes, adjust regionlist and update weigh_list
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
                           weights = NULL,
                           thin = 1,
                           maxSNP = Inf,
                           filter_z_ids = TRUE,
                           trim_by = c("random", "z"),
                           seed = 99,
                           minvar = 1,
                           mingene = 0,
                           adjust_boundary = TRUE) {

  trim_by <- match.arg(trim_by)

  loginfo("No. regions: %s", nrow(region_info))

  if (thin > 1 | thin <= 0){
    stop("thin needs to be in (0,1]")
  }

  # get gene info from z_gene and weights
  gene_info <- get_gene_info(z_gene, weights, region_info)

  region_info <- region_info[order(region_info$chrom, region_info$start),]

  regionlist <- list()
  boundary_genes <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(boundary_genes) <- c("gene","chrom","region1","region2")

  # get region data for regions in each chromosome
  for (b in unique(region_info$chrom)){
    # select regions in the chromosome
    regioninfo <- region_info[region_info$chrom == b, ]
    # select genes in the chromosome
    geneinfo <- gene_info[gene_info$chrom == b, ]

    # get snp info in LD in the chromosome
    snpinfo <- read_LD_SNP_files(regioninfo$SNP_info) #ctwas
    if (unique(snpinfo$chrom) != b){
      stop("Input genotype file not split by chromosome or not in correct order")
    }

    # select SNPs
    snpinfo$keep <- rep(1, nrow(snpinfo))
    if (isTRUE(filter_z_ids)){
      # remove SNPs not in z_snp
      snpinfo$keep[!(snpinfo$id %in% z_snp$id)] <- 0
    }

    # downsampling for SNPs
    snpinfo$thin_tag <- rep(0, nrow(snpinfo))
    nkept <- round(nrow(snpinfo) * thin)
    set.seed(seed)
    snpinfo$thin_tag[sample(1:nrow(snpinfo), nkept)] <- 1

    # check geneinfo
    if (nrow(geneinfo)!=0){
      if (unique(geneinfo$chrom) != b){
        stop("Imputed expression not by chromosome or not in correct order")
      }

      # select genes
      geneinfo$keep <- 1
      if (isTRUE(filter_z_ids)){
        # remove genes not in z_gene
        geneinfo[!(geneinfo$id %in% z_gene$id), "keep"] <- 0
      }
    }

    # get regionlist and boundary_genes for regions in the chromosome, and update weights
    res <- assign_region_ids(regioninfo,geneinfo,snpinfo,weights,minvar,adjust_boundary)
    loginfo("No. regions with at least one SNP/gene for chr%s: %d", b, length(res$regionlist))
    regionlist <- c(regionlist, res$regionlist)
    weights <- res$weights
    boundary_genes <- rbind(boundary_genes,res$boundary_genes)
  }

  # Trim regions with SNPs more than maxSNP
  regionlist <- trim_regionlist(regionlist, z_snp, trim_by = trim_by, maxSNP = maxSNP, seed = seed)

  return(list(regionlist=regionlist,
              weights=weights,
              boundary_genes=boundary_genes))
}

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
