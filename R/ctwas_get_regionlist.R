#' Get regionlist for each region in region_info
#'
#' @param region_info a data frame of region definition and associated file names
#'
#' @param gene_info a data frame of gene information obtained from \code{compute_gene_z}
#'
#' @param weight_list a list of weights
#'
#' @param select Default is NULL, all variants will be selected. Or a vector of variant IDs,
#' or a data frame with columns id and z (id is for gene or SNP id, z is for z scores).
#' z will be used for remove SNPs if the total number of SNPs exceeds limit. See
#' parameter `maxSNP` for more information.
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
#' @param minvar minimum number of SNPs (and genes) in a region
#'
#' @param adjust_boundary identify cross-boundary genes, adjust regionlist and update weigh_list
#'
#' @importFrom logging loginfo
#'
#' @return a list with regionlist, updated weight_list, and cross-bounary genes
#'
#' @export
#'
get_regionlist <- function(region_info,
                           gene_info,
                           weight_list,
                           select = NULL,
                           thin = 1,
                           maxSNP = Inf,
                           minvar = 1,
                           adjust_boundary = TRUE) {

  loginfo("No. LD regions: %s", nrow(region_info))

  if (thin > 1 | thin <= 0){
    stop("thin needs to be in (0,1]")
  }

  if (is.null(dim(select))) {
    selectid <- select
  } else {
    selectid <- select$id
  }

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
    if (isTRUE(unique(snpinfo$chrom) != b)){
      stop("Input genotype file not split by chromosome or not in correct order")
    }
    # select snps
    snpinfo$keep <- rep(1, nrow(snpinfo))
    if (!is.null(selectid)){
      snpinfo$keep[!(snpinfo$id %in% selectid)] <- 0
    }
    # downsampling for SNPs
    snpinfo$thin_tag <- rep(0, nrow(snpinfo))
    nkept <- round(nrow(snpinfo) * thin)
    set.seed(99)
    snpinfo$thin_tag[sample(1:nrow(snpinfo), nkept)] <- 1

    # check geneinfo
    if (nrow(geneinfo)!=0){
      if (unique(geneinfo$chrom) != b){
        stop("Imputed expression not by chromosome or not in correct order")
      }

      # select genes
      geneinfo$keep <- 1
      if (!is.null(selectid)){
        geneinfo[!(geneinfo$id %in% selectid), "keep"] <- 0
      }
    }

    # get regionlist and boundary_genes for regions in the chromosome, and update weight_list
    res <- assign_region_ids(regioninfo,geneinfo,snpinfo,weight_list,minvar,adjust_boundary)
    regionlist <- c(regionlist, res$regionlist)
    weight_list <- res$weight_list
    boundary_genes <- rbind(boundary_genes,res$boundary_genes)

    loginfo("No. regions with at least one SNP/gene for chr%s: %d",
            b, length(regionlist_chr))
  }

  if (maxSNP < Inf){
    loginfo("Trim regions with SNPs more than %s", maxSNP)

    if ("z" %in% colnames(select)) {
      # z score is given, trim snps with lower |z|
      for (region_tag in names(regionlist)){
        if (length(regionlist[[region_tag]][["sid"]]) > maxSNP){
          idx <- match(regionlist[[region_tag]][["sid"]], select[, "id"])
          z.abs <- abs(select[idx, "z"])
          ifkeep <- rank(-z.abs) <= maxSNP
          #regionlist[[region_tag]][["sidx"]] <- regionlist[[region_tag]][["sidx"]][ifkeep]
          regionlist[[region_tag]][["sid"]] <- regionlist[[region_tag]][["sid"]][ifkeep]
        }
      }
    } else{
      # if no z score information, randomly select snps
      for (region_tag in names(regionlist)){
        if (length(regionlist[[region_tag]][["sid"]]) > maxSNP){
          n.ori <- length(regionlist[[region_tag]][["sid"]])
          ifkeep <- rep(F, n.ori)
          set.seed <- 99
          ifkeep[sample.int(n.ori, size = maxSNP)] <- T
          #regionlist[[region_tag]][["sidx"]] <-  regionlist[[region_tag]][["sidx"]][ifkeep]
          regionlist[[region_tag]][["sid"]] <-  regionlist[[region_tag]][["sid"]][ifkeep]
        }
      }
    }
  }

  #colnames(boundary_genes) <- c("gene","chrom","region1","region2")
  #regionlsit <- compute_correlations(regionlist, ncore)
  return(list(regionlist=regionlist,
              weight_list=weight_list,
              boundary_genes=boundary_genes))
}
