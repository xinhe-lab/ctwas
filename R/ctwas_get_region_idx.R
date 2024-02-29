#' Get gene and SNP index for each region
#' @description For each region, get the index for snp and gene (index
#' is location/column number in .pgen file or .expr file) located within
#' this region.
#'
#' @param regionfile regions file. Has three columns: chr, start, end. The regions file
#' should provide non overlaping regions defining LD blocks. currently does not support
#' chromsome X/Y etc.
#'
#' @param select Default is NULL, all variants will be selected. Or a vector of variant IDs,
#' or a data frame with columns id and z (id is for gene or SNP id, z is for z scores).
#' z will be used for remove SNPs if the total number of SNPs exceeds limit. See
#' parameter `maxSNP` for more information.
#' 
#' @param thin  A scalar in (0,1]. The proportion of SNPs
#' left after down sampling. Only applied on SNPs after selecting variants.
#'  
#' @param maxSNP Default is Inf, no limit for the maximum number of SNPs in a region. Or an
#' integer indicating the maximum number of SNPs allowed in a region. This
#' parameter is useful when a region contains many SNPs and you don't have enough memory to
#' run the program. In this case, you can put a limit on the number of SNPs in the region.
#' If z scores are given in the parameter `select`, i.e. a data frame with columns id and z is
#' provided, SNPs are ranked based on |z| from high to low and only the top `maxSNP` SNPs
#' are kept. If only variant ids are provided, then `maxSNP` number of SNPs will be chosen
#' randomly.
#'  
#' @param minvar minimum number of variatns in a region
#'
#' @param merge TRUE/FALSE. If TRUE, merge regions when a gene spans a region boundary (i.e. belongs to multiple regions.)
#' 
#' @param outputdir a string, the directory to store output
#' 
#' @param outname a string, the output name
#' 
#' @param ncore the number of cores used to parallelize region indexing
#' 
#' @param reuse_R_gene an option to reuse the R_gene matrix when indexing for the final rerun step
#'
#' @return A list. Items correspond to each pvarf/exprvarf. Each Item is
#'  also a list, the items in this list are for each region.
#'
#' @importFrom logging loginfo
#'
get_region_idx <- function(region_info,
                       gene_info,
                       weight_list,
                       select = NULL,
                       thin = 1,
                       maxSNP = Inf,
                       minvar = 1
                       ) {

  loginfo("No. LD regions: %s", nrow(region_info))

  if (thin > 1 | thin <= 0){
    stop("thin needs to be in (0,1]")
  }

  if (is.null(dim(select))) {
    selectid <- select
  } else {
    selectid <- select$id
  }
  
  region_info <- region_info[order(region_info$chr, region_info$start),]
  
  regionlist <- list()
  boundary_genes <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(boundary_genes) <- c("gene","chr","region1","region2")
  for (b in unique(region_info$chr)){
    # select regions in chromosome b
    regions <- region_info[region_info$chr == b, ]
     # select genes in chromosome b
    geneinfo <- gene_info[gene_info$chrom == b, ]
    # get snp info in LD in the chromosome
    snpinfo <- ctwas:::read_LD_SNP_files(regions$SNP_info)
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

    #check geneinfo
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
    
    res <- index_region(regions,geneinfo,snpinfo,weight_list,minvar)
    regionlist[[as.character(b)]] <- res$regionlist
    weight_list <- res$weight_list
    boundary_genes <- rbind(boundary_genes,res$boundary_genes)
    
    loginfo("No. regions with at least one SNP/gene for chr%s: %s",
          b, length(regionlist[[as.character(b)]]))
  }

  loginfo("Trim regions with SNPs more than %s", maxSNP)

  if ("z" %in% colnames(select)) {
    # z score is given, trim snps with lower |z|
    for (b in names(regionlist)){
      for (rn in names(regionlist[[b]])) {
        if (length(regionlist[[b]][[rn]][["sid"]]) > maxSNP){
          idx <- match(regionlist[[b]][[rn]][["sid"]], select[, "id"])
          z.abs <- abs(select[idx, "z"])
          ifkeep <- rank(-z.abs) <= maxSNP
          regionlist[[b]][[rn]][["sidx"]] <-  regionlist[[b]][[rn]][["sidx"]][ifkeep]
          regionlist[[b]][[rn]][["sid"]] <-  regionlist[[b]][[rn]][["sid"]][ifkeep]
        }
      }
    }
  } else{
    # if no z score information, randomly select snps
    for (b in 1: length(regionlist)){
      for (rn in names(regionlist[[b]])) {
        if (length(regionlist[[b]][[rn]][["sid"]]) > maxSNP){
          n.ori <- length(regionlist[[b]][[rn]][["sid"]])
          ifkeep <- rep(F, n.ori)
          set.seed <- 99
          ifkeep[sample.int(n.ori, size = maxSNP)] <- T
          regionlist[[b]][[rn]][["sidx"]] <-  regionlist[[b]][[rn]][["sidx"]][ifkeep]
          regionlist[[b]][[rn]][["sid"]] <-  regionlist[[b]][[rn]][["sid"]][ifkeep]
        }
      }
    }
  }
  #colnames(boundary_genes) <- c("gene","chr","region1","region2")
  #regionlsit <- compute_correlations(regionlist, ncore)
  return(list(regionlist=regionlist,weight_list=weight_list, boundary_genes=boundary_genes))
}
