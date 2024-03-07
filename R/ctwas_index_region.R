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
index_region <- function(regions,
                         geneinfo,
                         snpinfo,
                         weight_list,
                         minvar){
  regionlist <- list()
  boundary_genes <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(boundary_genes) <- c("gene","chr","region1","region2")
  for (i in 1:nrow(regions)){

    rn.start <- regions$start[i]
    rn.stop <- regions$stop[i]
    rn <- paste0(rn.start,"-",rn.stop)
    gidx <- which(geneinfo$p0 >= rn.start & geneinfo$p0 < rn.stop
                  & geneinfo$keep == 1) # temporarily assign to the first region if its QTLs are across boundary

    sidx <- which(snpinfo$pos >= rn.start & snpinfo$pos < rn.stop
                  & snpinfo$keep == 1 & snpinfo$thin_tag == 1)

    if (length(gidx) + length(sidx) < minvar) {next}

    gid <- geneinfo$id[gidx]
    sid <- snpinfo$id[sidx]

    minpos <- min(c(geneinfo$p0[gidx], snpinfo$pos[sidx]))
    maxpos <- max(c(geneinfo$p1[gidx], snpinfo$pos[sidx]))

    regionlist[[rn]] <- list("gid"  = gid,
                             "sid" = sid,
                             "start" = rn.start,
                             "stop" = rn.stop,
                             "minpos" = minpos,
                             "maxpos" = maxpos,
                             "LD_matrix" = regions$LD_matrix[i],
                             "SNP_info" = regions$SNP_info[i],
                             "region_tag" = regions$region_tag[i])
  }

  if (nrow(regions) >=2){
    res <- adjust_boundary(regions, weight_list, regionlist)
    regionlist <- res$regionlist
    weight_list <- res$weight_list
    boundary_genes <- res$boundary_genes
  }

  return(list(regionlist=regionlist,weight_list=weight_list,boundary_genes=boundary_genes))
}
