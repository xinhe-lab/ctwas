# assign gene and SNP IDs for regions in the regioninfo
assign_region_ids <- function(regioninfo,
                              geneinfo,
                              snpinfo,
                              weight_list,
                              minvar = 1,
                              adjust_boundary = TRUE) {

  regionlist <- list()
  # boundary_genes <- data.frame(matrix(nrow = 0, ncol = 4))
  # colnames(boundary_genes) <- c("gene","chrom","region1","region2")

  # get gene IDs and SNP IDs for each region in regioninfo
  for (i in 1:nrow(regioninfo)){

    region_tag <- regioninfo$region_tag[i]
    region_chrom <- regioninfo$chrom[i]
    region_start <- regioninfo$start[i]
    region_stop <- regioninfo$stop[i]

    # find genes and SNPs in the region
    gidx <- which(geneinfo$chrom == region_chrom & geneinfo$p0 >= region_start & geneinfo$p0 < region_stop
                  & geneinfo$keep == 1) # temporarily assign to the first region if its QTLs are across boundary

    sidx <- which(snpinfo$chrom == region_chrom & snpinfo$pos >= region_start & snpinfo$pos < region_stop
                  & snpinfo$keep == 1 & snpinfo$thin_tag == 1)

    if (length(gidx) + length(sidx) < minvar) {next}

    gid <- geneinfo$id[gidx]
    sid <- snpinfo$id[sidx]

    minpos <- min(c(geneinfo$p0[gidx], snpinfo$pos[sidx]))
    maxpos <- max(c(geneinfo$p1[gidx], snpinfo$pos[sidx]))

    regionlist[[region_tag]] <- list("gid" = gid,
                                     "sid" = sid,
                                     "chrom" = region_chrom,
                                     "start" = region_start,
                                     "stop" = region_stop,
                                     "minpos" = minpos,
                                     "maxpos" = maxpos,
                                     "region_tag" = region_tag,
                                     "LD_matrix" = regioninfo$LD_matrix[i],
                                     "SNP_info" = regioninfo$SNP_info[i])
  }

  # identify genes across boundaries, update regionlist and weight_list according to new gene and snp assignment
  if (isTRUE(adjust_boundary) && nrow(regioninfo) >=2){
    res <- adjust_boundary(regioninfo, weight_list, regionlist)
    regionlist <- res$regionlist
    weight_list <- res$weight_list
    boundary_genes <- res$boundary_genes
  }

  return(list(regionlist=regionlist,weight_list=weight_list,boundary_genes=boundary_genes))

}
