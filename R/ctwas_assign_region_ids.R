# assign gene and SNP IDs for regions in the regioninfo
assign_region_ids <- function(regioninfo,
                              geneinfo,
                              snpinfo,
                              minvar = 1,
                              mingene = 0) {

  regionlist <- list()

  # assign genes and SNPs to the region
  for (i in 1:nrow(regioninfo)){

    region_tag <- regioninfo$region_tag[i]
    region_chrom <- regioninfo$chrom[i]
    region_start <- regioninfo$start[i]
    region_stop <- regioninfo$stop[i]

    # assign genes to regions based gene p0 positions
    # for genes across region boundaries, assign to the first region, and adjust later
    gidx <- which(geneinfo$chrom == region_chrom & geneinfo$p0 >= region_start & geneinfo$p0 < region_stop
                  & geneinfo$keep == 1)

    sidx <- which(snpinfo$chrom == region_chrom & snpinfo$pos >= region_start & snpinfo$pos < region_stop
                  & snpinfo$keep == 1 & snpinfo$thin_tag == 1)

    if (length(gidx) + length(sidx) < minvar) {next}

    if (length(gidx) < mingene) {next}

    gid <- geneinfo$id[gidx]
    sid <- snpinfo$id[sidx]

    minpos <- min(c(geneinfo$p0[gidx], snpinfo$pos[sidx]))
    maxpos <- max(c(geneinfo$p1[gidx], snpinfo$pos[sidx]))

    regionlist[[region_tag]] <- list("region_tag" = region_tag,
                                     "gid" = gid,
                                     "sid" = sid,
                                     "chrom" = region_chrom,
                                     "start" = region_start,
                                     "stop" = region_stop,
                                     "minpos" = minpos,
                                     "maxpos" = maxpos)
  }

  return(regionlist)

}
