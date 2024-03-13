#' Update regionlist with full SNPs
#'
#' @param regionlist a list of region gene IDs and SNP IDs and associated file names
#'
#' @param select Default is NULL, all variants will be selected.
#' Or a data frame with columns id and z (id is for gene or SNP id, z is for z scores).
#' z will be used for remove SNPs if the total number of SNPs exceeds limit. See
#' parameter `maxSNP` for more information.
#'
#' @param maxSNP Inf or integer. Maximum number of SNPs in a region. Default is
#' Inf, no limit. This can be useful if there are many SNPs in a region and you don't
#' have enough memory to run the program.
#'
#' @param minvar minimum number of SNPs (and genes) in a region
#'
#' @importFrom logging loginfo
#'
#' @return updated regionlist
#'
#' @export
#'
update_regionlist <- function(regionlist,
                              region_tags,
                              select = NULL,
                              maxSNP = Inf,
                              minvar = 1) {

  # update SNP IDs for each region
  for (region_tag in region_tags){

    region_chrom <- regionlist[[region_tag]][["chrom"]]
    region_start <- regionlist[[region_tag]][["start"]]
    region_stop <- regionlist[[region_tag]][["stop"]]
    minpos <- regionlist[[region_tag]][["minpos"]]
    maxpos <- regionlist[[region_tag]][["maxpos"]]
    gid <- regionlist[[region_tag]][["gid"]]

    # get snp info in LD in the chromosome
    snpinfo <- read_LD_SNP_files(regionlist[[region_tag]][["SNP_info"]]) #ctwas

    # find SNPs in the region
    sidx <- which(snpinfo$chrom == region_chrom & snpinfo$pos >= region_start & snpinfo$pos < region_stop)
    sid <- snpinfo$id[sidx]

    if (length(gidx) + length(sidx) < minvar) {next}

    # update sid in the region
    regionlist[[region_tag]][["sid"]] <- sid

    # update minpos and maxpos in the region
    regionlist[[region_tag]][["minpos"]] <- min(c(minpos, snpinfo$pos[sidx]))
    regionlist[[region_tag]][["maxpos"]] <- max(c(maxpos, snpinfo$pos[sidx]))

  }

  if (maxSNP < Inf){
    loginfo("Trim regions with SNPs more than %s", maxSNP)

    if ("z" %in% colnames(select)) {
      # z score is given, trim snps with lower |z|
      for (region_tag in region_tags){
        if (length(regionlist[[region_tag]][["sid"]]) > maxSNP){
          idx <- match(regionlist[[region_tag]][["sid"]], select[, "id"])
          z.abs <- abs(select[idx, "z"])
          ifkeep <- rank(-z.abs) <= maxSNP
          regionlist[[region_tag]][["sid"]] <- regionlist[[region_tag]][["sid"]][ifkeep]
        }
      }
    } else{
      # if no z score information, randomly select snps
      for (region_tag in region_tags){
        if (length(regionlist[[region_tag]][["sid"]]) > maxSNP){
          n.ori <- length(regionlist[[region_tag]][["sid"]])
          ifkeep <- rep(F, n.ori)
          set.seed <- 99
          ifkeep[sample.int(n.ori, size = maxSNP)] <- T
          regionlist[[region_tag]][["sid"]] <-  regionlist[[region_tag]][["sid"]][ifkeep]
        }
      }
    }
  }

  return(regionlist)

}
