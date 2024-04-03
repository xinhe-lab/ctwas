#' Expand regionlist with full SNPs
#'
#' @param regionlist a list of region gene IDs and SNP IDs and associated file names
#'
#' @param region_info a data frame of region definition and associated LD file names
#'
#' @param z_snp A data frame with columns: "id", "z", giving the z-scores for SNPs.
#'
#' @param trim_by remove SNPs if the total number of SNPs exceeds limit, options: "random",
#' or "z" (trim SNPs with lower |z|) See parameter `maxSNP` for more information.
#'
#' @param maxSNP Inf or integer. Maximum number of SNPs in a region. Default is
#' Inf, no limit. This can be useful if there are many SNPs in a region and you don't
#' have enough memory to run the program.
#'
#' @param seed seed for.random sampling
#'
#' @importFrom logging loginfo
#'
#' @return updated regionlist with full SNPs
#'
#' @export
#'
expand_regionlist <- function(regionlist,
                              region_info,
                              z_snp,
                              trim_by = c("z", "random"),
                              maxSNP = Inf,
                              seed = 99) {

  trim_by <- match.arg(trim_by)

  region_tags <- names(regionlist)
  loginfo("Update regionlist for %d regions with full SNPs", length(region_tags))
  pb <- txtProgressBar(min = 0, max = length(region_tags), initial = 0, style = 3)

  # update SNP IDs for each region
  for (i in 1:length(region_tags)){
    region_tag <- region_tags[i]

    # load all SNPs in the region
    snpinfo <- read_LD_SNP_files(region_info[region_info$region_tag == region_tag, "SNP_info"])

    # update sid in the region
    snpinfo$keep <- rep(1, nrow(snpinfo))
    # remove SNPs not in z_snp
    snpinfo$keep[!(snpinfo$id %in% z_snp$id)] <- 0

    sid <- snpinfo$id[snpinfo$keep == 1]
    sidx <- match(sid, snpinfo$id)
    regionlist[[region_tag]][["sid"]] <- sid

    # update minpos and maxpos in the region
    regionlist[[region_tag]][["minpos"]] <- min(c(regionlist[[region_tag]][["minpos"]], snpinfo$pos[sidx]))
    regionlist[[region_tag]][["maxpos"]] <- max(c(regionlist[[region_tag]][["maxpos"]], snpinfo$pos[sidx]))

    setTxtProgressBar(pb, i)
  }

  close(pb)

  # Trim regions with SNPs more than maxSNP
  regionlist <- trim_regionlist(regionlist, z_snp, trim_by = trim_by, maxSNP = maxSNP, seed = seed)

  return(regionlist)

}
