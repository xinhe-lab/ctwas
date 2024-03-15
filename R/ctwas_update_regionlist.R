#' Update regionlist with full SNPs
#'
#' @param regionlist a list of region gene IDs and SNP IDs and associated file names
#'
#' @param region_tags region tags to specify which regions to update. If NULL, update all regions in regionlist
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
#' @importFrom logging loginfo
#'
#' @return updated regionlist with full SNPs
#'
#' @export
#'
update_regionlist_fullSNPs <- function(regionlist,
                                       region_tags = NULL,
                                       select = NULL,
                                       maxSNP = Inf) {

  # if region_tags is null, update all regions in regionlist
  if (is.null(region_tags)){
    region_tags <- names(regionlist)
  }

  loginfo("Update regionlist for %d regions", length(region_tags))

  # update SNP IDs for each region
  for (region_tag in region_tags){

    region_chrom <- regionlist[[region_tag]][["chrom"]]
    region_start <- regionlist[[region_tag]][["start"]]
    region_stop <- regionlist[[region_tag]][["stop"]]

    # update sid in the region with all SNPs in the region
    snpinfo <- read_LD_SNP_files(regionlist[[region_tag]][["SNP_info"]]) #ctwas
    regionlist[[region_tag]][["sid"]] <- snpinfo$id

    # update minpos and maxpos in the region
    regionlist[[region_tag]][["minpos"]] <- min(snpinfo$pos)
    regionlist[[region_tag]][["maxpos"]] <- max(snpinfo$pos)

    # Trim regions with SNPs more than maxSNP
    if (maxSNP < Inf){
      if ("z" %in% colnames(select)) {
        # z score is given, trim snps with lower |z|
        if (length(regionlist[[region_tag]][["sid"]]) > maxSNP){
          idx <- match(regionlist[[region_tag]][["sid"]], select[, "id"])
          z.abs <- abs(select[idx, "z"])
          ifkeep <- rank(-z.abs) <= maxSNP
          regionlist[[region_tag]][["sid"]] <- regionlist[[region_tag]][["sid"]][ifkeep]
        }
      } else{
        # if no z score information, randomly select snps
        if (length(regionlist[[region_tag]][["sid"]]) > maxSNP){
          n.ori <- length(regionlist[[region_tag]][["sid"]])
          ifkeep <- rep(FALSE, n.ori)
          set.seed <- 99
          ifkeep[sample.int(n.ori, size = maxSNP)] <- TRUE
          regionlist[[region_tag]][["sid"]] <-  regionlist[[region_tag]][["sid"]][ifkeep]
        }
      }
    }
  }

  return(regionlist)

}
