#' Expand regionlist with full SNPs
#'
#' @param regionlist a list of region gene IDs and SNP IDs and associated file names
#'
#' @param filter_z_ids If TRUE, only keep SNPs and genes in z_snp and z_gene.
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
                              filter_z_ids = TRUE,
                              trim_by = c("random", "z"),
                              maxSNP = Inf,
                              seed = 99) {

  trim_by <- match.arg(trim_by)

  region_tags <- names(regionlist)

  # combine z-scores of SNPs and genes (for selecting SNPs and genes in regionlist)
  zdf <- combine_z(z_snp, z_gene)

  if (is.null(gene_info)){
    gene_info <- get_gene_info(z_gene, weights, region_info)
  }

  loginfo("Update regionlist for %d regions with full SNPs", length(region_tags))

  # update SNP IDs for each region
  for (region_tag in region_tags){

    # update sid in the region with all SNPs in the region
    snpinfo <- read_LD_SNP_files(regionlist[[region_tag]][["SNP_info"]]) #ctwas

    # select SNPs
    snpinfo$keep <- rep(1, nrow(snpinfo))
    if (isTRUE(filter_z_ids)){
      # remove SNPs not in z_snp
      snpinfo$keep[!(snpinfo$id %in% z_snp$id)] <- 0
    }
    sid <- snpinfo$id[snpinfo$keep == 1]
    sidx <- match(sid, snpinfo$id)

    # update sid in the region
    regionlist[[region_tag]][["sid"]] <- sid

    # update minpos and maxpos in the region
    regionlist[[region_tag]][["minpos"]] <- min(c(regionlist[[region_tag]][["minpos"]], snpinfo$pos[sidx]))
    regionlist[[region_tag]][["maxpos"]] <- max(c(regionlist[[region_tag]][["maxpos"]], snpinfo$pos[sidx]))

    # Trim regions with SNPs more than maxSNP
    if (maxSNP < Inf){
      if (trim_by == "z") {
        # trim SNPs with lower |z|
        for (region_tag in names(regionlist)){
          if (length(regionlist[[region_tag]][["sid"]]) > maxSNP){
            loginfo("Trim region %s with SNPs more than %s", region_tag, maxSNP)
            idx <- match(regionlist[[region_tag]][["sid"]], z_snp$id)
            z.abs <- abs(zdf[idx, "z"])
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
  }

  return(regionlist)

}
