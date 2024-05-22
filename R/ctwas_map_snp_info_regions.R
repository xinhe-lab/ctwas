
#' Read SNP info from all regions and map SNPs to regions.
#'
#' @param region_info a data frame of region definition and associated file names
#'
#' @param ref_snp_info a data frame of all variant info in the LD reference.
#'
#' @param ncore The number of cores used to parallelize over regions
#'
#' @importFrom parallel mclapply
#'
#' @export
#'
map_snp_info_regions <- function(region_info,
                                 ref_snp_info,
                                 ncore = 1){
  # check input
  ref_snp_info <- unique(as.data.frame(ref_snp_info))
  target_header <- c("chrom", "id", "pos", "alt", "ref")
  if (!all(target_header %in% colnames(ref_snp_info))){
    stop("SNP info needs to contain the following columns: ",
         paste(target_header, collapse = " "))
  }

  region_ids <- region_info$region_id
  snp_info <- parallel::mclapply(region_ids, function(region_id){
    regioninfo <- region_info[which(region_info$region_id == region_id), ]
    region_chrom <- regioninfo$chrom
    region_start <- regioninfo$start
    region_stop <- regioninfo$stop
    region_idx <- which(ref_snp_info$chrom == region_chrom & ref_snp_info$pos >= region_start & ref_snp_info$pos < region_stop)
    region_snp_info <- ref_snp_info[region_idx, ]
    region_snp_info
  }, mc.cores = ncore)
  names(snp_info) <- region_ids

  return(snp_info)
}

