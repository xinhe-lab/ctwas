
#' Read SNP info from all regions as a data frame and map SNPs to regions.
#' if snp_info_file is available, read SNP info from the snp_info_file;
#' otherwise, read from the "SNP_info" column of region_info.
#'
#' @param region_info a data frame of region definition and associated file names
#'
#' @param snp_info a data frame of all variant info in the LD reference.
#'
#' @param ncore The number of cores used to parallelize over regions
#'
#' @importFrom data.table rbindlist
#' @importFrom parallel mclapply
#'
#' @export
#'
map_snp_info_regions <- function(region_info, snp_info = NULL, ncore = 1){

  # check input
  if (!is.null(snp_info)){
    snp_info <- unique(as.data.frame(snp_info))
    target_header <- c("chrom", "id", "pos", "alt", "ref")
    if (!all(target_header %in% colnames(snp_info))){
      stop("SNP info needs to contain the following columns: ",
           paste(target_header, collapse = " "))
    }
  } else {
    if (!is.character(region_info$SNP_info) || is.null(region_info$SNP_info)){
      stop("SNP_info in region_info is needed when snp_info_file is not available")
    }
  }

  # if snp_info_file is available, read SNP_info from snp_info_file,
  # otherwise, use SNP info from SNP_info column of region_info.
  region_ids <- region_info$region_id

  snp_info_list <- parallel::mclapply(region_ids, function(region_id){
    regioninfo <- region_info[which(region_info$region_id == region_id), ]
    if (!is.null(snp_info)){
      region_chrom <- regioninfo$chrom
      region_start <- regioninfo$start
      region_stop <- regioninfo$stop
      region_idx <- which(snp_info$chrom == region_chrom & snp_info$pos >= region_start & snp_info$pos < region_stop)
      region_snp_info <- snp_info[region_idx, ]
    } else {
      snp_info_files <- unlist(strsplit(regioninfo$SNP_info, split = ";"))
      stopifnot(all(file.exists(snp_info_files)))
      region_snp_info <- read_snp_info_files(snp_info_files)
    }
    region_snp_info
  }, mc.cores = ncore)
  names(snp_info_list) <- region_ids

  snp_info <- data.table::rbindlist(snp_info_list, idcol = "region_id")
  snp_info <- unique(as.data.frame(snp_info))
  rownames(snp_info) <- NULL

  return(snp_info)
}

