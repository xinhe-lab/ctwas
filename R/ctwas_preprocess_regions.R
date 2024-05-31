
#' Preprocess region info, LD info, and SNP info from all regions and map SNPs to regions.
#'
#' @param region_info a data frame of region definition
#'
#' @param ref_snp_info a data frame of all variant info in the LD reference.
#'
#' @param use_LD TRUE/FALSE. If TRUE, use LD for finemapping. Otherwise, use "no-LD" version.
#'
#' @param chrom chromosomes to keep
#'
#' @param ncore The number of cores used to parallelize over regions
#'
#' @export
#'
preprocess_region_LD_snp_info <- function(region_info = NULL,
                                          ref_snp_info = NULL,
                                          use_LD = TRUE,
                                          chrom = 1:22,
                                          ncore = 1) {

  target_header <- c("chrom", "start", "stop")
  if (!all(target_header %in% colnames(region_info))){
    stop("region_info needs to contain the following columns: ",
         paste(target_header, collapse = " "))
  }
  region_info$chrom <- as.numeric(gsub("chr", "", region_info$chrom))
  region_info$start <- as.numeric(region_info$start)
  region_info$stop <- as.numeric(region_info$stop)

  # select example chromosome
  region_info <- region_info[region_info$chrom %in% chrom,]

  # assign region ids
  if (is.null(region_info$region_id)) {
    region_info$region_id <- paste0(region_info$chrom, ":", region_info$start, "-", region_info$stop)
  }
  rownames(region_info) <- NULL

  # LD info
  if (use_LD) {
    if (is.null(region_info$LD_matrix)){
      stop("Please provide paths to LD matrices in the 'LD_matrix' column of region_info.")
    }
    LD_matrix_files <- region_info$LD_matrix
    stopifnot(all(file.exists(LD_matrix_files)))
    LD_info <- data.frame(region_id = region_info$region_id, LD_matrix = LD_matrix_files)
  } else{
    LD_info <- NULL
  }

  # SNP info
  if (use_LD) {
    if (is.null(region_info$SNP_info)){
      stop("Please provide paths to SNP info tables in the 'SNP_info' column of region_info.")
    }
    snp_info_files <- region_info$SNP_info
    stopifnot(all(file.exists(snp_info_files)))
    snp_info <- lapply(snp_info_files, read_snp_info_file)
    names(snp_info) <- region_info$region_id
  } else {
    if (is.null(ref_snp_info)){
      stop("Please provide reference SNP info in 'ref_snp_info'.")
    }
    snp_info <- map_snp_info_regions(region_info, ref_snp_info, ncore = ncore)
  }

  region_info <- region_info[,c("chrom", "start", "stop", "region_id")]

  return(list(region_info = region_info,
              LD_info = LD_info,
              snp_info = snp_info))
}


#' Read SNP info from all regions and map SNPs to regions.
#'
#' @param region_info a data frame of region definition
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
    stop("Reference SNP info needs to contain the following columns: ",
         paste(target_header, collapse = " "))
  }

  ref_snp_info <- ref_snp_info[ref_snp_info$chrom %in% region_info$chrom, ]

  region_ids <- region_info$region_id
  snp_info <- parallel::mclapply(region_ids, function(region_id){
    regioninfo <- region_info[which(region_info$region_id == region_id), ]
    region_chrom <- regioninfo$chrom
    region_start <- regioninfo$start
    region_stop <- regioninfo$stop
    region_idx <- which(ref_snp_info$chrom == region_chrom & ref_snp_info$pos >= region_start & ref_snp_info$pos < region_stop)
    region_snp_info <- ref_snp_info[region_idx, ]
    rownames(region_snp_info) <- NULL
    region_snp_info
  }, mc.cores = ncore)
  names(snp_info) <- region_ids

  return(snp_info)
}
