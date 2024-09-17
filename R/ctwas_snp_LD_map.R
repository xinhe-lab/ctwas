
#' @title Map SNPs to regions using all the variants in the LD reference.
#'
#' @param region_info a data frame of region definitions, with columns:
#' "chrom", "start", "stop", "region_id".
#'
#' @param ref_snp_info a data frame of all variant info in the reference.
#'
#' @param ncore The number of cores used to parallelize over regions
#'
#' @importFrom readr parse_number
#'
#' @export
#'
create_snp_map <- function(region_info,
                           ref_snp_info,
                           ncore = 1) {

  # check inputs
  region_info <- as.data.frame(region_info)
  target_header <- c("chrom", "start", "stop")
  if (!all(target_header %in% colnames(region_info))){
    stop("region_info needs to contain the following columns: ",
         paste(target_header, collapse = " "))
  }

  ref_snp_info <- as.data.frame(ref_snp_info)
  target_header <- c("chrom", "id", "pos", "alt", "ref")
  if (!all(target_header %in% colnames(ref_snp_info))){
    stop("ref_snp_info needs to contain the following columns: ",
         paste(target_header, collapse = " "))
  }

  region_info$chrom <- parse_number(as.character(region_info$chrom))
  region_info$start <- as.numeric(region_info$start)
  region_info$stop <- as.numeric(region_info$stop)

  # assign region ids if not available
  if (is.null(region_info$region_id)) {
    region_info$region_id <- paste0(region_info$chrom, "_", region_info$start, "_", region_info$stop)
  }
  rownames(region_info) <- NULL

  region_info <- region_info[,c("chrom", "start", "stop", "region_id")]

  # map SNPs to regions
  ref_snp_info <- ref_snp_info[ref_snp_info$chrom %in% region_info$chrom, ]

  region_ids <- region_info$region_id
  snp_map <- mclapply_check(region_ids, function(region_id){
    regioninfo <- region_info[which(region_info$region_id == region_id), ]
    region_chrom <- regioninfo$chrom
    region_start <- regioninfo$start
    region_stop <- regioninfo$stop
    region_idx <- which(ref_snp_info$chrom == region_chrom & ref_snp_info$pos >= region_start & ref_snp_info$pos < region_stop)
    region_snp_info <- ref_snp_info[region_idx, ]
    rownames(region_snp_info) <- NULL
    region_snp_info
  }, mc.cores = ncore)
  names(snp_map) <- region_ids

  return(list("region_info" = region_info,
              "snp_map" = snp_map))
}


#' @title Map SNPs to regions using region meta table.
#'
#' @param region_metatable a data frame of region meta table, with columns:
#' "chrom", "start", "stop", "region_id", "LD_file", "SNP_file".
#'
#' @param snpinfo_loader_fun IF not NULL, use custom loader function to read SNP files.
#'
#' @importFrom readr parse_number
#'
#' @export
#'
create_snp_LD_map <- function(region_metatable,
                              snpinfo_loader_fun = NULL) {

  # check inputs
  region_metatable <- as.data.frame(region_metatable)
  target_header <- c("chrom", "start", "stop")
  if (!all(target_header %in% colnames(region_metatable))){
    stop("region_metatable needs to contain the following columns: ",
         paste(target_header, collapse = " "))
  }

  if (is.null(region_metatable$LD_file)){
    stop("Please provide filenames of LD matrices in the 'LD_file' column of region_metatable")
  }

  if (is.null(region_metatable$SNP_file)){
    stop("Please provide filenames of variant info in the 'SNP_file' column of region_metatable")
  }

  region_info <- region_metatable
  region_info$chrom <- parse_number(as.character(region_info$chrom))
  region_info$start <- as.numeric(region_info$start)
  region_info$stop <- as.numeric(region_info$stop)

  # assign region ids if not available
  if (is.null(region_info$region_id)) {
    region_info$region_id <- paste0(region_info$chrom, "_", region_info$start, "_", region_info$stop)
  }
  rownames(region_info) <- NULL
  region_info <- region_info[,c("chrom", "start", "stop", "region_id")]

  # map SNPs to regions
  stopifnot(all(file.exists(region_metatable$SNP_file)))
  if (!is.null(snpinfo_loader_fun)){
    snp_map <- lapply(region_metatable$SNP_file, snpinfo_loader_fun)
  } else {
    snp_map <- lapply(region_metatable$SNP_file, read_snp_info_file)
  }
  names(snp_map) <- region_metatable$region_id

  # map LD to regions
  stopifnot(all(file.exists(region_metatable$LD_file)))
  LD_map <- data.frame(region_id = region_metatable$region_id,
                       LD_file = region_metatable$LD_file,
                       SNP_file = region_metatable$SNP_file)

  return(list("region_info" = region_info,
              "snp_map" = snp_map,
              "LD_map" = LD_map))
}

