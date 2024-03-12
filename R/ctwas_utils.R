#' Function to identify consecutive regions and label them
label_merged_regions <- function(df) {
  # Identify consecutive regions
  df$merged <- c(FALSE, df$start[-1] == df$stop[-length(df$stop)])
  df$label <- cumsum(!df$merged)

  # Remove the consecutive column as it's no longer needed
  df$merged <- NULL

  return(df)
}

#' read all LD SNP info files as a data frame
read_LD_SNP_files <- function(files){
  if (length(files)>0){
    LD_SNP <- do.call(rbind, lapply(files, read_LD_SNP_file))
  } else {
    LD_SNP <- data.table::data.table(chrom=as.integer(), id=as.character(),
                                     pos=as.integer(), alt=as.character(),
                                     ref=as.character(), variance=as.numeric())
  }
  LD_SNP
}

#' read a single LD SNP info file as a data frame
read_LD_SNP_file <- function(file){
  LD_SNP <- data.table::fread(file, header = T)
  target_header <- c("chrom", "id", "pos", "alt", "ref")

  if (!all(target_header %in% colnames(LD_SNP))){
    stop("The .Rvar file needs to contain the following columns: ",
         paste(target_header, collapse = " "))
  }
  if (length(unique(LD_SNP$chrom)) != 1){
    stop("LD region needs to be on only one chromosome.")
  }

  return(LD_SNP)
}

# # check SNP_info files from the region_info table, and return a list of region_info tables
# # adapted from old write_ld_Rf() function
# get_region_info_list <- function(region_info){
#
#   # read and check SNP info (Rvar) data and extract min and max positions in each region
#   for (i in 1:nrow(region_info)){
#     region_snp_info <- read_LD_SNP_file(region_info$SNP_info[i])
#     region_info$SNP_start[i] <- min(region_snp_info$pos)
#     region_info$SNP_stop[i] <- max(region_snp_info$pos) + 1
#   }
#
#   region_info <- data.frame(region_info, stringsAsFactors = F)
#   region_info <- transform(region_info,
#                            chr = as.integer(chr),
#                            SNP_start = as.integer(SNP_start),
#                            SNP_stop = as.integer(SNP_stop))
#
#   # save a data frame for each chromosome with region_name
#   region_info_list <- vector("list", length = 22)
#
#   for (b in sort(unique(region_info$chr))) {
#     region_info_chr <- region_info[region_info$chr == b, , drop = F]
#     region_info_chr <- region_info_chr[order(region_info_chr$start), ]
#     region_info_chr$region_name <- 1:nrow(region_info_chr)
#     region_info_chr$region_tag <- paste0(region_info_chr$start, "-", region_info_chr$stop)
#     region_info_list[[b]] <- region_info_chr
#   }
#
#   return(region_info_list)
# }


#' read LD by file format
read_LD <- function(file, format = c("rds", "rdata", "csv", "txt", "tsv")) {
  format <- match.arg(format)

  # if format is missing, try to guess format by file extension
  if (missing(format)) {
    file_ext_lower <- tolower(tools::file_ext(file))

    if (file_ext_lower == "rds"){
      format <- "rds"
    } else if (file_ext_lower %in% c("rdata", "rd", "rda", "rdat")){
      format <- "rdata"
    } else if (file_ext_upper %in% c("csv", "csv.gz")) {
      format <- "csv"
    } else if (file_ext_lower %in% c("txt", "txt.gz")){
      format <- "txt"
    } else if (file_ext_lower %in% c("tsv", "tsv.gz")){
      format <- "tsv"
    } else {
      stop("Unknown LD file format!")
    }
  }

  if (format == "rds"){
    res <- readRDS(file)
  } else if (format == "rdata"){
    res <- get(load(file))
  } else if (format == "csv"){
    res <- as.matrix(read.csv(file, sep=",", row.names=1))
  } else if (format %in% c("txt", "tsv")){
    res <- as.matrix(data.table::fread(file))
  } else {
    stop("Unknown file format!")
  }

  return(res)
}

