
# Read a single SNP info file as a data frame
#' @importFrom data.table fread
read_snp_info_file <- function (file){
  snp_info <- as.data.frame(fread(file, header = TRUE))
  target_header <- c("chrom", "id", "pos", "alt", "ref")
  if (!all(target_header %in% colnames(snp_info))){
    stop("The SNP info file needs to contain the following columns: ",
         paste(target_header, collapse = " "))
  }
  return(snp_info)
}

# Read all SNP info files as a data frame
read_snp_info_files <- function (files){
  snp_info <- do.call(rbind, lapply(files, read_snp_info_file))
  snp_info <- unique(as.data.frame(snp_info))
  return(snp_info)
}

#' Load LD matrix
#' @param file path to LD matrix
#' @param format file format for LD matrix. If "custom", use a user defined
#' \code{LD_loader_fun()} function to load LD matrix.
#' @param LD_loader_fun a user defined function to load LD matrix
#'
#' @importFrom utils read.csv
#' @importFrom Matrix readMM
#' @importFrom data.table fread
#' @importFrom tools file_ext
#' @export
load_LD <- function (file,
                     format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
                     LD_loader_fun) {

  format <- match.arg(format)

  if (format == "rds"){
    R <- readRDS(file)
  } else if (format == "rdata"){
    R <- get(load(file))
  } else if (format == "mtx"){
    R <- readMM(file)
  } else if (format %in% c("txt", "csv")){
    R <- fread(file)
  } else if (format == "custom") {
    # use LD_loader_fun() function to load LD matrix
    R <- LD_loader_fun(file)
  }

  return(as.matrix(R))
}

# Prepare .pvar file
#
# @param pgenf pgen file
# .pvar file format: https://www.cog-genomics.org/plink/2.0/formats#pvar
#
# @param outputdir a string, the directory to store output
#
# @return corresponding pvar file
#
#' @importFrom utils read.table
#' @importFrom data.table fread fwrite
#' @importFrom tools file_ext file_path_sans_ext
prep_pvar <- function(pgenf, outputdir = getwd()){

  if (file_ext(pgenf) == "pgen"){
    pvarf <- paste0(file_path_sans_ext(pgenf), ".pvar")
    pvarf2 <-  paste0(outputdir, basename(file_path_sans_ext(pgenf)), ".hpvar")

    # pgenlib can't read pvar without header, check if header present
    firstl <- read.table(file = pvarf, header = F, comment.char = '',
                         nrows = 1, stringsAsFactors = F)

    if (substr(firstl[1,1],1,1) == '#') {
      pvarfout <- pvarf
    } else {
      pvarfout <- pvarf2

      if (!file.exists(pvarf2)) {
        pvar <- fread(pvarf, header = F)

        if (ncol(pvar) == 6) {
          colnames(pvar) <- c('#CHROM', 'ID', 'CM', 'POS', 'ALT', 'REF')
        } else if (ncol(pvar) == 5){
          colnames(pvar) <- c('#CHROM', 'ID', 'POS', 'ALT', 'REF')
        } else {
          stop(".pvar file has incorrect format")
        }

        fwrite(pvar, file = pvarf2 , sep="\t", quote = F)
      }
    }

  } else if (file_ext(pgenf) == "bed"){
    # .bim file has no header
    pvarf <- paste0(file_path_sans_ext(pgenf), ".bim")
    pvarf2 <-  file.path(outputdir, paste0(basename(file_path_sans_ext(pgenf)), ".hbim"))

    if (!file.exists(pvarf2)){
      pvar <- fread(pvarf, header = F)
      colnames(pvar) <- c('#CHROM', 'ID', 'CM', 'POS', 'ALT', 'REF')
      fwrite(pvar, file = pvarf2 , sep="\t", quote = F)
    }
    pvarfout <- pvarf2
  } else {
    stop("Unrecognized genotype input format")
  }

  return(pvarfout)
}

# Read .pgen file into R
#
# @param pgenf .pgen file or .bed file
#
# @param pvarf .pvar file or .bim file with have proper
#  header.  Matching `pgenf`.
#
# @return  A matrix of allele count for each variant (columns) in each sample
#  (rows). ALT allele in pvar file is counted (A1 allele in .bim file is the ALT
#   allele).
#
#' @importFrom data.table fread
#' @importFrom pgenlibr NewPgen NewPvar
#' @importFrom tools file_ext file_path_sans_ext
prep_pgen <- function(pgenf, pvarf){

  pvar <- NewPvar(pvarf)

  if (file_ext(pgenf) == "pgen"){
    pgen <- NewPgen(pgenf, pvar = pvar)

  } else if (file_ext(pgenf) == "bed"){
    famf <- paste0(file_path_sans_ext(pgenf), ".fam")
    fam <- fread(famf, header = F)
    raw_s_ct <- nrow(fam)
    pgen <- NewPgen(pgenf, pvar = pvar, raw_sample_ct = raw_s_ct)

  } else{
    stop("unrecognized input")
  }

  return(pgen)
}

# Read pgen file into R
#
# @param pgen .pgen file or .bed file
#
# @param variantidx variant index. If NULL, all variants will be extracted.
#
# @return A matrix, columns are allele count for each SNP, rows are
#  for each sample.
#
#' @importFrom pgenlibr GetVariantCt
#' @importFrom pgenlibr ReadList
read_pgen <- function(pgen, variantidx = NULL, meanimpute = F ){
  if (is.null(variantidx)){
    variantidx <- 1: GetVariantCt(pgen)}

  ReadList(pgen,
           variant_subset = variantidx,
           meanimpute = meanimpute)
}


# Read .pvar file into R
# @param pvarf .pvar file with proper format: https://www.cog-genomics.org/plink/2.0/formats#pvar
#
# @return A data.table. variant info
#' @importFrom data.table fread
#' @importFrom dplyr rename
read_pvar <- function(pvarf){
  pvar <- fread(pvarf, skip = "#CHROM")
  pvar <- rename(pvar, "chrom" = "#CHROM", "pos" = "POS",
                 "alt" = "ALT", "ref" = "REF", "id" = "ID")
  pvar <- pvar[, c("chrom", "id", "pos", "alt", "ref")]
  return(pvar)
}

# Read .bim file into R
# @param bimf .bim file with proper format: https://www.cog-genomics.org/plink/2.0/formats#bim
#
# @return A data.table. variant info
#' @importFrom data.table fread
read_bim <- function(bimf) {
  bim <- fread(bimf)
  colnames(bim) <- c("chrom", "id", "cm", "pos", "alt", "ref")
  return(bim)
}

# Read variant information from .pvar or .bim file into R
# @param var_info_file .pvar or .bim file with proper format:
# .pvar: https://www.cog-genomics.org/plink/2.0/formats#pvar
# .bim: https://www.cog-genomics.org/plink/2.0/formats#bim
#
# @return A data.table. variant info
#' @importFrom tools file_ext
read_var_info <- function(var_info_file){
  if (file_ext(var_info_file) == "pvar"){
    var_info <- read_pvar(var_info_file)
  } else if (file_ext(var_info_file) == "bim"){
    var_info <- read_bim(var_info_file)
  } else{
    stop("unrecognized input")
  }
  return(var_info)
}


# check mclapply result
check_mc_res <- function(x){
  if (any(sapply(x, is.null))) {
    stop("Not all cores returned results. Try rerun with bigger memory or fewer cores")
  }
}

