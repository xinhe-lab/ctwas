
#' Prepare .pvar file
#'
#' @param pgenf pgen file
#' .pvar file format: https://www.cog-genomics.org/plink/2.0/formats#pvar
#'
#' @param outputdir a string, the directory to store output
#'
#' @return corresponding pvar file
#'
#' @importFrom tools file_ext file_path_sans_ext
#'
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
        pvar <- data.table::fread(pvarf, header = F)

        if (ncol(pvar) == 6) {
          colnames(pvar) <- c('#CHROM', 'ID', 'CM', 'POS', 'ALT', 'REF')
        } else if (ncol(pvar) == 5){
          colnames(pvar) <- c('#CHROM', 'ID', 'POS', 'ALT', 'REF')
        } else {
          stop(".pvar file has incorrect format")
        }

        data.table::fwrite(pvar, file = pvarf2 , sep="\t", quote = F)
      }
    }

  } else if (file_ext(pgenf) == "bed"){
    # .bim file has no header
    pvarf <- paste0(file_path_sans_ext(pgenf), ".bim")
    pvarf2 <-  file.path(outputdir, paste0(basename(file_path_sans_ext(pgenf)), ".hbim"))

    if (!file.exists(pvarf2)){
      pvar <- data.table::fread(pvarf, header = F)
      colnames(pvar) <- c('#CHROM', 'ID', 'CM', 'POS', 'ALT', 'REF')
      data.table::fwrite(pvar, file = pvarf2 , sep="\t", quote = F)
    }
    pvarfout <- pvarf2
  } else {
    stop("Unrecognized genotype input format")
  }

  return(pvarfout)
}

#' Read .pgen file into R
#'
#' @param pgenf .pgen file or .bed file
#'
#' @param pvarf .pvar file or .bim file with have proper
#'  header.  Matching `pgenf`.
#'
#' @return  A matrix of allele count for each variant (columns) in each sample
#'  (rows). ALT allele in pvar file is counted (A1 allele in .bim file is the ALT
#'   allele).
#'
#' @importFrom pgenlibr NewPvar
#' @importFrom pgenlibr NewPgen
#' @importFrom tools file_ext file_path_sans_ext
#'
prep_pgen <- function(pgenf, pvarf){

  pvar <- pgenlibr::NewPvar(pvarf)

  if (file_ext(pgenf) == "pgen"){
    pgen <- pgenlibr::NewPgen(pgenf, pvar = pvar)

  } else if (file_ext(pgenf) == "bed"){
    famf <- paste0(file_path_sans_ext(pgenf), ".fam")
    fam <- data.table::fread(famf, header = F)
    raw_s_ct <- nrow(fam)
    pgen <- pgenlibr::NewPgen(pgenf, pvar = pvar, raw_sample_ct = raw_s_ct)

  } else{
    stop("unrecognized input")
  }

  return(pgen)
}

#' Read pgen file into R
#'
#' @param pgen .pgen file or .bed file
#'
#' @param variantidx variant index. If NULL, all variants will be extracted.
#'
#' @return A matrix, columns are allele count for each SNP, rows are
#'  for each sample.
#'
#' @importFrom pgenlibr GetVariantCt
#' @importFrom pgenlibr ReadList
#'
read_pgen <- function(pgen, variantidx = NULL, meanimpute = F ){
  if (is.null(variantidx)){
    variantidx <- 1: pgenlibr::GetVariantCt(pgen)}

  pgenlibr::ReadList(pgen,
                     variant_subset = variantidx,
                     meanimpute = meanimpute)
}


#' Read .pvar file into R
#' @param pvarf .pvar file with proper format: https://www.cog-genomics.org/plink/2.0/formats#pvar
#'
#' @return A data.table. variant info
#'
read_pvar <- function(pvarf){
  pvar <- data.table::fread(pvarf, skip = "#CHROM")
  pvar <- dplyr::rename(pvar, "chrom" = "#CHROM", "pos" = "POS",
                        "alt" = "ALT", "ref" = "REF", "id" = "ID")
  pvar <- pvar[, c("chrom", "id", "pos", "alt", "ref")]
  return(pvar)
}

#' Read .bim file into R
#' @param bimf .bim file with proper format: https://www.cog-genomics.org/plink/2.0/formats#bim
#'
#' @return A data.table. variant info
#'
read_bim <- function(bimf) {
  bim <- data.table::fread(bimf)
  colnames(bim) <- c("chr", "id", "cm", "pos", "alt", "ref")
  return(bim)
}

#' Read variant information from .pvar or .bim file into R
#' @param var_info_file .pvar or .bim file with proper format:
#' .pvar: https://www.cog-genomics.org/plink/2.0/formats#pvar
#' .bim: https://www.cog-genomics.org/plink/2.0/formats#bim
#'
#' @return A data.table. variant info
#' @importFrom tools file_ext
#'
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
