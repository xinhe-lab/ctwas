#' Prepare .pvar file
#' @param pgenf pgen file
#'  .pvar file format: https://www.cog-genomics.org/plink/2.0/formats#pvar
#' @return corresponding pvar file
#'
#' @importFrom tools file_ext file_path_sans_ext
#' @export
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

  pvarfout
}

#' Read .pvar file into R
#' @param pvarf .pvar file or .bim file with have proper
#'  .pvar file format: https://www.cog-genomics.org/plink/2.0/formats#pvar
#'
#' @return A data.table. variant info
#' @export
read_pvar <- function(pvarf){

  pvardt <- data.table::fread(pvarf, skip = "#CHROM")
  pvardt <- dplyr::rename(pvardt, "chrom" = "#CHROM", "pos" = "POS",
                "alt" = "ALT", "ref" = "REF", "id" = "ID")
  pvardt <- pvardt[, c("chrom", "id", "pos", "alt", "ref")]
  pvardt
}

#' Read .pgen file into R
#' @param pgenf .pgen file or .bed file
#' @param pvarf .pvar file or .bim file with have proper
#'  header.  Matching `pgenf`.
#' @return  A matrix of allele count for each variant (columns) in each sample
#'  (rows). ALT allele in pvar file is counted (A1 allele in .bim file is the ALT
#'   allele).
#'
#' @importFrom pgenlibr NewPvar
#' @importFrom pgenlibr NewPgen
#' @importFrom tools file_ext file_path_sans_ext
#'
#' @export
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

  pgen
}

#' Read pgen file into R
#' @param variantidx variant index. If NULL, all variants will be extracted.
#' @return A matrix, columns are allele count for each SNP, rows are
#'  for each sample.
#' @importFrom pgenlibr GetVariantCt
#' @importFrom pgenlibr ReadList
#' @export
read_pgen <- function(pgen, variantidx = NULL, meanimpute = F ){
  if (is.null(variantidx)){
    variantidx <- 1: pgenlibr::GetVariantCt(pgen)}

  pgenlibr::ReadList(pgen,
                     variant_subset = variantidx,
                     meanimpute = meanimpute)
}


#' Prepare .exprvar file
#' @return corresponding exprvar file
#'
#' @importFrom tools file_ext file_path_sans_ext
#' @export
prep_exprvar <- function(exprf){
  if (file_ext(exprf) == "gz"){
    exprf <- file_path_sans_ext(exprf)
  }
  exprvarf <- paste0(exprf, "var")
  exprvarf
}

#' Read .exprvar file into R
#' @return A data.table. variant info
#' @export
read_exprvar <- function(exprvarf){

  exprvar <- try(data.table::fread(exprvarf, header = T))

  if (inherits(exprvar, "try-error")){
    exprvar <-  setNames(data.table(matrix(nrow = 0, ncol = 4)),
                         c("chrom", "id", "p0", "p1"))
  }
  exprvar
}

#' Read .expr file into R
#' @param variantidx variant index. If NULL, all variants will be extracted.
#' @return A matrix, columns are imputed expression for each gene, rows are
#'  for each sample.
#' @export
read_expr <- function(exprf, variantidx = NULL){
  if (!is.null(variantidx) & length(variantidx)==0){
    return(NULL)
  } else{
    return(as.matrix(data.table::fread(exprf, header = F,
                                       select = variantidx)))
  }
}

#' @return a data frame with columns: chr region_name start stop RDS_file
read_ld_Rinfo <- function(ld_Rf){
   Rinfo <- data.table::fread(ld_Rf, header = T)
   target_header <- c("chrom", "region_name", "start", "stop", "RDS_file")
   if (all(target_header %in% colnames(Rinfo))){
     if (length(unique(Rinfo$chrom)) == 1){
       return(Rinfo)
     } else{
       stop("Each ld_Rf needs to contain regions from only one chromosome")
     }
   } else {
     stop("The ld_Rf file needs to contain the following columns: ",
          paste(target_header, collapse = " "))
   }
}

#' read variant information associated with a LD R matrix .RDS file.
#'
#' @return a data frame with columns: "chrom", "id", "pos", "alt", "ref". "alt" is
#' the coded allele
#'
#' @importFrom tools file_ext file_path_sans_ext
read_ld_Rvar_RDS <- function(ld_RDSf){
  ld_Rvarf <- paste0(file_path_sans_ext(ld_RDSf), ".Rvar")
  ld_Rvar <- data.table::fread(ld_Rvarf, header = T)
  target_header <- c("chrom", "id", "pos", "alt", "ref")
  if (all(target_header %in% colnames(ld_Rvar))){
      return(ld_Rvar)
  } else {
    stop("The .Rvar file needs to contain the following columns: ",
         paste(target_header, collapse = " "))
  }
}

#' read variant information for all regions in ld_Rf.
#' @return a data frame with columns: "chrom", "id", "pos", "alt", "ref"
read_ld_Rvar <- function(ld_Rf){
  Rinfo <- read_ld_Rinfo(ld_Rf)
  ld_Rvar <- do.call(rbind, lapply(Rinfo$RDS_file, read_ld_Rvar_RDS))
  ld_Rvar <- unique(ld_Rvar)
  ld_Rvar
}



