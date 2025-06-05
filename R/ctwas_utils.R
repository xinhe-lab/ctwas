
#' @title Read a single SNP info file as a data frame
#'
#' @param file SNP info file name.
#'
#' @param snpinfo_loader_fun a user defined function to load SNP information file,
#' if SNP information files are not in standard cTWAS reference format.
#'
#' @importFrom data.table fread
#'
#' @return a data frame of SNP info
#'
#' @export
read_snp_info_file <- function (file,
                                snpinfo_loader_fun = NULL){

  if (!is.null(snpinfo_loader_fun)){
    snp_info <- snpinfo_loader_fun(file)
  } else {
    snp_info <- fread(file, header = TRUE)
  }
  snp_info <- unique(as.data.frame(snp_info))

  target_header <- c("chrom", "id", "pos", "alt", "ref")
  if (!all(target_header %in% colnames(snp_info))){
    stop("The SNP file needs to contain the following columns: ",
         paste(target_header, collapse = " "))
  }
  return(snp_info)
}

#' @title Read multiple single SNP info files as a data frame
#'
#' @param files a vector of file names for SNP info in all regions
#'
#' @param snpinfo_loader_fun a user defined function to load SNP information file,
#' if SNP information files are not in standard cTWAS reference format.
#'
#' @importFrom data.table fread
#'
#' @return a data frame of SNP info from all files
#'
#' @export
read_snp_info_files <- function (files,
                                 snpinfo_loader_fun = NULL){
  if (!is.null(snpinfo_loader_fun)){
    snp_info <- do.call(rbind, lapply(files, snpinfo_loader_fun))
  } else {
    snp_info <- do.call(rbind, lapply(files, read_snp_info_file))
  }
  snp_info <- unique(as.data.frame(snp_info))

  target_header <- c("chrom", "id", "pos", "alt", "ref")
  if (!all(target_header %in% colnames(snp_info))){
    stop("The SNP file needs to contain the following columns: ",
         paste(target_header, collapse = " "))
  }
  return(snp_info)
}

#' @title Load LD matrix
#'
#' @param file LD matrix file name.
#'
#' @param format file format for LD matrix. If "custom", use a user defined
#' \code{LD_loader_fun()} function to load LD matrix.
#'
#' @param LD_loader_fun a user defined function to load LD matrix
#'
#' @importFrom utils read.csv
#' @importFrom Matrix readMM
#' @importFrom data.table fread
#' @importFrom tools file_ext
#'
#' @export
#'
load_LD <- function (file,
                     format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
                     LD_loader_fun = NULL) {

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


#' @title Convert variant IDs from Open GWAS format or PredictDB weight format
# to UKB reference variant format.
#'
#' @param varIDs a vector of variant IDs from Open GWAS format
#' ("chr_pos_ref_alt") or PredictDB weight format
#' (“chr_pos_ref_alt_build” or “chr:pos_ref_alt_build”)
#'
#' @return a vector of variant IDs in the format of "chr:pos_ref_alt"
#' @export
#'
convert_to_ukb_varIDs <- function(varIDs){

  varID_list <- strsplit(varIDs, split = "_|:")
  new_varIDs <- unlist(lapply(varID_list, function(x){
    if (length(x) >= 4) {
      sprintf("%s:%s_%s_%s", x[1], x[2], x[3], x[4])
    } else {
      x
    }
  }))

  return(new_varIDs)
}


#' @title Check the numbers of SNPs in snp_map, z_snp and weights
#'
#' @param snp_map a list of data frames with SNP-to-region map for the reference.
#'
#' @param z_snp A data frame with columns: "id", "z", giving the z-scores for SNPs.
#'
#' @param weights a list of pre-processed prediction weights.
#'
#' @export
check_n_snps <- function(snp_map, z_snp, weights){

  if (!inherits(snp_map,"list"))
    stop("'snp_map' should be a list.")

  if (anyNA(z_snp))
    stop("z_snp contains missing values!")

  if (!inherits(weights,"list"))
    stop("'weights' should be a list object.")

  if (any(sapply(weights, is.null)))
    stop("weights contain NULL, remove empty weights!")

  # count the number of SNPs in snp_map
  # combine snp_map into a data frame of snp_info
  snp_info <- as.data.frame(rbindlist(snp_map, idcol = "region_id"))
  n_snps_snp_map <- length(unique(snp_info$id))
  loginfo("Number of SNPs in snp_map = %d.", n_snps_snp_map)

  # count the number of SNPs in z_snp
  n_snps_z_snp <- length(unique(z_snp$id))
  loginfo("Number of SNPs in z_snp = %d.", n_snps_z_snp)

  # count the number of SNPs in weights

  # get gene info from weights
  gene_info <- get_gene_info(weights)

  n_wgt <- sapply(weights, "[[", "n_wgt")
  loginfo("Average number of SNPs in weights per gene = %.1f", mean(n_wgt))

  # count by chromosome, avoid memory issues when there are too many SNPs in weights
  n_snps_weights <- 0
  for(chrom in unique(gene_info$chrom)){
    tmp_gids <- gene_info[gene_info$chrom == chrom, "id"]
    tmp_weights <- weights[tmp_gids]
    n_snps_tmp_weights <- length(unique(unlist(lapply(tmp_weights, function(x){rownames(x[["wgt"]])}))))
    n_snps_weights <- n_snps_weights + n_snps_tmp_weights
  }
  loginfo("Total number of SNPs in weights = %d.", n_snps_weights)

  frac_snps_in_weights <- n_snps_weights/n_snps_z_snp

  if (frac_snps_in_weights == 1)
    stop("Error: all SNPs (from z_snp) are in weights! ")

  if (frac_snps_in_weights > 0.5)
    logwarn("%.2f%% SNPs (from z_snp) are in weights!", frac_snps_in_weights*100)

}


# run mclapply and check to see if NULL was found in mclapply output.
#' @importFrom parallel mclapply
#' @importFrom logging logwarn
mclapply_check <- function(X, FUN, mc.cores = 1, stop_if_missing = FALSE){
  if (length(X) <= 1 || mc.cores == 1) {
    res <- lapply(X, FUN)
  } else {
    res <- mclapply(X, FUN, mc.cores = mc.cores)
    if (any(sapply(res, is.null))) {
      msg <- paste("'NULL' found in mclapply output. Results may be incomplete!",
                   "Try more memory or ncore = 1.")
      if (stop_if_missing) {
        stop(msg)
      } else {
        logwarn(msg)
      }
    }
  }
  return(res)
}
