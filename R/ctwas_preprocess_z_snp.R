
#' @title Read GWAS summary statistics
#'
#' @param gwas A data frame of GWAS summary statistics.
#'
#' @param id Column name of the variant IDs.
#'
#' @param A1 Column name of the alternative alleles.
#'
#' @param A2 Column name of the reference alleles.
#'
#' @param z Column name of z-scores.
#'
#' @param beta Column name of effect sizes.
#'
#' @param se Column name of the standard errors.
#'
#' @return A data frame of processed GWAS summary statistics.
#'
#' @export
read_gwas <- function(gwas,
                      id = 'rsID',
                      A1 = 'ALT',
                      A2 = 'REF',
                      z = 'Z',
                      beta = 'ES',
                      se = 'SE'){

  # Extract relevant columns
  z_snp <- data.frame(id = gwas[,id], A1 = gwas[, A1], A2 = gwas[,A2])

  # Convert alleles to upper case
  z_snp$A1 <- toupper(z_snp$A1)
  z_snp$A2 <- toupper(z_snp$A2)

  if (z %in% colnames(gwas)){
    z_snp$z <- gwas[,z]
  } else {
    z_snp$z <- gwas[,beta]/gwas[,se]
    z_snp <- z_snp[!is.na(z_snp$z),]
  }

  z_snp <- z_snp[,c("id", "A1", "A2", "z")]

  return(z_snp)
}


#' @title Preprocess GWAS z-scores, harmonize GWAS z-scores with LD reference
#'
#' @param z_snp A data frame with two columns: "id", "A1", "A2", "z".
#' giving the z scores for snps. "A1" is effect allele. "A2" is the other allele.
#'
#' @param snp_map a list of SNP-to-region map for the reference.
#'
#' @param drop_multiallelic If TRUE, multiallelic variants will be dropped from the summary statistics.
#'
#' @param drop_strand_ambig If TRUE remove strand ambiguous variants (A/T, G/C).
#'
#' @param varID_converter_fun a user defined function to convert
#' GWAS variant IDs to the reference variant format.
#'
#' @param logfile The log filename. If NULL, print log info on screen.
#'
#' @return a data frame preprocessed GWAS z-scores
#'
#' @importFrom logging addHandler loginfo writeToFile
#' @importFrom data.table rbindlist
#'
#' @export
#'
preprocess_z_snp <- function(z_snp,
                             snp_map,
                             drop_multiallelic = TRUE,
                             drop_strand_ambig = TRUE,
                             varID_converter_fun = NULL,
                             logfile = NULL){

  if (!is.null(logfile)) {
    addHandler(writeToFile, file = logfile, level = "DEBUG")
  }

  loginfo("Preprocessing z_snp...")

  if (!inherits(snp_map,"list")){
    stop("'snp_map' should be a list object.")
  }

  snp_info <- as.data.frame(rbindlist(snp_map, idcol = "region_id"))

  loginfo("z_snp has %d variants in total", length(z_snp$id))

  # remove SNPs with missing values
  if (anyNA(z_snp)){
    missing.idx <- which(is.na(z_snp$z))
    loginfo("Remove variants with %d missing z-scores", length(missing.idx))
    z_snp <- z_snp[-missing.idx,]
  }

  # convert format of variant IDs
  if (!is.null(varID_converter_fun)){
    loginfo("Convert variant IDs")
    z_snp$id <- varID_converter_fun(z_snp$id)
  }

  # remove SNPs not in LD reference
  z_snp <- z_snp[z_snp$id %in% snp_info$id,]
  loginfo("%d variants left after filtering by the reference SNPs.", length(z_snp$id))

  # drop multiallelic variants (id not unique)
  if (drop_multiallelic) {
    duplicated.idx <- which(z_snp$id %in% z_snp$id[duplicated(z_snp$id)])
    if(length(duplicated.idx) > 0){
      loginfo("Remove %d multiallelic variants", length(duplicated.idx))
      z_snp <- z_snp[-duplicated.idx,]
    }
  }

  # harmonize alleles between z_snp and LD reference
  z_snp <- harmonize_z(z_snp, snp_info, drop_strand_ambig)

  if (length(z_snp$id) == 0){
    stop("No variants left after preprocessing and harmonization!")
  } else{
    loginfo("%d variants left after preprocessing and harmonization.", length(z_snp$id))
  }

  rownames(z_snp) <- NULL

  return(z_snp)
}




