
#' Preprocess GWAS z-scores, harmonize GWAS z-scores with LD reference
#'
#' @param z_snp A data frame with two columns: "id", "A1", "A2", "z". giving the z scores for
#' snps. "A1" is effect allele. "A2" is the other allele.
#'
#' @param snp_info a data frame, SNP info for LD reference,
#'  with columns "chrom", "id", "pos", "alt", "ref".
#'
#' @param drop_multiallelic TRUE/FALSE. If TRUE, multiallelic variants will be dropped from the summary statistics.
#'
#' @param drop_strand_ambig TRUE/FALSE, if TRUE remove strand ambiguous variants (A/T, G/C).
#'
#' @param logfile the log file, if NULL will print log info on screen.
#'
#' @importFrom logging addHandler loginfo writeToFile
#'
#' @return a data frame preprocessed GWAS z-scores
#'
#' @export
#'
preprocess_z_snp <- function(z_snp,
                             snp_info,
                             drop_multiallelic = TRUE,
                             drop_strand_ambig = TRUE,
                             logfile = NULL){

  if (!is.null(logfile)) {
    addHandler(writeToFile, file = logfile, level = "DEBUG")
  }

  # Check LD reference SNP info
  target_header <- c("chrom", "id", "pos", "alt", "ref")
  if (!all(target_header %in% colnames(snp_info))){
    stop("SNP info needs to contain the following columns: ",
         paste(target_header, collapse = " "))
  }

  loginfo("Preprocessing z_snp...")
  loginfo("z_snp has %d variants in total", length(z_snp$id))

  # remove SNPs not in LD reference
  z_snp <- z_snp[z_snp$id %in% snp_info$id,]

  # drop multiallelic variants (id not unique)
  if (drop_multiallelic) {
    duplicated.idx <- which(z_snp$id %in% z_snp$id[duplicated(z_snp$id)])
    if(length(duplicated.idx) > 0){
      loginfo("Drop %d multiallelic variants", length(duplicated.idx))
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

  return(z_snp)
}

