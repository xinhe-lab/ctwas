
#' @title Preprocess GWAS z-scores, harmonize GWAS z-scores with LD reference
#'
#' @param z_snp A data frame with two columns: "id", "A1", "A2", "z". giving the z scores for
#' snps. "A1" is effect allele. "A2" is the other allele.
#'
#' @param snp_map a list of SNP-to-region map for the reference.
#'
#' @param drop_multiallelic If TRUE, multiallelic variants will be dropped from the summary statistics.
#'
#' @param drop_strand_ambig If TRUE remove strand ambiguous variants (A/T, G/C).
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

  # remove SNPs not in LD reference
  z_snp <- z_snp[z_snp$id %in% snp_info$id,]
  loginfo("%d variants left after filtering by snp_map", length(z_snp$id))

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

