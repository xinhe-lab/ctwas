
#' Preprocess GWAS z-scores, harmonize GWAS z-scores with LD reference
#'
#' @param z_snp A data frame with two columns: "id", "A1", "A2", "z". giving the z scores for
#' snps. "A1" is effect allele. "A2" is the other allele.
#'
#' @param region_info a data frame of region definition and associated file names.
#'
#' @param gwas_n integer, GWAS sample size
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
                             region_info,
                             gwas_n = NULL,
                             drop_multiallelic = TRUE,
                             drop_strand_ambig = TRUE,
                             logfile = NULL){

  if (!is.null(logfile)) {
    addHandler(writeToFile, file = logfile, level = "DEBUG")
  }

  loginfo("Preprocessing z_snp...")
  loginfo("z_snp has %d variants in total", length(z_snp$id))

  # remove SNPs not in LD reference
  ld_snpinfo <- do.call(rbind, lapply(region_info$SNP_info, read_LD_SNP_file))
  z_snp <- z_snp[z_snp$id %in% ld_snpinfo$id,]

  # drop multiallelic variants (id not unique)
  if (isTRUE(drop_multiallelic)) {
    duplicated.idx <- which(z_snp$id %in% z_snp$id[duplicated(z_snp$id)])
    if(length(duplicated.idx) > 0){
      loginfo("Drop %d multiallelic variants", length(duplicated.idx))
      z_snp <- z_snp[-duplicated.idx,]
    }
  }

  # harmonize alleles between z_snp and LD reference
  z_snp <- harmonize_z_ld(z_snp, ld_snpinfo, drop_strand_ambig)

  if (length(z_snp$id) == 0){
    stop("No variants left after preprocessing and harmonization!")
  } else{
    loginfo("%d variants left after preprocessing and harmonization.", length(z_snp$id))
  }

  return(z_snp)
}


#' Detect LD mismatches using SuSiE RSS
#'
#' @param z_snp A data frame with two columns: "id", "A1", "A2", "z". giving the z scores for
#' snps. "A1" is effect allele. "A2" is the other allele.
#' @param ld_Rinfo a vector of paths to the variant information for all LD matrices
#' @param gwas_n integer, GWAS sample size
#' @param ncore integer, number of cores for parallel computing.
#' @param p_diff_thresh numeric, p-value threshold for identifying problematic SNPs
#' with significant difference between observed z-scores and estimated values
#'
#' @importFrom logging addHandler loginfo
#' @importFrom foreach %dopar% foreach
#'
#' @export
#'
detect_ld_mismatch_susie <- function(z_snp,
                                     region_info,
                                     gwas_n = NULL,
                                     ncore = 1,
                                     p_diff_thresh = 5e-8){

  region_tags <- region_info$region_tag
  loginfo("Run LD mismatch diagnosis in %d regions", length(region_tags))

  nregions <- length(region_tags)
  corelist <- lapply(1:ncore, function(core){
    njobs <- ceiling(nregions/ncore);
    jobs <- ((core-1)*njobs+1):(core*njobs);
    jobs[jobs<=nregions]
  })
  names(corelist) <- 1:ncore

  cl <- parallel::makeCluster(ncore, outfile = "")
  doParallel::registerDoParallel(cl)

  outlist <- foreach(core = 1:ncore, .combine = "c", .packages = c("ctwas", "stats")) %dopar% {

    region_tags_core <- region_tags[corelist[[core]]]

    outlist_core <- list()
    for(region_tag in region_tags_core) {

      # Load reference LD matrix and SNP info in the region
      R_snp <- load_LD(region_info$LD_matrix[region_info$region_tag == region_tag])
      ld_snpinfo <- read_LD_SNP_file(region_info$SNP_info[region_info$region_tag == region_tag])

      # Match GWAS sumstats with LD reference files. Only keep variants included in LD reference.
      z_snp.region <- z_snp[z_snp$id %in% ld_snpinfo$id,]
      R_snp.idx <- match(z_snp.region$id, ld_snpinfo$id)
      R_snp.region <- R_snp[R_snp.idx, R_snp.idx]
      stopifnot(nrow(z_snp.region) == nrow(R_snp.region))

      # # Estimate lambda (consistency) between the z-scores and LD matrix
      # lambda <- estimate_s_rss(z = z.locus$z, R = R.locus, n = gwas_n)

      # Compute expected z-scores based on conditional distribution of z-scores
      condz_dist <- kriging_rss(z = z_snp.region$z, R = R_snp.region, n = gwas_n)$conditional_dist
      condz_dist <- cbind(z_snp.region[,c("id", "A1", "A2")], condz_dist)

      # compute p-values for the significance of z-score difference between observed and estimated values
      condz_dist$p_diff <- pchisq(condz_dist$z_std_diff^2, df = 1, lower.tail=F)

      outlist_core[[as.character(region_tag)]] <- condz_dist
    }
    outlist_core
  }
  parallel::stopCluster(cl)
  stopifnot(length(outlist) == length(region_tags))

  # return problematic variants and flipped variants
  condz_dist <- data.table::rbindlist(outlist, idcol = "region_tag")
  problematic_snps <- condz_dist$id[which(condz_dist$p_diff < p_diff_thresh)]
  flipped_snps <- condz_dist$id[which(condz_dist$logLR > 2 & abs(condz_dist$z) > 2)]

  return(list(condz_dist = condz_dist,
              problematic_snps = problematic_snps,
              flipped_snps = flipped_snps))
}

