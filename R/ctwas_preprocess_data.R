
#' Preprocess GWAS summary statistics, harmonize GWAS z-scores and LD reference, and filter LD mismatches using susie_rss
#'
#' @param z_snp A data frame with two columns: "id", "A1", "A2", "z". giving the z scores for
#' snps. "A1" is effect allele. "A2" is the other allele. If `harmonize= False`, A1 and A2 are not required.
#'
#' @param LD_R_dir a string, pointing to a directory containing all LD matrix files and variant information. Expects .RDS files which contain LD correlation matrices for a region/block.
#' For each RDS file, a file with same base name but ended with .Rvar needs to be present in the same folder. the .Rvar file has 5 required columns: "chrom", "id", "pos", "alt", "ref".
#' If using PredictDB format weights and \code{scale_by_ld_variance=T}, a 6th column is also required: "variance", which is the variance of the each SNP.
#' The order of rows needs to match the order of rows in .RDS file.
#'
#' @param chr a vector containing the chromosome numbers to process.
#'
#' @param ld_regions a string, population of LD reference, "EUR", "ASN", or "AFR".
#'
#' @param ld_regions_version a string, version LD reference, "b37", "b38".
#'
#' @param ld_regions_custom a string, path of LD regions file.
#'
#' @param filestem a string, filestem of reference LD matrix.
#'
#' @param gwas_n integer, GWAS sample size
#' @param outputdir a string, the directory to store output.
#'
#' @param outname a string, the output name.
#'
#' @param logfile the log file, if NULL will print log info on screen.
#'
#' @param drop_strand_ambig TRUE/FALSE. If TRUE, removes the ambiguous variant from the z scores.
#'
#' @param drop_multiallelic TRUE/FALSE. If TRUE, multiallelic variants will be dropped from the summary statistics.
#'
#' @param filter_ld_mismatch_susie TRUE/FALSE. If TRUE, detect LD mismatches using susie_rss,
#' flip the sign for variants detected with allele flipping, and filter problematic variants with p-value < 5e-8.
#'
#' @param ncore integer, number of cores for parallel computing.
#'
#' @importFrom logging addHandler loginfo
#' @importFrom tools file_ext
#' @importFrom foreach %dopar% foreach
#'
#' @export
#'
preprocess_z_ld <- function (z_snp,
                             ld_R_dir,
                             chr=NULL,
                             ld_regions = c("EUR", "ASN", "AFR"),
                             ld_regions_version = c("b37", "b38"),
                             ld_regions_custom = NULL,
                             filestem,
                             gwas_n,
                             outputdir = getwd(),
                             outname = NULL,
                             logfile = NULL,
                             drop_strand_ambig = F,
                             drop_multiallelic = T,
                             filter_ld_mismatch_susie = T,
                             ncore = 1){

  dir.create(outputdir, showWarnings = F, recursive=T)

  if (!is.null(logfile)) {
    addHandler(writeToFile, file = logfile, level = "DEBUG")
  }

  ld_Rfs <- write_ld_Rf(ld_R_dir, outname = outname, outputdir = outputdir)

  ld_snplist <- c()

  loginfo("GWAS summary stats (z_snp) has %d variants", length(z_snp$id))

  # drop multiallelic variants (id not unique)
  if (drop_multiallelic){
    duplicated.idx <- which(z_snp$id %in% z_snp$id[duplicated(z_snp$id)])
    if(length(duplicated.idx) > 0){
      loginfo("Remove %d multiallelic variants", length(duplicated.idx))
      z_snp <- z_snp[-duplicated.idx,]
    }
  }

  if(is.null(chr)){
    chr <- 1:22
  }else{
    outname <- paste0(outname, ".chr", chr)
  }
  loginfo("Process summary statistics for chromosomes: %s", chr)

  for (b in chr){
    loginfo("Harmonizing summary statistics for chromosome %s", b)

    ld_Rf <- ld_Rfs[b]
    ld_Rinfo <- data.table::fread(ld_Rf, header = T)
    ld_snpinfo <- read_ld_Rvar(ld_Rf)

    chrom <- unique(ld_snpinfo$chrom)
    if (length(chrom) > 1) {
      stop("Input LD reference not split by chromosome")
    }
    ld_snplist <- c(ld_snplist, ld_snpinfo$id) #store names of snps in ld reference

    z_snp <- match_allele_z_ld(z_snp, ld_snpinfo, drop_strand_ambig = drop_strand_ambig)
  }

  z_snp <- z_snp[z_snp$id %in% ld_snplist,]
  loginfo("%d variants left after allele matching", length(z_snp$id))

  if(filter_ld_mismatch_susie) {
    # select LD regions to run
    regions_df <- load_ld_regions(ld_regions = ld_regions,
                                  ld_regions_version = ld_regions_version,
                                  ld_regions_custom = ld_regions_custom)
    head(regions_df)
    regions_df <- regions_df[regions_df$chr %in% paste0("chr", chr), ]
    loginfo("Detect LD mismatches in %d LD regions", nrow(regions_df))

    # detect LD mismatches using susie_rss
    condz_dist <- detect_ld_mismatch_susie_rss(z_snp,
                                               regions_df,
                                               ld_R_dir = ld_R_dir,
                                               filestem = filestem,
                                               gwas_n = gwas_n,
                                               ncore = ncore)

    condz_dist.df <- do.call(rbind.data.frame, condz_dist)

    # flip and filter variants with LD mismatches
    detected_snps <- condz_dist.df$id[which(condz_dist.df$p_diff < 5e-8)]
    flipped_snps <- condz_dist.df$id[which(condz_dist.df$logLR > 2 & abs(condz_dist.df$z) > 2)]
    loginfo("%d variants (%.3f%%) with SuSiE_RSS diagnosis p-value < 5e-8", length(detected_snps),
                length(detected_snps)/length(condz_dist.df$id)*100)
    loginfo("%d variants with allele flipped.\n", length(flipped_snps))
    loginfo("%d variants before SuSiE_RSS filtering \n", length(z_snp$id))
    flip.idx <- which(z_snp$id %in% flipped_snps)
    filter.idx <- which(z_snp$id %in% setdiff(detected_snps, flipped_snps))

    loginfo("Flip %d variants detected with allele flipping.\n", length(flip.idx))
    z_snp$z[flip.idx] <- -z_snp$z[flip.idx]
    loginfo("Remove %d variants \n", length(filter.idx))
    z_snp <- z_snp[-filter.idx, ]
    loginfo("%d variants left after SuSiE_RSS filtering \n", length(z_snp$id))

    res <- list(z_snp = z_snp, regions = regions_df, condz_dist = condz_dist)
  }else{
    res <- list(z_snp = z_snp)
  }

  loginfo("Save processed result to %s", outputdir)
  saveRDS(res, file.path(outputdir, paste0(outname, ".res.RDS")))

  return(res)
}


#' Match alleles between GWAS and LD reference genotypes.
#' Flip signs when reverse complement matches.
#'
#' @param z_snp a data frame, with columns "id", "A1", "A2" and "z".
#'     Z scores for every SNP. "A1" is the effect allele.
#'
#' @param ld_snpinfo a data frame, snp info for LD reference,
#'  with columns "chrom", "id", "pos", "alt", "ref".
#'
#' @param drop_strand_ambig If TRUE, removes the ambiguous variant from the z scores.
#'
#' @return a data frame, z_snp with the "z" columns flipped to match LD ref.
#'
match_allele_z_ld <- function(z_snp, ld_snpinfo, drop_strand_ambig = F){
  # Remove variants in GWAS, but not in LD reference
  snpnames <- intersect(z_snp$id, ld_snpinfo$id)
  loginfo("%d variants in both GWAS and LD reference", length(snpnames))

  if (length(snpnames) != 0) {
    z.idx <- match(snpnames, z_snp$id)
    ld.idx <- match(snpnames, ld_snpinfo$id)
    # Allele flipping of unambiguous variants by checking ref/alt. alleles.
    qc <- allele.qc(z_snp[z.idx,]$A1, z_snp[z.idx,]$A2, ld_snpinfo[ld.idx,]$alt, ld_snpinfo[ld.idx,]$ref)
    ifflip <- qc[["flip"]]
    ifremove <- !qc[["keep"]]
    flip.idx <- z.idx[ifflip]
    loginfo("Flip alleles and z-scores for %d (%.2f%%) unambiguous variants", length(flip.idx), length(flip.idx)/length(z.idx)*100)
    z_snp[flip.idx, c("A1", "A2")] <- z_snp[flip.idx, c("A2", "A1")]
    z_snp[flip.idx, "z"] <- -z_snp[flip.idx, "z"]

    if (drop_strand_ambig & any(ifremove)) {
      # Remove strand ambiguous variants
      remove.idx <- z.idx[ifremove]
      loginfo("Remove %d strand ambiguous variants", length(remove.idx))
      z_snp <- z_snp[-remove.idx, , drop = F]
    }

  }
  return(z_snp)
}

# Load LD Regions (ldetect blocks)
load_ld_regions <- function(ld_regions = c("EUR", "ASN", "AFR"),
                            ld_regions_version = c("b37", "b38"),
                            ld_regions_custom = NULL) {

  if (is.null(ld_regions_custom)){
    ld_regions <- match.arg(ld_regions)
    ld_regions_version <- match.arg(ld_regions_version)
    regionfile <- system.file("extdata", "ldetect",
                              paste0(ld_regions, "." , ld_regions_version, ".bed"),
                              package="ctwas")
  } else {
    regionfile <- ld_regions_custom
  }

  loginfo("Load LD region from: %s", regionfile)
  regions_df <- read.table(regionfile, header = T)
  regions_df$locusID <- 1:nrow(regions_df)

  return(regions_df)
}


# Load reference LD matrix and SNP info
load_R_snp_info <- function(region_df, ld_R_dir, filestem = "ukb_b38_0.1"){
  filename <- sprintf("%s_chr%s.R_snp.%d_%d", filestem,
                      gsub("chr", "", region_df$chr), region_df$start, region_df$stop)
  if(!file.exists(file.path(ld_R_dir, paste0(filename, ".RDS"))) || !file.exists(file.path(ld_R_dir, paste0(filename, ".Rvar")))){
    stop("LD Reference files not exist!")
  }
  R_snp <- readRDS(file.path(ld_R_dir, paste0(filename, ".RDS")))
  R_snp_info <- read.table(file.path(ld_R_dir, paste0(filename, ".Rvar")), header=TRUE)
  return(list(R_snp = R_snp, R_snp_info = R_snp_info))
}

# Match GWAS sumstats with LD reference files. Only keep variants included in LD reference.
match_gwas_R_snp <- function(sumstats, R, snp_info){
  sumstats <- sumstats[sumstats$id %in% snp_info$id,]
  R_snp_index <- na.omit(match(sumstats$id, snp_info$id))
  sumstats$R_snp_index <- R_snp_index
  R <- R[R_snp_index, R_snp_index]
  stopifnot(nrow(sumstats) == nrow(R))
  return(list(sumstats = sumstats, R = R))
}

# Detect LD mismatches using SuSiE RSS
detect_ld_mismatch_susie_rss <- function (z_snp,
                                          regions_df,
                                          locusIDs=NULL,
                                          gwas_n,
                                          ld_R_dir = "/project2/mstephens/wcrouse/UKB_LDR_0.1",
                                          filestem = "ukb_b38_0.1",
                                          ncore = 1){

  if(is.null(locusIDs)){
    locusIDs <- regions_df$locusID
  }
  loginfo("Run LD mismatch diagnosis using susie_rss in %d loci", length(locusIDs))

  # registerDoParallel
  # cl <- parallel::makeCluster(ncore, outfile = "")
  # doParallel::registerDoParallel(cl)
  doParallel::registerDoParallel(cores=ncore)

  # condz_dist <- foreach(i=1:length(loci), .packages = "ctwas") %dopar%{
  condz_dist <- foreach(locusID = locusIDs) %dopar%{

    region_df <- regions_df[regions_df$locusID == locusID,]
    ldref_res <- load_R_snp_info(region_df, ld_R_dir, filestem)
    matched.z.ld <- match_gwas_R_snp(z_snp, ldref_res$R_snp, ldref_res$R_snp_info)
    z.locus <- matched.z.ld$sumstats
    R.locus <- matched.z.ld$R
    rm(matched.z.ld)

    # # Estimate lambda (consistency) between the z-scores and LD matrix
    # lambda <- estimate_s_rss(z = z.locus$z, R = R.locus, n = gwas_n)
    # loginfo("Estimated lambda between the z-scores and LD matrix: %.2f", lambda)

    # Compute expected z-scores based on conditional distribution of z-scores
    condz <- kriging_rss(z = z.locus$z, R = R.locus, n = gwas_n)
    condz$conditional_dist <- cbind(z.locus[,1:3], condz$conditional_dist)
    # compute p-values
    condz$conditional_dist$p_diff <- pchisq(condz$conditional_dist$z_std_diff^2, df = 1, lower.tail=FALSE)

    condz$conditional_dist
  }
  names(condz_dist) <- as.character(locusIDs)

  # parallel::stopCluster(cl)

  return(condz_dist)
}
