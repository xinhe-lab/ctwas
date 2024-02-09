
#' Preprocess GWAS summary statistics, harmonize GWAS z-scores and LD reference, and filter LD mismatches using susie_rss
#'
#' @param z_snp A data frame with two columns: "id", "A1", "A2", "z". giving the z scores for
#' snps. "A1" is effect allele. "A2" is the other allele.
#'
#' @param LD_R_dir a string, pointing to a directory containing all LD matrix files and variant information. Expects .RDS files which contain LD correlation matrices for a region/block.
#' For each RDS file, a file with same base name but ended with .Rvar needs to be present in the same folder. the .Rvar file has 5 required columns: "chrom", "id", "pos", "alt", "ref".
#' If using PredictDB format weights and \code{scale_by_ld_variance=T}, a 6th column is also required: "variance", which is the variance of the each SNP.
#' The order of rows needs to match the order of rows in .RDS file.
#'
#' @param chrom a vector containing the chromosome numbers to process.
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
#'
#' @param outputdir a string, the directory to store output.
#'
#' @param outname a string, the output name. It does not save output file if outname is NULL.
#'
#' @param logfile the log file, if NULL will print log info on screen.
#'
#' @param drop_multiallelic TRUE/FALSE. If TRUE, multiallelic variants will be dropped from the summary statistics.
#'
#' @param strand_ambig_action the action to take to harmonize strand ambiguous variants (A/T, G/C) between
#' the z scores and LD reference.
#' "drop" removes the ambiguous variant from the z scores.
#' "none" takes no additional action.
#' "recover" imputes the sign of ambiguous z scores using unambiguous z scores and the LD reference and
#' flips the z scores if there is a mismatch between the imputed sign and the observed sign of the z score.
#' This option is computationally intensive.
#'
#' @param detect_ld_mismatch TRUE/FALSE. If TRUE, detect LD mismatches by susie_rss,
#' and report problematic variants, and variants with allele flipping.
#'
#' @param flip_allele TRUE/FALSE. If TRUE, flip the sign for variants detected with allele flipping.
#'
#' @param remove_ld_mismatch TRUE/FALSE. If TRUE, remove problematic variants with LD mismatches.
#'
#' @param ncore integer, number of cores for parallel computing when detecting LD mismatches.
#'
#' @param save_result TRUE/FALSE. If TRUE, save results in the outputdir directory.
#'
#' @importFrom logging addHandler loginfo writeToFile
#' @importFrom tools file_ext
#' @importFrom foreach %dopar% foreach
#'
#' @export
#'
preprocess_z_ld <- function (z_snp,
                             ld_R_dir,
                             chrom=1:22,
                             ld_regions = c("EUR", "ASN", "AFR"),
                             ld_regions_version = c("b37", "b38"),
                             ld_regions_custom = NULL,
                             filestem,
                             gwas_n,
                             outputdir = getwd(),
                             outname = NULL,
                             logfile = NULL,
                             drop_multiallelic = T,
                             strand_ambig_action = c("none", "drop", "recover"),
                             detect_ld_mismatch = F,
                             flip_allele = T,
                             remove_ld_mismatch = F,
                             ncore = 1,
                             save_result = F){

  if (!is.null(logfile)) {
    addHandler(writeToFile, file = logfile, level = "DEBUG")
  }

  strand_ambig_action <- match.arg(strand_ambig_action)

  dir.create(outputdir, showWarnings = F, recursive=T)

  ld_Rfs <- write_ld_Rf(ld_R_dir, outname = outname, outputdir = outputdir)

  ld_snplist <- c()

  loginfo("GWAS summary stats (z_snp) has %d variants", length(z_snp$id))

  # drop multiallelic variants (id not unique)
  if (isTRUE(drop_multiallelic)) {
    duplicated.idx <- which(z_snp$id %in% z_snp$id[duplicated(z_snp$id)])
    if(length(duplicated.idx) > 0){
      loginfo("Drop %d multiallelic variants", length(duplicated.idx))
      z_snp <- z_snp[-duplicated.idx,]
    }
  }

  loginfo("Process summary statistics for chromosomes: %s", chrom)

  for (b in chrom){
    loginfo("Harmonizing summary statistics for chromosome %s", b)

    ld_Rf <- ld_Rfs[b]
    ld_Rinfo <- data.table::fread(ld_Rf, header = T)
    ld_snpinfo <- read_ld_Rvar(ld_Rf)

    if (length(unique(ld_snpinfo$chrom)) > 1) {
      stop("Input LD reference not split by chromosome")
    }
    ld_snplist <- c(ld_snplist, ld_snpinfo$id) #store names of snps in ld reference

    # harmonize alleles between z_snp and LD reference
    z_snp <- harmonize_z_ld(z_snp, ld_snpinfo,
                            strand_ambig_action = strand_ambig_action,
                            ld_Rinfo = ld_Rinfo)
  }

  z_snp <- z_snp[z_snp$id %in% ld_snplist,]
  loginfo("%d variants left after allele harmonization.", length(z_snp$id))

  if(length(z_snp$id) == 0)
    stop("No variants left!")

  if( isTRUE(detect_ld_mismatch) ) {
    # select LD regions to run
    regions_df <- load_ld_regions(ld_regions = ld_regions,
                                  ld_regions_version = ld_regions_version,
                                  ld_regions_custom = ld_regions_custom)
    regions_df <- regions_df[regions_df$chr %in% paste0("chr", chrom), ]
    loginfo("Detect LD mismatches in %d LD regions", nrow(regions_df))

    # detect LD mismatches using susie_rss
    ld_mismatch_res <- detect_ld_mismatch_susie_rss(z_snp,
                                                    regions_df,
                                                    ld_R_dir = ld_R_dir,
                                                    filestem = filestem,
                                                    gwas_n = gwas_n,
                                                    ncore = ncore)

    condz_dist <- ld_mismatch_res$condz_dist
    problematic_snps <- ld_mismatch_res$problematic_snps
    flipped_snps <- ld_mismatch_res$flipped_snps

    if (isTRUE(flip_allele)) {
      flip.idx <- which(z_snp$id %in% flipped_snps)
      loginfo("Flip %d variants.", length(flip.idx))
      z_snp$z[flip.idx] <- -z_snp$z[flip.idx]
      problematic_snps <- setdiff(problematic_snps, flipped_snps)
    }

    loginfo("Number of problematic variants: %d", length(problematic_snps))

    if (isTRUE(remove_ld_mismatch)) {
      filter.idx <- which(z_snp$id %in% problematic_snps)
      loginfo("Remove %d variants with LD mismatches.", length(filter.idx))
      z_snp <- z_snp[-filter.idx, ]
      loginfo("%d variants left after LD mismatch filtering.", length(z_snp$id))
    }

    ld_mismatch_res <- list(condz_dist = condz_dist,
                            problematic_snps = problematic_snps,
                            flipped_snps = flipped_snps)

    if (isTRUE(save_result)) {
      saveRDS(ld_mismatch_res, file.path(outputdir, paste0(outname, ".ld_mismatch_res.RDS")))
    }
  }else{
    ld_mismatch_res <- NULL
  }

  if (isTRUE(save_result)) {
    saveRDS(z_snp, file.path(outputdir, paste0(outname, ".harmonized.z_snp.RDS")))
  }

  return(list(z_snp = z_snp, ld_mismatch_res = ld_mismatch_res))
}


#' Detect LD mismatches using SuSiE RSS
#'
#' @param z_snp A data frame with two columns: "id", "A1", "A2", "z". giving the z scores for
#' snps. "A1" is effect allele. "A2" is the other allele.
#' @param regions_df A data frame with regions to be analyzed
#' @param LD_R_dir a string, pointing to a directory containing all LD matrix files and variant information. Expects .RDS files which contain LD correlation matrices for a region/block.
#' For each RDS file, a file with same base name but ended with .Rvar needs to be present in the same folder. the .Rvar file has 5 required columns: "chrom", "id", "pos", "alt", "ref".
#' If using PredictDB format weights and \code{scale_by_ld_variance=T}, a 6th column is also required: "variance", which is the variance of the each SNP.
#' The order of rows needs to match the order of rows in .RDS file.
#' @param filestem a string, filestem of reference LD matrix.
#' @param gwas_n integer, GWAS sample size
#' @param ncore integer, number of cores for parallel computing.
#'
#' @importFrom logging addHandler loginfo
#' @importFrom tools file_ext
#' @importFrom foreach %dopar% foreach
#'
#' @export
#'
detect_ld_mismatch_susie_rss <- function (z_snp,
                                          regions_df,
                                          ld_R_dir = NULL,
                                          filestem = NULL,
                                          gwas_n = NULL,
                                          ncore = 1){

  locusIDs <- regions_df$locusID
  loginfo("Run LD mismatch diagnosis in %d loci", length(locusIDs))

  nregions <- length(locusIDs)
  corelist <- lapply(1:ncore, function(core){
    njobs <- ceiling(nregions/ncore);
    jobs <- ((core-1)*njobs+1):(core*njobs);
    jobs[jobs<=nregions]
  })
  names(corelist) <- 1:ncore

  cl <- parallel::makeCluster(ncore, outfile = "")
  doParallel::registerDoParallel(cl)

  outlist <- foreach(core = 1:ncore, .combine = "c", .packages = c("ctwas", "stats")) %dopar% {

    locuslist_core <- locusIDs[corelist[[core]]]

    outlist_core <- list()
    for(locusID in locuslist_core) {

      # Load reference LD matrix and SNP info in the locus
      region_df <- regions_df[regions_df$locusID == locusID,]
      filename <- sprintf("%s_chr%s.R_snp.%d_%d", filestem,
                          gsub("chr", "", region_df$chr), region_df$start, region_df$stop)
      if(!file.exists(file.path(ld_R_dir, paste0(filename, ".RDS"))) || !file.exists(file.path(ld_R_dir, paste0(filename, ".Rvar")))){
        stop("LD Reference files not exist!")
      }
      R_snp <- readRDS(file.path(ld_R_dir, paste0(filename, ".RDS")))
      R_snp_info <- read.table(file.path(ld_R_dir, paste0(filename, ".Rvar")), header=T)

      # Match GWAS sumstats with LD reference files. Only keep variants included in LD reference.
      z.locus <- z_snp[z_snp$id %in% R_snp_info$id,]
      R_snp_index <- na.omit(match(z.locus$id, R_snp_info$id))
      z.locus$R_snp_index <- R_snp_index
      R.locus <- R_snp[R_snp_index, R_snp_index]
      stopifnot(nrow(z.locus) == nrow(R.locus))

      # # Estimate lambda (consistency) between the z-scores and LD matrix
      # lambda <- estimate_s_rss(z = z.locus$z, R = R.locus, n = gwas_n)

      # Compute expected z-scores based on conditional distribution of z-scores
      condz_dist <- kriging_rss(z = z.locus$z, R = R.locus, n = gwas_n)$conditional_dist
      condz_dist <- cbind(z.locus[,c("id", "A1", "A2")], condz_dist)
      # compute p-values
      condz_dist$p_diff <- pchisq(condz_dist$z_std_diff^2, df = 1, lower.tail=F)

      outlist_core[[as.character(locusID)]] <- condz_dist
    }
    outlist_core
  }
  parallel::stopCluster(cl)
  stopifnot(length(outlist) == length(locusIDs))

  # return problematic variants and flipped variants
  condz_dist <- data.table::rbindlist(outlist, idcol = "locus")
  problematic_snps <- condz_dist$id[which(condz_dist$p_diff < 5e-8)]
  flipped_snps <- condz_dist$id[which(condz_dist$logLR > 2 & abs(condz_dist$z) > 2)]

  return(list(condz_dist = condz_dist,
              problematic_snps = problematic_snps,
              flipped_snps = flipped_snps))
}

#' Preprocess PredictDB weights and harmonize with LD reference
#' (adapted from preharmonize_wgt_ld)
#'
#' @param weight a string, pointing to a directory with the fusion/twas format of weights, or a .db file in predictdb format.
#' A vector of multiple sets of weights in PredictDB format can also be specified; genes will have their filename appended
#' to their gene name to ensure IDs are unique.
#'
#' @param LD_R_dir a string, pointing to a directory containing all LD matrix files and variant information. Expects .RDS files which contain LD correlation matrices for a region/block.
#' For each RDS file, a file with same base name but ended with .Rvar needs to be present in the same folder. the .Rvar file has 5 required columns: "chrom", "id", "pos", "alt", "ref".
#' If using PredictDB format weights and \code{scale_by_ld_variance=T}, a 6th column is also required: "variance", which is the variance of the each SNP.
#' The order of rows needs to match the order of rows in .RDS file.
#'
#' @param outputdir a string, the directory to store output
#'
#' @param outname a string, the output name.
#'
#' @param strand_ambig_action the action to take to harmonize strand ambiguous variants (A/T, G/C) between
#' the weights and LD reference. "drop" removes the ambiguous variant from the prediction models. "none" treats the variant
#' as unambiguous, flipping the weights to match the LD reference and then taking no additional action. "recover" uses a procedure
#' to recover strand ambiguous variants. This procedure compares correlations between variants in the
#' LD reference and prediction models, and it can only be used with PredictDB format prediction models, which include this
#' information.
#'
#' @importFrom logging addHandler loginfo
#' @importFrom tools file_ext
#'
#' @export
#'
preprocess_wgt_ld <- function (weight,
                               ld_R_dir,
                               outputdir = getwd(),
                               outname,
                               strand_ambig_action = c("drop", "none", "recover")){

  strand_ambig_action <- match.arg(strand_ambig_action)

  dir.create(outputdir, showWarnings = F)

  # read the PredictDB weights
  sqlite <- RSQLite::dbDriver("SQLite")
  db = RSQLite::dbConnect(sqlite, weight)
  query <- function(...) RSQLite::dbGetQuery(db, ...)
  weight_table <- query("select * from weights")
  extra_table <- query("select * from extra")
  RSQLite::dbDisconnect(db)

  gnames <- unique(weight_table$gene)
  loginfo("Number of genes with weights provided: %s", length(gnames))

  # load LD information and subset to variants in weights
  ld_Rfs <- write_ld_Rf(ld_R_dir, outname = outname, outputdir = outputdir)

  ld_Rinfo <- lapply(ld_Rfs, data.table::fread, header=T)
  ld_Rinfo <- as.data.frame(do.call(rbind, ld_Rinfo))

  ld_snpinfo <- lapply(ld_Rfs, ctwas:::read_ld_Rvar)
  ld_snpinfo <- as.data.frame(do.call(rbind, ld_snpinfo))

  # remove variants in weight table, but not in LD reference
  loginfo("Number of variants in weights: %s", length(unique(weight_table$rsid)))
  loginfo("Remove %s variants in weights but not in LD reference", length(setdiff(weight_table$rsid, ld_snpinfo$id)))
  weight_table <- weight_table[weight_table$rsid %in% ld_snpinfo$id, ]

  # remove genes with no variants in LD reference
  loginfo("Remove %s genes with no variants in LD reference", length(setdiff(gnames, weight_table$gene)))
  gnames <- unique(weight_table$gene)
  loginfo("Number of genes left after removing variants not in LD reference: %s", length(gnames))

  # subset to variants in weight table
  ld_snpinfo <- ld_snpinfo[ld_snpinfo$id %in% weight_table$rsid,]

  if (strand_ambig_action=="recover"){
    # load covariances for variants in each gene (accompanies .db file)
    R_wgt_all <- read.table(gzfile(paste0(file_path_sans_ext(weight), ".txt.gz")), header=T)
    loginfo("Harmonizing strand ambiguous weights using correlations with unambiguous variants")
  }

  weight_table_harmonized <- list()

  loginfo("Processing weights for %s genes to match LD reference", length(gnames))

  for (i in 1:length(gnames)){

    if (i %% 1000 == 0){
      loginfo("Current gene: %s", i)
    }

    gname <- gnames[i]

    wgt <- weight_table[weight_table$gene==gname,]
    wgt.matrix <- as.matrix(wgt[, "weight", drop = F])

    rsid_varID <- wgt[,c("rsid", "varID")]

    rownames(wgt.matrix) <- wgt$rsid
    chrpos <- do.call(rbind, strsplit(wgt$varID, "_"))

    snps <- data.frame(gsub("chr", "", chrpos[, 1]), wgt$rsid,
                       "0", chrpos[, 2], wgt$eff_allele, wgt$ref_allele,
                       stringsAsFactors = F)
    colnames(snps) <- c("chrom", "id", "cm", "pos", "alt", "ref")
    snps$chrom <- as.integer(snps$chrom)
    snps$pos <- as.integer(snps$pos)

    chrom <- unique(snps$chrom)
    ld_Rinfo_chrom <- ld_Rinfo[ld_Rinfo$chrom==chrom,]

    if (strand_ambig_action=="recover"){
      #subset R_wgt_all to current weight
      R_wgt <- R_wgt_all[R_wgt_all$GENE == gname,]

      #convert covariance to correlation
      R_wgt_stdev <- R_wgt[R_wgt$RSID1==R_wgt$RSID2,]
      R_wgt_stdev <- setNames(sqrt(R_wgt_stdev$VALUE), R_wgt_stdev$RSID1)
      R_wgt$VALUE <- R_wgt$VALUE/(R_wgt_stdev[R_wgt$RSID1]*R_wgt_stdev[R_wgt$RSID2])

      #discard variances
      R_wgt <- R_wgt[R_wgt$RSID1!=R_wgt$RSID2,]

      #fix edge case where variance=0; treat correlations with these variants as uninformative (=0) for harmonization
      R_wgt$VALUE[is.nan(R_wgt$VALUE)] <- 0
    } else {
      R_wgt <- NULL
    }

    w <- harmonize_wgt_ld(wgt.matrix,
                          snps,
                          ld_snpinfo,
                          strand_ambig_action=strand_ambig_action,
                          ld_Rinfo=ld_Rinfo_chrom,
                          R_wgt=R_wgt,
                          wgt=wgt)

    wgt.matrix <- w[["wgt"]]
    snps <- w[["snps"]]

    wgt.matrix <- wgt.matrix[abs(wgt.matrix[, "weight"]) > 0, , drop = F]
    wgt.matrix <- wgt.matrix[complete.cases(wgt.matrix),, drop = F]

    snpnames <- intersect(rownames(wgt.matrix), ld_snpinfo$id)

    wgt.idx <- match(snpnames, rownames(wgt.matrix))
    wgt <- wgt.matrix[wgt.idx, "weight", drop = F]

    snps.idx <- match(snpnames, snps$id)
    snps <- snps[snps.idx,]

    if (length(snpnames)>0){
      weight_table_current <- data.frame(gene=gname,
                                         rsid=snps$id,
                                         varID=rsid_varID$varID[match(snps$id, rsid_varID$rsid)],
                                         ref_allele=snps$ref,
                                         eff_allele=snps$alt,
                                         weight=wgt[,"weight"])

      weight_table_harmonized[[gname]] <- weight_table_current
    }
  }
  weight_table_harmonized <- do.call(rbind, weight_table_harmonized)

  loginfo("Number of genes with weights harmonized with LD reference: %s", length(unique(weight_table_harmonized$gene)))

  extra_table <- extra_table[extra_table$gene %in% weight_table_harmonized$gene,]

  if (file.exists(paste0(outputdir, outname, ".db"))){
    invisible(file.remove(paste0(outputdir, outname, ".db")))
  }

  db <- RSQLite::dbConnect(sqlite, file.path(outputdir, paste0(outname, ".db")))
  RSQLite::dbWriteTable(db, "extra", extra_table)
  RSQLite::dbWriteTable(db, "weights", weight_table_harmonized)
  RSQLite::dbDisconnect(db)

  invisible(file.remove(file.path(outputdir, paste0(outname, "_ld_R_chr", 1:22, ".txt"))))

  loginfo("Save processed weights to: %s", file.path(outputdir, paste0(outname, ".db")))

  return(list(weight_table = weight_table_harmonized, extra_table = extra_table))
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

  # loginfo("Load LD region from: %s", regionfile)
  regions_df <- read.table(regionfile, header = T)
  regions_df$locusID <- 1:nrow(regions_df)

  return(regions_df)
}
