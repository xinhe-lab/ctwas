
#' @title Converts PLINK genotype data to LD matrices and SNP info files,
#' saves LD matrices as .RDS files and SNP info as .Rvar files
#'
#' @param region_info a data frame of region definitions, with columns:
#' chrom, start, stop, and region_id.
#'
#' @param genotype_files Reference genotype files in PLINK binary genotype data
#' in .pgen or .bed format. It should contain files for all chromosomes
#' (from 1 to 22), one file per chromosome.
#'
#' @param varinfo_files Reference variant information files in PLINK
#' .pvar or .bim format. It could have one file per chromosome or
#' have one big file for all chromosomes.
#' The output will use the genome positions in \code{varinfo_files}.
#'
#' @param chrom a vector of chromosome numbers to process genotype data.
#'
#' @param outputdir Output directory.
#'
#' @param outname Output filestem.
#'
#' @param include_variance If TRUE, include variance in .Rvar output.
#'
#' @param include_allele_freq If TRUE, include allele frequency in .Rvar output.
#'
#' @param show_progress_bar If TRUE, print progress bar.
#'
#' @param verbose If TRUE, print detail messages.
#'
#' @param logfile The log filename. If NULL, print log info on screen.
#'
#' @return a data frame of region_metatable, with region definitions and
#' filenames of LD matrices and variant information.
#'
#' @importFrom logging addHandler loginfo writeToFile
#' @importFrom utils write.table
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom readr parse_number
#' @importFrom data.table fwrite
#'
#' @export
#'
convert_geno_to_LD_matrix <- function(region_info,
                                      genotype_files,
                                      varinfo_files,
                                      chrom = 1:22,
                                      outputdir = getwd(),
                                      outname = "",
                                      include_variance = TRUE,
                                      include_allele_freq = TRUE,
                                      show_progress_bar = TRUE,
                                      verbose = FALSE,
                                      logfile = NULL) {

  if (!is.null(logfile)){
    addHandler(writeToFile, file=logfile, level='DEBUG')
  }

  loginfo("Convert genotype data to LD matrices")

  if (!requireNamespace("Rfast", quietly = TRUE)){
    stop("Rfast package is required for converting genotype to LD matrix")
  }

  if (!dir.exists(outputdir)){
    dir.create(outputdir, showWarnings=FALSE, recursive = TRUE)
  }

  if (is.character(region_info$chrom)) {
    region_info$chrom <- parse_number(region_info$chrom)
  }

  region_info <- region_info
  region_info$start <- as.numeric(region_info$start)
  region_info$stop <- as.numeric(region_info$stop)
  region_info <- region_info[order(region_info$chrom, region_info$start), ]

  if (is.null(region_info$region_id)) {
    region_info$region_id <- paste(region_info$chrom, region_info$start, region_info$stop, sep = "_")
  }

  # load variant positions from varinfo_files
  loginfo("Load variant information ...")
  ref_snpinfo_all <- do.call(rbind, lapply(varinfo_files, read_var_info))

  # read genotype files and prepare pvar files accompanying the genotype table
  # pvar_files <- sapply(genotype_files, prep_pvar, outputdir = outputdir)

  region_metatable <- region_info
  region_metatable$LD_file <- NA
  region_metatable$SNP_file <- NA

  for (chr in chrom){

    region_info_chr <- region_info[region_info$chrom == chr, ]
    loginfo("No. regions in chr%s: %d", chr, nrow(region_info_chr))

    ref_snpinfo_chr <- ref_snpinfo_all[ref_snpinfo_all$chrom == chr,]

    # read genotype file and prepare pvar file accompanying the genotype table
    pgen_file <- genotype_files[chr]
    pvar_file <- prep_pvar(pgen_file, outputdir)
    pgen <- prep_pgen(pgen_file, pvar_file)
    snpinfo <- read_pvar(pvar_file)

    if (unique(snpinfo$chrom) != chr) {
      stop("Input genotype file not split by chromosome or not in correct order")
    }

    if (show_progress_bar) {
      pb <- txtProgressBar(min = 0, max = nrow(region_info_chr), initial = 0, style = 3)
    }

    for (i in 1:nrow(region_info_chr)) {
      regioninfo <- region_info_chr[i, ]
      region_id <- regioninfo$region_id

      region_start <- regioninfo$start
      region_stop <- regioninfo$stop

      LD_matrix_file <- file.path(outputdir,
                                  paste0(outname, "_chr", chr, ".R_snp.", region_start, "_", region_stop, ".RDS"))
      LD_Rvar_file <- file.path(outputdir,
                                paste0(outname, "_chr", chr, ".R_snp.", region_start, "_", region_stop, ".Rvar"))

      outfile_temp <- paste0(LD_Rvar_file, "-temp")

      if (!(file.exists(LD_matrix_file)) | !(file.exists(LD_Rvar_file)) | file.exists(outfile_temp)) {

        if (verbose) {
          loginfo("Processing region %s ...", region_id)
        }

        file.create(outfile_temp)

        # make LD matrix
        # use the snp positions and allele information in the region from the LD reference
        sidx_ref <- which(ref_snpinfo_chr$pos >= region_start & ref_snpinfo_chr$pos < region_stop)
        sid_ref <- ref_snpinfo_chr$id[sidx_ref]

        # select snp ids available in LD reference
        sid <- intersect(snpinfo$id, sid_ref)
        sidx <- match(sid, snpinfo$id)

        X.g <- read_pgen(pgen, variantidx = sidx)
        R_snp <- Rfast::cora(X.g)
        rownames(R_snp) <- sid
        colnames(R_snp) <- sid

        # convert to the idx in LD reference
        R_snp <- R_snp[sid_ref, sid_ref]
        rownames(R_snp) <- NULL
        colnames(R_snp) <- NULL
        saveRDS(R_snp, file=LD_matrix_file)

        # make Rvar file
        R_var <- ref_snpinfo_chr[sidx_ref, c("chrom","id","pos","alt","ref")]

        if (include_variance) {
          R_snp_variances <- Rfast::colVars(X.g)
          names(R_snp_variances) <- sid
          R_snp_variances <- R_snp_variances[sid_ref]
          names(R_snp_variances) <- NULL
          R_var <- cbind(R_var, variance = R_snp_variances)
        }

        if (include_allele_freq) {
          R_snp_AFs <- colSums(X.g)/(2*nrow(X.g))
          names(R_snp_AFs) <- sid
          R_snp_AFs <- R_snp_AFs[sid_ref]
          names(R_snp_AFs) <- NULL
          R_var <- cbind(R_var, allele_freq = R_snp_AFs)
        }

        fwrite(R_var, file=LD_Rvar_file, sep="\t", col.names=TRUE, row.names=FALSE)

        file.remove(outfile_temp)
      }

      region_idx <- which(region_metatable$region_id == region_id)
      region_metatable[region_idx, "LD_file"] <- LD_matrix_file
      region_metatable[region_idx, "SNP_file"] <- LD_Rvar_file

      if (show_progress_bar) {
        setTxtProgressBar(pb, i)
      }
    }

    if (show_progress_bar) {
      close(pb)
    }
  }

  return(region_metatable)
}
