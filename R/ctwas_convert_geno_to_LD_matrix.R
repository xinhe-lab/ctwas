
#' Convert PLINK genotype data to LD matrices and SNP info files,
#' save LD matrices as .RDS files and SNP info as .Rvar files
#'
#' @param region_info a data frame of region information, with columns: chr, start, stop positions
#'
#' @param genotype_files Reference genotype files in PLINK binary genotype data in .pgen or .bed format
#'
#' @param varinfo_files Reference variant information files in PLINK .pvar or .bim format.
#'
#' The output will use the genome positions in \code{varinfo_files}.
#'
#' @param chrom a vector of chromosomes to process genotype data
#'
#' @param outputdir Output directory
#'
#' @param outname Output filestem
#'
#' @param include_variance TRUE/FALSE, if TRUE, include variance in .Rvar output
#'
#' @param include_allele_freq TRUE/FALSE, if TRUE, include allele frequency in .Rvar output
#'
#' @importFrom logging loginfo
#'
#' @return a data frame of updated region info, with paths of LD matrices and variance information files
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
                                      include_allele_freq = TRUE) {

  loginfo("Convert genotype data to LD matrices")

  if (!dir.exists(outputdir)){
    dir.create(outputdir, showWarnings=FALSE, recursive = TRUE)
  }

  if (is.character(region_info$chr)) {
    region_info$chr <- readr::parse_number(region_info$chr)
  }
  region_info <- region_info[order(region_info$chr, region_info$start), ]
  region_info$region_id <- paste0(region_info$chr, ":", region_info$start, "-", region_info$stop)

  # load SNP positions from varinfo_files
  ld_snpinfo_all <- do.call(rbind, lapply(varinfo_files, read_var_info))

  # read genotype files and prepare pvar files accompanying the genotype table
  pvar_files <- sapply(genotype_files, prep_pvar, outputdir = outputdir)

  for (b in chrom){
    region_info_chr <- region_info[region_info$chr == b, ]
    loginfo("No. regions in chr%s: %d", b, nrow(region_info_chr))

    ld_snpinfo_chr <- ld_snpinfo_all[ld_snpinfo_all$chr == b,]

    pgen <- prep_pgen(genotype_files[b], pvar_files[b])
    snpinfo <- read_pvar(pvar_files[b])

    if (unique(snpinfo$chrom) != b) {
      stop("Input genotype file not split by chromosome or not in correct order")
    }

    pb <- txtProgressBar(min = 0, max = nrow(region_info_chr), initial = 0, style = 3)

    for (rn in 1:nrow(region_info_chr)) {

      # loginfo("Region %s", rn)
      regioninfo <- region_info_chr[rn, ]
      region_id <- paste0(regioninfo$chr, ":", regioninfo$start, "-", regioninfo$stop)
      region_idx <- which(region_info$region_id == region_id)
      region_start <- regioninfo$start
      region_stop <- regioninfo$stop

      LD_matrix_file <- file.path(outputdir, paste0(outname, "_chr", b, ".R_snp.", region_start, "_", region_stop, ".RDS"))
      LD_Rvar_file <- file.path(outputdir, paste0(outname, "_chr", b, ".R_snp.", region_start, "_", region_stop, ".Rvar"))

      outfile_temp <- paste0(LD_Rvar_file, "-temp")

      if (!(file.exists(LD_matrix_file)) | !(file.exists(LD_Rvar_file)) | file.exists(outfile_temp)){

        file.create(outfile_temp)

        # use the position and allele information from LD reference SNP info
        sidx_ldref <- which(ld_snpinfo_chr$pos >= region_start & ld_snpinfo_chr$pos < region_stop)
        sid_ldref <- ld_snpinfo_chr$id[sidx_ldref]

        sidx <- which(snpinfo$id %in% sid_ldref)
        sid <- snpinfo[sidx, "id"]

        X.g <- read_pgen(pgen, variantidx = sidx)
        R_snp <- Rfast::cora(X.g)
        rownames(R_snp) <- sid
        colnames(R_snp) <- sid
        R_snp <- R_snp[sid_ldref, sid_ldref]
        rownames(R_snp) <- NULL
        colnames(R_snp) <- NULL
        saveRDS(R_snp, file=LD_matrix_file)

        R_var <- ld_snpinfo_chr[sidx_ldref, c("chr","id","pos","alt","ref")]
        colnames(R_var) <- c("chrom","id","pos","alt","ref")

        if (isTRUE(include_variance)) {
          R_snp_variances <- Rfast::colVars(X.g)
          names(R_snp_variances) <- sid
          R_snp_variances <- R_snp_variances[sid_ldref]
          names(R_snp_variances) <- NULL
          R_var <- cbind(R_var, variance = R_snp_variances)
        }

        if (isTRUE(include_allele_freq)){
          R_snp_AFs <- colSums(X.g)/(2*nrow(X.g))
          names(R_snp_AFs) <- sid
          R_snp_AFs <- R_snp_AFs[sid_ldref]
          names(R_snp_AFs) <- NULL
          R_var <- cbind(R_var, allele_freq = R_snp_AFs)
        }

        write.table(R_var, file=LD_Rvar_file, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

        file.remove(outfile_temp)
      }

      region_info[region_idx, "LD_matrix"] <- LD_matrix_file
      region_info[region_idx, "SNP_info"] <- LD_Rvar_file

      setTxtProgressBar(pb, rn)
    }

    close(pb)

  }

  file.remove(pvar_files)

  return(region_info)
}
