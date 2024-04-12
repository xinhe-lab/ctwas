
#' Convert PLINK genotype data to LD matrices and SNP info files,
#' save LD matrices as .RDS files and SNP info as .Rvar files
#'
#' @param region_info a data frame of region information, with columns: chr, start, stop positions
#' @param genotype_files Reference genotype files in PLINK binary genotype format (BED or PGEN format)
#' @param snpinfo_files Reference variant information files for the genome version in PLINK BIM format
#' @param chrom a vector of chromosomes to run
#' @param outputdir Output directory
#' @param outname Output filestem
#' @param print_variance TRUE/FALSE, if TRUE, print variance
#' @param print_allele_freq TRUE/FALSE, if TRUE, print allele frequency
#'
#' @importFrom logging loginfo
#'
#' @export
#'
convert_geno_to_LD_matrix <- function(region_info,
                                      genotype_files,
                                      snpinfo_files,
                                      chrom = 1:22,
                                      outputdir = getwd(),
                                      outname = "LD",
                                      print_variance = TRUE,
                                      print_allele_freq = TRUE) {

  if (!dir.exists(outputdir)){
    dir.create(outputdir, showWarnings=FALSE, recursive = TRUE)
  }

  if (is.character(region_info$chr)) {
    region_info$chr <- readr::parse_number(region_info$chr)
  }
  region_info <- region_info[order(region_info$chr, region_info$start), ]

  # load positions
  ld_snpinfo_all <- do.call(rbind, lapply(snpinfo_files, data.table::fread))
  colnames(ld_snpinfo_all) <- c("chr", "id", "cm", "pos", "alt", "ref")

  ld_pvarfs <- sapply(genotype_files, prep_pvar, outputdir = outputdir)

  for (b in chrom){
    loginfo("Convert reference genotype data to LD matrices for chr%s", b)
    regioninfo <- region_info[region_info$chr == b, ]
    ld_snpinfo_chr <- ld_snpinfo_all[ld_snpinfo_all$chr == b,]

    ld_pvarf <- ld_pvarfs[b]
    ld_pgen <- prep_pgen(genotype_files[b], ld_pvarf)
    snpinfo <- read_pvar(ld_pvarf)

    if (unique(snpinfo$chrom) != b) {
      stop("Input genotype file not split by chromosome or not in correct order")
    }

    loginfo("No. regions in chr%s: %d", b, nrow(regioninfo))

    pb <- txtProgressBar(min = 0, max = nrow(regioninfo), initial = 0, style = 3)

    for (rn in 1:nrow(regioninfo)) {

      # loginfo("Region %s", rn)
      p0 <- regioninfo[rn, "start"]
      p1 <- regioninfo[rn, "stop"]

      LD_matrix_file <- file.path(outputdir, paste0(outname, "_chr", b, ".R_snp.", p0, "_", p1, ".RDS"))
      LD_Rvar_file <- file.path(outputdir, paste0(outname, "_chr", b, ".R_snp.", p0, "_", p1, ".Rvar"))
      outfile_temp <- paste0(LD_Rvar_file, "-temp")

      if (!(file.exists(LD_matrix_file)) | !(file.exists(LD_Rvar_file)) | file.exists(outfile_temp)){

        file.create(outfile_temp)

        # use the position and allele information from LD reference SNP info
        sidx_ldref <- which(ld_snpinfo_chr$pos >= p0 & ld_snpinfo_chr$pos < p1)
        sid_ldref <- ld_snpinfo_chr$id[sidx_ldref]

        sidx <- which(snpinfo$id %in% sid_ldref)
        sid <- snpinfo[sidx, "id"]

        X.g <- read_pgen(ld_pgen, variantidx = sidx)
        R_snp <- Rfast::cora(X.g)
        rownames(R_snp) <- sid
        colnames(R_snp) <- sid
        R_snp <- R_snp[sid_ldref, sid_ldref]
        rownames(R_snp) <- NULL
        colnames(R_snp) <- NULL
        saveRDS(R_snp, file=LD_matrix_file)

        R_var <- ld_snpinfo_chr[sidx_ldref, c("chr","id","pos","alt","ref")]
        colnames(R_var) <- c("chrom","id","pos","alt","ref")

        if (isTRUE(print_variance)) {
          R_snp_variances <- Rfast::colVars(X.g)
          names(R_snp_variances) <- sid
          R_snp_variances <- R_snp_variances[sid_ldref]
          names(R_snp_variances) <- NULL
          R_var <- cbind(R_var, variance = R_snp_variances)
        }

        if (isTRUE(print_allele_freq)){
          R_snp_AFs <- colSums(X.g)/(2*nrow(X.g))
          names(R_snp_AFs) <- sid
          R_snp_AFs <- R_snp_AFs[sid_ldref]
          names(R_snp_AFs) <- NULL
          R_var <- cbind(R_var, allele_freq = R_snp_AFs)
        }

        write.table(R_var, file=LD_Rvar_file, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

        file.remove(outfile_temp)
      }

      setTxtProgressBar(pb, rn)
    }

    close(pb)

  }

}
