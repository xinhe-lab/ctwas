
#' Preprocess PredictDB weights and harmonize with LD reference
#'
#' @param weight_file a string, pointing to a directory with the fusion/twas format of weights, or a .db file in predictdb format.
#' A vector of multiple sets of weights in PredictDB format can also be specified; genes will have their filename appended
#' to their gene name to ensure IDs are unique.
#'
#' @param region_info a data frame of region definition and associated file names.
#'
#' @param z_snp A data frame with columns: "id", "A1", "A2", "z". giving the z scores for
#' snps. "A1" is effect allele. "A2" is the other allele. For harmonized data, A1 and A2 are not required.
#'
#' @param weight_format weight format: PredictDB, or Fusion
#'
#' @param filter_protein_coding_genes TRUE/FALSE. If TRUE, keep protein coding genes only. This option is only for PredictDB weights
#'
#' @param drop_strand_ambig TRUE/FALSE, if TRUE remove strand ambiguous variants (A/T, G/C).
#'
#' @param scale_by_ld_variance TRUE/FALSE. If TRUE, PredictDB weights are scaled by genotype variance, which is the default
#' behavior for PredictDB
#'
#' @return a list of processed weights
#'
#' @export
#'
preprocess_weights <- function(weight_file,
                               region_info,
                               z_snp,
                               weight_format = c("PredictDB", "Fusion"),
                               drop_strand_ambig = TRUE,
                               filter_protein_coding_genes = FALSE,
                               scale_by_ld_variance = TRUE){

  weight_format <- match.arg(weight_format)
  # load LD SNPs information
  region_info <- region_info[order(region_info$chrom, region_info$start),]
  ld_snpinfo <- read_LD_SNP_files(region_info$SNP_info)
  weights <- list()
  for(weight in weight_file){
    loginfo("Load weight: %s", weight)
    loaded_weight <- load_weights(weight, weight_format = weight_format, filter_protein_coding_genes = filter_protein_coding_genes)
    weight_table <- loaded_weight$weight_table
    weight_name <- loaded_weight$weight_name

    gnames <- unique(weight_table$gene)
    loginfo("Number of genes with weights provided: %d in %s", length(gnames), weight_name)
    # remove variants in weight table, but not in LD reference and GWAS
    loginfo("Number of variants in weights: %d", length(unique(weight_table$rsid)))
    # take the intersect of SNPs in weights, LD reference and z_snp
    snpnames <- Reduce(intersect, list(weight_table$rsid, ld_snpinfo$id, z_snp$id))
    # loginfo("Remove %d variants after intersecting with LD reference and GWAS", length(setdiff(weight_table$rsid, snpnames)))
    weight_table <- weight_table[weight_table$rsid %in% snpnames, ]
    # loginfo("Remove %s genes after intersecting with LD reference and GWAS", length(setdiff(gnames, weight_table$gene)))
    gnames <- unique(weight_table$gene)
    loginfo("%d variants and %d genes left after intersecting with LD reference and GWAS z_snp", length(snpnames), length(gnames))
    # subset to variants in weight table
    ld_snpinfo_wgt <- ld_snpinfo[ld_snpinfo$id %in% weight_table$rsid,]
    loginfo("Harmonize weights with LD reference")

    pb <- txtProgressBar(min = 0, max = length(gnames), initial = 0, style = 3)
    for (i in 1:length(gnames)){
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
      chrom <- unique(snps$chrom)
      snps$chrom <- as.integer(snps$chrom)
      snps$pos <- as.integer(snps$pos)
      w <- harmonize_wgt_ld(wgt.matrix,
                            snps,
                            ld_snpinfo_wgt,
                            drop_strand_ambig = drop_strand_ambig)
      wgt.matrix <- w[["wgt"]]
      snps <- w[["snps"]]
      wgt.matrix <- wgt.matrix[abs(wgt.matrix[, "weight"]) > 0, , drop = F]
      wgt.matrix <- wgt.matrix[complete.cases(wgt.matrix),, drop = F]

      snpnames <- intersect(rownames(wgt.matrix), ld_snpinfo_wgt$id)
      wgt.idx <- match(snpnames, rownames(wgt.matrix))
      wgt <- wgt.matrix[wgt.idx, "weight", drop = F]

      snps.idx <- match(snpnames, snps$id)
      snps <- snps[snps.idx,]

      if (isTRUE(scale_by_ld_variance)){
        ld_snpinfo_wgt.idx <- match(snpnames, ld_snpinfo_wgt$id)
        wgt <- wgt*sqrt(ld_snpinfo_wgt$variance[ld_snpinfo_wgt.idx])
      }

      nwgt <- nrow(wgt.matrix)

      if(nwgt>0){
        p0 <- min(snps[snps[, "id"] %in% snpnames, "pos"])
        p1 <- max(snps[snps[, "id"] %in% snpnames, "pos"])
        weight_id <- paste0(gname, "|", weight_name)
        weights[[weight_id]] <- list(chrom=chrom, p0=p0, p1=p1, wgt=wgt, gene_name=gname, weight_name=weight_name, n=nwgt)
      }
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }

  return(weights)
}
