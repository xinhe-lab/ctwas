
#' Preprocess PredictDB weights and harmonize with LD reference
#' (adapted from preharmonize_wgt_ld)
#'
#' @param weight_file a string, pointing to a directory with the fusion/twas format of weights, or a .db file in predictdb format.
#' A vector of multiple sets of weights in PredictDB format can also be specified; genes will have their filename appended
#' to their gene name to ensure IDs are unique.
#'
#' @param region_info a data frame of region definition and associated file names.
#'
#' @param drop_strand_ambig TRUE/FALSE, if TRUE remove strand ambiguous variants (A/T, G/C).
#'
#' @return a list of processed weight table and extra table
#'
#' @export
#'
process_weight <- function (weight_file,
                            region_info,
                            weight_format = c("PredictDB", "Fusion"),
                            drop_strand_ambig = TRUE,
                            scale_by_ld_variance = F){

  weight_format <- match.arg(weight_format)
  # load LD SNPs information
  region_info <- region_info[order(region_info$chrom, region_info$start),]
  ld_snpinfo <- read_LD_SNP_files(region_info$SNP_info)
  outlist <- list()
  for(weight in weight_file){
    loaded_weight <- load_weight(weight, weight_format = weight_format)
    weight_table <- loaded_weight$weight_table
    weight_name <- loaded_weight$weight_name
    gnames <- unique(weight_table$gene)
    loginfo("Number of genes with weights provided: %s in %s", length(gnames), weight_name)

    # remove variants in weight table, but not in LD reference
    loginfo("Number of variants in weights: %s", length(unique(weight_table$rsid)))
    loginfo("Remove %s variants in weights but not in LD reference", length(setdiff(weight_table$rsid, ld_snpinfo$id)))
    weight_table <- weight_table[weight_table$rsid %in% ld_snpinfo$id, ]

    # remove genes with no variants in LD reference
    loginfo("Remove %s genes with no variants in LD reference", length(setdiff(gnames, weight_table$gene)))
    gnames <- unique(weight_table$gene)
    loginfo("Number of genes left after removing variants not in LD reference: %s", length(gnames))

    # subset to variants in weight table
    ld_snpinfo_wgt <- ld_snpinfo[ld_snpinfo$id %in% weight_table$rsid,]
    loginfo("Processing weights for %s genes to match LD reference", length(gnames))

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
                            #ld_snpinfo,
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

      if (scale_by_ld_variance){
        ld_snpinfo_wgt.idx <- match(snpnames, ld_snpinfo_wgt$id)
        wgt <- wgt*sqrt(ld_snpinfo_wgt$variance[ld_snpinfo_wgt.idx])
      }

      nwgt <- nrow(wgt.matrix)
      if(nwgt>0){
        p0 <- min(snps[snps[, "id"] %in% snpnames, "pos"])
        p1 <- max(snps[snps[, "id"] %in% snpnames, "pos"])
        gname_weight <- paste0(gname, "|", weight_name)
        outlist[[gname_weight]] <- list(chrom = chrom, p0 = p0, p1 = p1, wgt = wgt, gname=gname, weight_name=weight_name, n = nwgt)
      }
    }
  }

  weight_list <- lapply(outlist, "[[", "wgt")
  names(weight_list) <- names(outlist)

  weight_info <- lapply(names(outlist), function(x){
    as.data.frame(outlist[[x]][c("chrom", "p0","p1", "gname", "weight_name", "n")])})
  weight_info <- do.call(rbind, weight_info)
  weight_info$id <- names(outlist)
  rownames(weight_info) <- names(outlist)

  return(list(weight_list = weight_list, weight_info = weight_info))
}
