
#' @title Preprocess PredictDB/FUSION weights and harmonize with LD reference
#'
#' @param weight_path path to the '.db' file for PredictDB weights;
#' or the directory containing '.wgt.RDat' files for FUSION weights.
#'
#' @param region_info a data frame of region definition and associated file names.
#'
#' @param gwas_snp_ids a vector of SNP IDs in GWAS summary statistics (z_snp$id).
#'
#' @param snp_info a list of SNP info data frames for LD reference,
#'  with columns "chrom", "id", "pos", "alt", "ref", and "region_id".
#'
#' @param LD_info a list of paths to LD matrices for each of the regions. Required when \code{load_predictdb_LD = FALSE}.
#'
#' @param type a string, specifying QTL type of each weight file, e.g. expression, splicing, protein.
#'
#' @param context a string, specifying tissue/cell type/condition of each weight file, e.g. Liver, Lung, Brain.
#'
#' @param weight_format a string, specifying format of each weight file, e.g. PredictDB, FUSION.
#'
#' @param filter_protein_coding_genes TRUE/FALSE. If TRUE, keep protein coding genes only. This option is only for PredictDB weights
#'
#' @param load_predictdb_LD TRUE/FALSE. If TRUE, load pre-computed LD among weight SNPs. This option is only for PredictDB weights
#'
#' @param drop_strand_ambig TRUE/FALSE, if TRUE remove strand ambiguous variants (A/T, G/C).
#'
#' @param scale_by_ld_variance TRUE/FALSE, if TRUE scale PredictDB weights by the variance.
#' This is because PredictDB weights assume that variant genotypes are not
#' standardized, but our implementation assumes standardized variant genotypes.
#'
#' @param fusion_method a string, specifying the method to choose in FUSION models
#'
#' @param fusion_genome_version a string, specifying the genome version of FUSION models
#'
#' @param fusion_top_n_snps a number, specifying the top n weight SNPs included in FUSION models. If NULL, using all weight SNPs
#'
#' @param LD_format file format for LD matrix. If "custom", use a user defined
#' \code{LD_loader()} function to load LD matrix.
#'
#' @param LD_loader a user defined function to load LD matrix when \code{LD_format = "custom"}.
#'
#' @param ncore The number of cores used to parallelize computation.
#'
#' @param logfile the log file, if NULL will print log info on screen.
#'
#' @return a list of processed weights
#'
#' @importFrom utils head
#' @importFrom logging addHandler loginfo writeToFile
#' @importFrom data.table rbindlist
#' @importFrom tools file_path_sans_ext
#' @importFrom stats complete.cases
#'
#' @export
#'
preprocess_weights <- function(weight_path,
                               region_info,
                               gwas_snp_ids,
                               type,
                               context,
                               snp_info,
                               LD_info = NULL,
                               weight_format = c("PredictDB", "FUSION"),
                               drop_strand_ambig = TRUE,
                               scale_by_ld_variance = FALSE,
                               filter_protein_coding_genes = FALSE,
                               load_predictdb_LD = FALSE,
                               fusion_method = c("lasso","enet","top1","blup"),
                               fusion_genome_version = c("b38","b37"),
                               fusion_top_n_snps = NULL,
                               LD_format = c("rds", "rdata", "csv", "txt", "custom"),
                               LD_loader = NULL,
                               ncore = 1,
                               logfile = NULL){
  if (!is.null(logfile)) {
    addHandler(writeToFile, file = logfile, level = "DEBUG")
  }
  # check input arguments
  weight_format <- match.arg(weight_format)
  fusion_method <- match.arg(fusion_method)
  fusion_genome_version <- match.arg(fusion_genome_version)
  LD_format <- match.arg(LD_format)

  if (length(weight_path) != 1) {
    stop("Please provide only one weight path in `weight_path`.")
  }

  # Check LD reference SNP info
  if (inherits(snp_info,"list")) {
    snp_info_df <- as.data.frame(rbindlist(snp_info, idcol = "region_id"))
  }
  target_header <- c("chrom", "id", "pos", "alt", "ref")
  if (!all(target_header %in% colnames(snp_info_df))){
    stop("SNP info needs to contain the following columns: ",
         paste(target_header, collapse = " "))
  }

  # set default type and context
  if (missing(type)) {
    type <- "gene"
  }

  if (missing(context)) {
    context <- file_path_sans_ext(basename(weight_path))
  }

  loginfo("Load weight: %s", weight_path)
  loginfo("type: %s", type)
  loginfo("context: %s", context)

  weights <- list()
  loaded_weights_res <- load_weights(weight_path,
                                     weight_format,
                                     filter_protein_coding_genes = filter_protein_coding_genes,
                                     load_predictdb_LD = load_predictdb_LD,
                                     fusion_method = fusion_method,
                                     fusion_genome_version = fusion_genome_version,
                                     ncore=ncore)
  weight_table <- loaded_weights_res$weight_table
  weight_name <- loaded_weights_res$weight_name
  R_wgt_all <- loaded_weights_res$R_wgt
  if (!is.null(R_wgt_all)) {
    #remove genes without predictdb LD
    weight_table <- weight_table[weight_table$gene %in% unique(R_wgt_all$GENE), ]
  }
  gnames <- unique(weight_table$gene)
  loginfo("Number of genes with weights provided: %d in %s", length(gnames), weight_name)
  # remove variants in weight table, but not in LD reference and GWAS
  loginfo("Number of variants in weights: %d", length(unique(weight_table$rsid)))
  # take the intersect of SNPs in weights, LD reference and SNPs in z_snp
  snpnames <- Reduce(intersect, list(weight_table$rsid, snp_info_df$id, gwas_snp_ids))
  weight_table <- weight_table[weight_table$rsid %in% snpnames, ]
  gnames <- unique(weight_table$gene)
  loginfo("%d variants and %d genes left after filtering by GWAS and reference SNPs",
          length(snpnames), length(gnames))
  # subset to variants in weight table
  snp_info_wgt <- snp_info_df[snp_info_df$id %in% weight_table$rsid,]

  loginfo("Harmonizing weights with LD reference ...")
  rsid_varID <- weight_table[,c("rsid", "varID")]
  for (i in 1:length(gnames)) {
    gname <- gnames[i]
    wgt <- weight_table[weight_table$gene==gname,]
    wgt.matrix <- as.matrix(wgt[, "weight", drop = F])
    rownames(wgt.matrix) <- wgt$rsid
    chrpos <- do.call(rbind, strsplit(wgt$varID, "_"))
    chrom <- unique(as.integer(gsub("chr", "", chrpos[, 1])))
    if (length(chrom) > 1) {
      stop("More than one chrom in weight for %s!", gname)
    }
    snp_pos <- as.integer(chrpos[, 2])
    wgt_ld_idx <- match(wgt$rsid, snp_info_wgt$id)
    snp_info_pos <- as.integer(snp_info_wgt$pos[wgt_ld_idx])
    if (any(snp_pos != snp_info_pos)) {
      # loginfo("Positions in weights are different from those in snp_info for %s.
      # Use the positions in snp_info.", gname)
      snp_pos <- snp_info_pos
    }
    snps <- data.frame(chrom = chrom,
                       id = wgt$rsid,
                       cm = "0",
                       pos = snp_pos,
                       alt = wgt$eff_allele,
                       ref = wgt$ref_allele,
                       stringsAsFactors = F)

    w <- harmonize_weights(wgt.matrix,
                           snps,
                           snp_info_wgt,
                           drop_strand_ambig = drop_strand_ambig)

    wgt.matrix <- w[["wgt"]]
    snps <- w[["snps"]]
    wgt.matrix <- wgt.matrix[wgt.matrix[, "weight"] != 0, , drop = F]
    wgt.matrix <- wgt.matrix[complete.cases(wgt.matrix), , drop = F]

    snpnames <- intersect(rownames(wgt.matrix), snp_info_wgt$id)
    wgt.idx <- match(snpnames, rownames(wgt.matrix))
    wgt <- wgt.matrix[wgt.idx, "weight", drop = F]

    if (weight_format == "FUSION") {
      wgt <- wgt[order(-abs(wgt[,"weight"])), , drop = F]
      if (!is.null(fusion_top_n_snps)) {
        wgt <- head(wgt,fusion_top_n_snps)
      }
      snpnames <- intersect(rownames(wgt), snp_info_wgt$id)
    }

    snps.idx <- match(snpnames, snps$id)
    snps <- snps[snps.idx,]

    if (scale_by_ld_variance){
      snp_info_wgt.idx <- match(snpnames, snp_info_wgt$id)
      wgt <- wgt*sqrt(snp_info_wgt$variance[snp_info_wgt.idx])
    }

    n_wgt <- nrow(wgt)
    if (n_wgt > 0) {
      p0 <- min(snps[snps[, "id"] %in% snpnames, "pos"])
      p1 <- max(snps[snps[, "id"] %in% snpnames, "pos"])
      # weight_id <- paste0(gname, "|", type, "|", context)
      weight_id <- paste0(gname, "|", weight_name)

      # Add LD matrix of weights
      if(!is.null(R_wgt_all)){
        R_wgt <- get_weight_LD(R_wgt_all, gname, rsid_varID)
        R_wgt <- R_wgt[snps$id, snps$id, drop=F]
      }
      else{
        R_wgt <- NULL
      }
      weights[[weight_id]] <- list(chrom=chrom, p0=p0, p1=p1,
                                   wgt=wgt, R_wgt=R_wgt,
                                   gene_name=gname, weight_name=weight_name,
                                   type = type, context = context,
                                   n_wgt=n_wgt)
    }
  }

  if (!load_predictdb_LD) {
    loginfo("Computing LD among variants in weights ...")
    weights <- compute_weight_LD_from_ref(weights,
                                          weight_name,
                                          region_info = region_info,
                                          LD_info = LD_info,
                                          snp_info = snp_info,
                                          LD_format = LD_format,
                                          LD_loader = LD_loader,
                                          ncore = ncore)
  }

  return(weights)
}


