
#' @title Preprocess PredictDB/FUSION weights and harmonize with LD reference
#'
#' @param weight_file filename of the '.db' file for PredictDB weights;
#' or the directory containing '.wgt.RDat' files for FUSION weights.
#'
#' @param region_info a data frame of region definitions.
#'
#' @param gwas_snp_ids a vector of SNP IDs in GWAS summary statistics (z_snp$id).
#'
#' @param snp_map a list of SNP-to-region map for the reference.
#'
#' @param LD_map a data frame with filenames of LD matrices for the regions.
#' Required when \code{load_predictdb_LD = FALSE}.
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
#' @param scale_predictdb_weights TRUE/FALSE, if TRUE scale PredictDB weights by the variance.
#' This is because PredictDB weights assume that variant genotypes are not
#' standardized, but our implementation assumes standardized variant genotypes.
#'
#' @param fusion_method a string, specifying the method to choose in FUSION models.
#' "best.cv" option will use the best model (smallest p-value) under cross-validation.
#'
#' @param fusion_genome_version a string, specifying the genome version of FUSION models
#'
#' @param fusion_top_n_snps a number, specifying the top n weight SNPs included in FUSION models.
#' By default, use all weight SNPs.
#'
#' @param LD_format file format for LD matrix. If "custom", use a user defined
#' \code{LD_loader_fun()} function to load LD matrix.
#'
#' @param LD_loader_fun a user defined function to load LD matrix when \code{LD_format = "custom"}.
#'
#' @param ncore The number of cores used to parallelize computation.
#'
#' @param logfile the log file, if NULL will print log info on screen.
#'
#' @return a list of processed weights
#'
#' @importFrom logging addHandler loginfo writeToFile
#' @importFrom data.table rbindlist
#' @importFrom tools file_path_sans_ext
#'
#' @export
#'
preprocess_weights <- function(weight_file,
                               region_info,
                               gwas_snp_ids,
                               type,
                               context,
                               snp_map,
                               LD_map,
                               weight_format = c("PredictDB", "FUSION"),
                               drop_strand_ambig = TRUE,
                               scale_predictdb_weights = TRUE,
                               filter_protein_coding_genes = TRUE,
                               load_predictdb_LD = TRUE,
                               fusion_method = c("lasso","enet","top1","blup","bslmm","best.cv"),
                               fusion_genome_version = c("b38","b37"),
                               fusion_top_n_snps,
                               LD_format = c("rds", "rdata", "csv", "txt", "custom"),
                               LD_loader_fun = NULL,
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

  if (length(weight_file) != 1) {
    stop("Please provide only one weight path in `weight_file`.")
  }

  snp_info <- as.data.frame(rbindlist(snp_map, idcol = "region_id"))

  # set default type and context
  if (missing(type)) {
    type <- "gene"
  }

  if (missing(context)) {
    context <- file_path_sans_ext(basename(weight_file))
  }

  loginfo("Load weight: %s", weight_file)
  loginfo("type: %s", type)
  loginfo("context: %s", context)

  loaded_weights_res <- load_weights(weight_file,
                                     weight_format,
                                     filter_protein_coding_genes = filter_protein_coding_genes,
                                     load_predictdb_LD = load_predictdb_LD,
                                     fusion_method = fusion_method,
                                     fusion_genome_version = fusion_genome_version,
                                     ncore = ncore)
  weight_name <- loaded_weights_res$weight_name
  weight_table <- loaded_weights_res$weight_table
  cov_table <- loaded_weights_res$cov_table
  loginfo("weight_name: %s", weight_name)
  if (!is.null(cov_table)) {
    # remove genes without predictdb LD
    genes_with_predictdb_LD <- unique(cov_table$GENE)
    loginfo("Remove %d genes without predictdb LD", length(setdiff(weight_table$gene, genes_with_predictdb_LD)))
    weight_table <- weight_table[weight_table$gene %in% genes_with_predictdb_LD, ]
  }
  gene_names <- unique(weight_table$gene)
  loginfo("Number of genes with weights: %d", length(gene_names))
  snpnames <- unique(weight_table$rsid)
  loginfo("Number of variants in weights: %d", length(snpnames))
  # take the intersect of SNPs in weights, LD reference and GWAS SNPs
  snpnames <- Reduce(intersect, list(snpnames, snp_info$id, gwas_snp_ids))
  weight_table <- weight_table[weight_table$rsid %in% snpnames, ]
  gene_names <- unique(weight_table$gene)
  loginfo("%d genes and %d variants left after filtering by GWAS and LD reference",
          length(snpnames), length(gene_names))
  # subset to variants in weight table
  snp_info <- snp_info[snp_info$id %in% weight_table$rsid,]

  loginfo("Harmonizing and processing weights ...")

  weights <- mclapply_check(gene_names, function(gene_name){
    process_weight(gene_name,
                   type = type,
                   context = context,
                   weight_name = weight_name,
                   weight_table = weight_table,
                   cov_table = cov_table,
                   snp_info = snp_info,
                   weight_format = weight_format,
                   fusion_top_n_snps = fusion_top_n_snps,
                   drop_strand_ambig = drop_strand_ambig,
                   scale_predictdb_weights = scale_predictdb_weights)
  }, mc.cores = ncore)

  names(weights) <- paste0(gene_names, "|", weight_name)

  empty_wgt_idx <- which(sapply(weights, "[[", "n_wgt") == 0)
  if (any(empty_wgt_idx)) {
    loginfo("Remove %d empty weights after harmonization", length(empty_wgt_idx))
    weights[empty_wgt_idx] <- NULL
  }

  if (!load_predictdb_LD) {
    loginfo("Computing LD among variants in weights ...")
    weights <- compute_weight_LD_from_ref(weights,
                                          weight_name,
                                          region_info = region_info,
                                          LD_map = LD_map,
                                          snp_map = snp_map,
                                          LD_format = LD_format,
                                          LD_loader_fun = LD_loader_fun,
                                          ncore = ncore)
  }

  loginfo("Number of genes with weights after preprocessing: %d", length(weights))

  return(weights)
}

# Process weights for a single gene
#' @importFrom readr parse_number
#' @importFrom stats complete.cases
#' @importFrom utils head
process_weight <- function(gene_name,
                           type,
                           context,
                           weight_name,
                           weight_table,
                           cov_table,
                           snp_info,
                           weight_format = c("PredictDB", "FUSION"),
                           fusion_top_n_snps = NULL,
                           drop_strand_ambig = TRUE,
                           scale_predictdb_weights = TRUE) {

  # Check LD reference SNP info
  target_header <- c("chrom", "id", "pos", "alt", "ref")
  if (!all(target_header %in% colnames(snp_info))){
    stop("snp_info needs to contain the following columns: ",
         paste(target_header, collapse = " "))
  }

  g_weight_table <- weight_table[weight_table$gene==gene_name,]
  wgt.matrix <- as.matrix(g_weight_table[, "weight", drop = FALSE])
  rownames(wgt.matrix) <- g_weight_table$rsid
  chrom <- sapply(strsplit(g_weight_table$varID, "_"), "[[", 1)
  chrom <- unique(parse_number(chrom))
  if (length(chrom) > 1) {
    stop("More than one chrom in weight for %s!", gene_name)
  }

  snp_pos <- as.integer(snp_info$pos[match(g_weight_table$rsid, snp_info$id)])

  snps <- data.frame(chrom = chrom,
                     id = g_weight_table$rsid,
                     cm = 0,
                     pos = snp_pos,
                     alt = g_weight_table$eff_allele,
                     ref = g_weight_table$ref_allele,
                     stringsAsFactors = FALSE)

  harmonized_wgt_res <- harmonize_weights(wgt.matrix,
                                          snps,
                                          snp_info,
                                          drop_strand_ambig = drop_strand_ambig)
  wgt.matrix <- harmonized_wgt_res$wgt.matrix
  snps <- harmonized_wgt_res$snps
  wgt.matrix <- wgt.matrix[wgt.matrix[, "weight"] != 0, , drop = FALSE]
  wgt.matrix <- wgt.matrix[complete.cases(wgt.matrix), , drop = FALSE]

  snpnames <- intersect(rownames(wgt.matrix), snp_info$id)
  wgt.idx <- match(snpnames, rownames(wgt.matrix))
  wgt <- wgt.matrix[wgt.idx, "weight", drop = FALSE]

  if (weight_format == "FUSION") {
    wgt <- wgt[order(-abs(wgt[,"weight"])), , drop = FALSE]
    if (!is.null(fusion_top_n_snps)) {
      wgt <- head(wgt,fusion_top_n_snps)
    }
    snpnames <- rownames(wgt)
  }

  if (weight_format == "PredictDB") {
    if (scale_predictdb_weights){
      wgt_snp_var <- snp_info$variance[match(snpnames, snp_info$id)]
      wgt <- wgt*sqrt(wgt_snp_var)
    }
  }

  snps <- snps[match(snpnames, snps$id),]

  n_wgt <- nrow(wgt)
  if (n_wgt > 0) {
    p0 <- min(snps[snps[, "id"] %in% snpnames, "pos"])
    p1 <- max(snps[snps[, "id"] %in% snpnames, "pos"])
    # Add LD matrix of weights
    if (!is.null(cov_table)) {
      # get pre-computed LD matrix from predictedDB weights
      g_cov_table <- cov_table[cov_table$GENE == gene_name,]
      R_wgt <- get_LD_matrix_from_predictdb(g_cov_table,
                                            g_weight_table,
                                            convert_cov_to_cor = TRUE)
      R_wgt <- R_wgt[snps$id, snps$id, drop=FALSE]
    } else{
      R_wgt <- NULL
    }
    return(list(chrom = chrom,
                p0 = p0,
                p1 = p1,
                wgt = wgt,
                R_wgt = R_wgt,
                gene_name = gene_name,
                weight_name = weight_name,
                type = type,
                context = context,
                n_wgt = n_wgt))
  }

}
