
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
#' @param LD_map a data frame with filenames of LD matrices and SNP information for the regions.
#' Required when \code{load_predictdb_LD = FALSE}.
#'
#' @param type a string, specifying QTL type of the weight file, e.g. expression, splicing, protein.
#'
#' @param context a string, specifying context (tissue/cell type) of the weight file, e.g. Liver, Lung, Brain.
#'
#' @param weight_name a string, specifying name of the weight file.
#' By default, it is \code{weight_name = paste0(context, "_", type)}
#'
#' @param weight_format a string, specifying format of the weight file, e.g. PredictDB, FUSION.
#'
#' @param drop_strand_ambig If TRUE remove strand ambiguous variants (A/T, G/C).
#'
#' @param filter_protein_coding_genes If TRUE, keep protein coding genes only.
#' This option is only for PredictDB weights.
#'
#' @param scale_predictdb_weights If TRUE, scale PredictDB weights by the variance.
#' This is because PredictDB weights assume that variant genotypes are not
#' standardized, but our implementation assumes standardized variant genotypes.
#' This option is only for PredictDB weights.
#'
#' @param load_predictdb_LD If TRUE, load pre-computed LD among weight SNPs.
#' This option is only for PredictDB weights.
#'
#' @param fusion_method a string, specifying the method to choose in FUSION models.
#' "best.cv" option will use the best model (smallest p-value) under cross-validation.
#'
#' @param fusion_genome_version a string, specifying the genome version of FUSION models
#'
#' @param top_n_snps a number, specifying the top n SNPs included in weight models.
#' By default, use all SNPs in weights.
#'
#' @param LD_format file format for LD matrix. If "custom", use a user defined
#' \code{LD_loader_fun()} function to load LD matrix.
#'
#' @param LD_loader_fun a user defined function to load LD matrix when \code{LD_format = "custom"}.
#'
#' @param snpinfo_loader_fun a user defined function to load SNP information file,
#' if SNP information files are not in standard cTWAS reference format.
#'
#' @param ncore The number of cores used to parallelize computation.
#'
#' @param logfile The log filename. If NULL, will print log info on screen.
#'
#' @return a list of processed weights
#'
#' @importFrom logging addHandler loginfo logwarn writeToFile
#' @importFrom data.table rbindlist
#' @importFrom tools file_path_sans_ext
#'
#' @export
#'
preprocess_weights <- function(weight_file,
                               region_info,
                               gwas_snp_ids,
                               snp_map,
                               LD_map = NULL,
                               type,
                               context,
                               weight_name = paste0(context, "_", type),
                               weight_format = c("PredictDB", "FUSION"),
                               drop_strand_ambig = TRUE,
                               filter_protein_coding_genes = TRUE,
                               scale_predictdb_weights = TRUE,
                               load_predictdb_LD = TRUE,
                               fusion_method = c("lasso","enet","top1","blup","bslmm","best.cv"),
                               fusion_genome_version = NA,
                               top_n_snps = NULL,
                               LD_format = c("rds", "rdata", "csv", "txt", "custom"),
                               LD_loader_fun = NULL,
                               snpinfo_loader_fun = NULL,
                               ncore = 1,
                               logfile = NULL){
  if (!is.null(logfile)) {
    addHandler(writeToFile, file = logfile, level = "DEBUG")
  }

  # check input arguments
  weight_format <- match.arg(weight_format)
  fusion_method <- match.arg(fusion_method)
  LD_format <- match.arg(LD_format)

  if (length(weight_file) != 1) {
    stop("Please provide only one weight file in `weight_file`.")
  }

  if (weight_format != "PredictDB") {
    load_predictdb_LD <- FALSE
    filter_protein_coding_genes <- FALSE
    scale_predictdb_weights <- FALSE
  }

  if (!load_predictdb_LD) {
    if (is.null(LD_map))
      stop("'LD_map' is required when load_predictdb_LD = FALSE!")

    if (!inherits(LD_map,"data.frame"))
      stop("'LD_map' should be a data frame!")
  }

  if (!inherits(snp_map,"list")){
    stop("'snp_map' should be a list object.")
  }

  snp_info <- as.data.frame(rbindlist(snp_map, idcol = "region_id"))

  loginfo("Load weight: %s", weight_file)
  loginfo("weight_name: %s", weight_name)
  loginfo("type: %s", type)
  loginfo("context: %s", context)

  loaded_weights_res <- load_weights(weight_file,
                                     weight_format,
                                     filter_protein_coding_genes = filter_protein_coding_genes,
                                     load_predictdb_LD = load_predictdb_LD,
                                     fusion_method = fusion_method,
                                     fusion_genome_version = fusion_genome_version,
                                     ncore = ncore)
  weight_table <- loaded_weights_res$weight_table
  cov_table <- loaded_weights_res$cov_table

  if (is.null(cov_table)) {
    load_predictdb_LD <- FALSE
  }

  # remove genes without predictdb LD
  if (load_predictdb_LD && !is.null(cov_table)) {
    genes_without_predictdb_LD <- setdiff(weight_table$gene, cov_table$GENE)
    if (length(genes_without_predictdb_LD) > 0){
      loginfo("Remove %d molecular traits without predictdb LD", length(genes_without_predictdb_LD))
      weight_table <- weight_table[!weight_table$gene %in% genes_without_predictdb_LD, ]
    }
  }
  molecular_ids <- unique(weight_table$gene)
  loginfo("Number of molecular traits in weights: %d", length(molecular_ids))
  snp_ids <- unique(weight_table$rsid)
  loginfo("Number of variants in weights: %d", length(snp_ids))

  # take the intersection of SNPs in weights, LD reference and GWAS SNPs
  snp_ids <- Reduce(intersect, list(snp_ids, snp_info$id, gwas_snp_ids))
  weight_table <- weight_table[weight_table$rsid %in% snp_ids, ]
  molecular_ids <- unique(weight_table$gene)
  loginfo("%d molecular traits and %d variants left after filtering by GWAS and the reference.",
          length(molecular_ids), length(snp_ids))
  # subset to variants in weight table
  snp_info <- snp_info[snp_info$id %in% weight_table$rsid,]

  loginfo("Harmonizing and processing weights ...")

  weights <- mclapply_check(molecular_ids, function(molecular_id){
    process_weight(molecular_id,
                   type = type,
                   context = context,
                   weight_name = weight_name,
                   weight_table = weight_table,
                   cov_table = cov_table,
                   snp_info = snp_info,
                   weight_format = weight_format,
                   top_n_snps = top_n_snps,
                   drop_strand_ambig = drop_strand_ambig,
                   scale_predictdb_weights = scale_predictdb_weights)
  }, mc.cores = ncore)
  names(weights) <- paste0(molecular_ids, "|", weight_name)

  empty_wgt_idx <- which(sapply(weights, "[[", "n_wgt") == 0)
  if (any(empty_wgt_idx)) {
    loginfo("Remove %d molecular traits with no weights after harmonization", length(empty_wgt_idx))
    weights[empty_wgt_idx] <- NULL
  }

  if (!load_predictdb_LD) {
    loginfo("Computing LD for variants in weights using reference LD matrices ...")
    weights <- compute_weight_LD_from_ref(weights,
                                          region_info = region_info,
                                          LD_map = LD_map,
                                          snp_map = snp_map,
                                          LD_format = LD_format,
                                          LD_loader_fun = LD_loader_fun,
                                          snpinfo_loader_fun = snpinfo_loader_fun,
                                          ncore = ncore)
  }

  loginfo("Number of molecular traits in weights after preprocessing: %d", length(weights))

  return(weights)
}

# Process weights for a single gene
#' @importFrom readr parse_number
#' @importFrom stats complete.cases
#' @importFrom utils head
process_weight <- function(molecular_id,
                           type,
                           context,
                           weight_name,
                           weight_table,
                           cov_table,
                           snp_info,
                           weight_format = c("PredictDB", "FUSION"),
                           top_n_snps = NULL,
                           drop_strand_ambig = TRUE,
                           scale_predictdb_weights = TRUE) {

  if (weight_format != "PredictDB") {
    scale_predictdb_weights <- FALSE
  }

  # Check LD reference SNP info
  target_header <- c("chrom", "id", "pos", "alt", "ref")
  if (!all(target_header %in% colnames(snp_info))){
    stop("snp_info needs to contain the following columns: ",
         paste(target_header, collapse = " "))
  }

  g_weight_table <- weight_table[weight_table$gene==molecular_id,]
  wgt.matrix <- as.matrix(g_weight_table[, "weight", drop = FALSE])
  rownames(wgt.matrix) <- g_weight_table$rsid
  chrom <- sapply(strsplit(g_weight_table$varID, "_"), "[[", 1)
  chrom <- unique(parse_number(chrom))
  if (length(chrom) > 1) {
    stop("More than one chrom in weight for %s!", molecular_id)
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

  snp_ids <- intersect(rownames(wgt.matrix), snp_info$id)
  wgt.idx <- match(snp_ids, rownames(wgt.matrix))
  wgt <- wgt.matrix[wgt.idx, "weight", drop = FALSE]

  # use top n snps in weights
  if (!is.null(top_n_snps)) {
    wgt <- wgt[order(-abs(wgt[,"weight"])), ]
    wgt <- head(wgt,top_n_snps)
    snp_ids <- rownames(wgt)
  }

  # scale weights by variance from LD reference
  if (weight_format == "PredictDB") {
    if (scale_predictdb_weights){
      wgt_snp_var <- snp_info$variance[match(snp_ids, snp_info$id)]
      wgt <- wgt*sqrt(wgt_snp_var)
    }
  }

  snps <- snps[match(snp_ids, snps$id),]

  n_wgt <- nrow(wgt)
  if (n_wgt > 0) {
    p0 <- min(snps[snps[, "id"] %in% snp_ids, "pos"])
    p1 <- max(snps[snps[, "id"] %in% snp_ids, "pos"])
    # Add LD matrix of weights
    if (!is.null(cov_table)) {
      # get pre-computed LD matrix from predictedDB weights
      g_cov_table <- cov_table[cov_table$GENE == molecular_id,]
      R_wgt <- get_LD_matrix_from_predictdb(g_cov_table,
                                            g_weight_table,
                                            convert_cov_to_cor = TRUE)
      R_wgt <- R_wgt[snps$id, snps$id, drop=FALSE]
    } else{
      R_wgt <- NULL
    }
  } else {
    p0 <- p1 <- NA
    wgt <- R_wgt <- NULL
  }

  return(list("chrom" = chrom,
              "p0" = p0,
              "p1" = p1,
              "wgt" = wgt,
              "R_wgt" = R_wgt,
              "molecular_id" = molecular_id,
              "weight_name" = weight_name,
              "type" = type,
              "context" = context,
              "n_wgt" = n_wgt))
}
