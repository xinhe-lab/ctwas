
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
#' @param include_weight_LD If TRUE, include LD of variants in weights (R_wgt) in the weights object.
#' R_wgt is used for computing gene Z-scores.
#' If FALSE, will skip computing R_wgt. This could save running time
#' when using precomputed gene Z-scores.
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
#' @param varID_converter_fun a user defined function to convert
#' weight variant IDs to the reference variant format.
#'
#' @param ncore The number of cores used to parallelize computation.
#'
#' @param logfile The log filename. If NULL, will print log info on screen.
#'
#' @param verbose If TRUE, print detail messages.
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
                               include_weight_LD = TRUE,
                               fusion_method = c("lasso","enet","top1","blup","bslmm","best.cv"),
                               fusion_genome_version = NA,
                               top_n_snps = NULL,
                               LD_format = c("rds", "rdata", "csv", "txt", "custom"),
                               LD_loader_fun = NULL,
                               snpinfo_loader_fun = NULL,
                               varID_converter_fun = NULL,
                               ncore = 1,
                               logfile = NULL,
                               verbose = FALSE){

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

  if (!inherits(gwas_snp_ids,"character")){
    stop("'gwas_snp_ids' should be a character object.")
  }

  if (!is.null(LD_map)){
    if (!inherits(LD_map,"data.frame"))
      stop("'LD_map' should be a data frame!")
  } else {
    if (!load_predictdb_LD && include_weight_LD)
      stop("'LD_map' is required when load_predictdb_LD = FALSE & include_weight_LD = TRUE!")
  }

  if (weight_format == "FUSION") {
    load_predictdb_LD <- FALSE
    filter_protein_coding_genes <- FALSE
    scale_predictdb_weights <- FALSE
    if (is.null(LD_map) && include_weight_LD)
      stop("'LD_map' is required when using weight_format = 'FUSION' & include_weight_LD = TRUE!")
  }

  if (!inherits(snp_map,"list")){
    stop("'snp_map' should be a list object.")
  }

  if (!all(region_info$region_id %in% names(snp_map))){
    stop("Not all region IDs in 'region_info' are found in 'snp_map'!")
  }

  if (length(region_info$region_id) < length(snp_map)){
    loginfo("Select %d regions in region_info", length(region_info$region_id))
    snp_map <- snp_map[region_info$region_id]
  }

  loginfo("Load weight: %s", weight_file)
  loginfo("weight_name: %s", weight_name)
  loginfo("type: %s", type)
  loginfo("context: %s", context)

  if (!include_weight_LD) {
    loginfo("Do not include weight_LD (R_wgt) in weights")
    load_predictdb_LD <- FALSE
  }

  res <- load_weights(weight_file,
                      weight_format,
                      filter_protein_coding_genes = filter_protein_coding_genes,
                      load_predictdb_LD = load_predictdb_LD,
                      fusion_method = fusion_method,
                      fusion_genome_version = fusion_genome_version,
                      ncore = ncore)
  weight_table <- res$weight_table
  cov_table <- res$cov_table
  rm(res)

  # convert variant IDs in weights to reference variant format
  if (!is.null(varID_converter_fun)) {
    loginfo("Convert variant IDs")
    weight_table$rsid <- varID_converter_fun(weight_table$rsid)
    weight_table$varID <- varID_converter_fun(weight_table$varID)
    if (!is.null(cov_table)) {
      cov_table$RSID1 <- varID_converter_fun(cov_table$RSID1)
      cov_table$RSID2 <- varID_converter_fun(cov_table$RSID2)
    }
  }

  if (is.null(cov_table)) {
    load_predictdb_LD <- FALSE
  }

  # remove genes without predictdb LD
  if (load_predictdb_LD && !is.null(cov_table)) {
    genes_without_predictdb_LD <- setdiff(weight_table$gene, cov_table$GENE)
    if (length(genes_without_predictdb_LD) > 0){
      loginfo("Remove %d molecular traits without PredictDB LD", length(genes_without_predictdb_LD))
      weight_table <- weight_table[!weight_table$gene %in% genes_without_predictdb_LD, ]
    }
  }
  weight_molecular_ids <- unique(weight_table$gene)
  weight_snp_ids <- unique(weight_table$rsid)
  loginfo("Number of molecular traits in weights: %d", length(weight_molecular_ids))
  loginfo("Number of variants in weights: %d", length(weight_snp_ids))

  # take the intersection of SNPs in weights, LD reference and GWAS SNPs
  snp_info <- as.data.frame(rbindlist(snp_map, idcol = "region_id"))
  weight_snp_ids <- Reduce(intersect, list(weight_snp_ids, snp_info$id, gwas_snp_ids))
  weight_table <- weight_table[weight_table$rsid %in% weight_snp_ids, ]
  weight_molecular_ids <- unique(weight_table$gene)
  loginfo("%d molecular traits and %d variants in weights after filtering by GWAS and the reference.",
          length(weight_molecular_ids), length(weight_snp_ids))
  # subset to variants in weight table
  snp_info <- snp_info[snp_info$id %in% weight_table$rsid,]
  # subset to genes in weight table
  cov_table <- cov_table[cov_table$GENE %in% weight_table$gene, ]

  loginfo("Harmonizing and processing weights ...")
  if (load_predictdb_LD && include_weight_LD) {
    loginfo("Compute LD for variants in weights (R_wgt) using PredictDB LD ...")
  }
  weights <- mclapply_check(weight_molecular_ids, function(molecular_id){
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
                   scale_predictdb_weights = scale_predictdb_weights,
                   verbose = verbose)
  }, mc.cores = ncore)
  names(weights) <- paste0(weight_molecular_ids, "|", weight_name)

  empty_wgt_idx <- which(sapply(weights, "[[", "n_wgt") == 0)
  if (length(empty_wgt_idx) > 0) {
    loginfo("Remove %d molecular traits with no weights after harmonization", length(empty_wgt_idx))
    weights[empty_wgt_idx] <- NULL
  }

  if (!load_predictdb_LD && include_weight_LD) {
    loginfo("Computing LD for variants in weights (R_wgt) using reference LD ...")
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
#' @importFrom logging loginfo logwarn
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
                           scale_predictdb_weights = TRUE,
                           verbose = FALSE) {

  if (weight_format == "FUSION") {
    scale_predictdb_weights <- FALSE
  }

  # Check LD reference SNP info
  target_header <- c("chrom", "id", "pos", "alt", "ref")
  if (!all(target_header %in% colnames(snp_info))){
    stop("snp_info needs to contain the following columns: ",
         paste(target_header, collapse = " "))
  }

  if (verbose)
    loginfo("Processing weight for %s ...", molecular_id)

  g.weight_table <- weight_table[weight_table$gene==molecular_id,]
  wgt.matrix <- as.matrix(g.weight_table[, "weight", drop = FALSE])
  rownames(wgt.matrix) <- g.weight_table$rsid

  chrom <- sapply(strsplit(g.weight_table$varID, "_"), "[[", 1)
  chrom <- unique(parse_number(chrom))
  if (length(chrom) > 1) {
    stop("More than one chrom in weight for %s!", molecular_id)
  }

  snp_pos <- as.integer(snp_info$pos[match(g.weight_table$rsid, snp_info$id)])
  wgt.snpinfo <- data.frame(chrom = chrom,
                            id = g.weight_table$rsid,
                            cm = 0,
                            pos = snp_pos,
                            alt = g.weight_table$eff_allele,
                            ref = g.weight_table$ref_allele,
                            stringsAsFactors = FALSE)

  if (verbose)
    loginfo("Harmonize weight")
  harmonized_wgt_res <- harmonize_weights(wgt.matrix,
                                          wgt.snpinfo,
                                          snp_info,
                                          drop_strand_ambig = drop_strand_ambig)
  wgt.matrix <- harmonized_wgt_res$wgt.matrix
  wgt.snpinfo <- harmonized_wgt_res$wgt.snpinfo
  wgt.matrix <- wgt.matrix[wgt.matrix[, "weight"] != 0, , drop = FALSE]
  wgt.matrix <- wgt.matrix[complete.cases(wgt.matrix), , drop = FALSE]

  wgt.snp_ids <- intersect(rownames(wgt.matrix), snp_info$id)
  wgt <- wgt.matrix[match(wgt.snp_ids, rownames(wgt.matrix)), "weight", drop = FALSE]

  # use top n SNPs in weights
  if (!is.null(top_n_snps)) {
    wgt <- wgt[order(-abs(wgt[,"weight"])), ]
    wgt <- head(wgt,top_n_snps)
    wgt.snp_ids <- rownames(wgt)
  }

  # scale weights by variance from LD reference
  if (weight_format == "PredictDB") {
    if (scale_predictdb_weights){
      if (verbose) {
        loginfo("Scale weights by variance from LD reference")
      }
      wgt_snp_var <- snp_info$variance[match(wgt.snp_ids, snp_info$id)]
      wgt <- wgt*sqrt(wgt_snp_var)
    }
  }

  wgt.snpinfo <- wgt.snpinfo[match(wgt.snp_ids, wgt.snpinfo$id),]
  g.weight_table <- g.weight_table[match(wgt.snp_ids, g.weight_table$rsid),]

  n_wgt <- nrow(wgt)
  if (n_wgt > 0) {
    p0 <- min(wgt.snpinfo$pos, na.rm = TRUE)
    p1 <- max(wgt.snpinfo$pos, na.rm = TRUE)

    # Add LD matrix of weights (R_wgt)
    if (!is.null(cov_table)) {
      # compute LD of variants in weights using PredictedDB cov_table
      if (verbose) {
        loginfo("Compute LD of variants in weights (R_wgt) from PredictedDB LD")
      }
      g.cov_table <- cov_table[cov_table$GENE == molecular_id,]
      R_wgt <- get_weight_LD_from_predictdb(g.cov_table,
                                        g.weight_table,
                                        convert_cov_to_cor = TRUE)
      R_wgt <- R_wgt[wgt.snp_ids, wgt.snp_ids, drop=FALSE]
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
