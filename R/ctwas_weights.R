
#' @title Loads weights in PredictDB or FUSION format
#'
#' @param weight_file filename of the '.db' file for PredictDB weights;
#' or the directory containing '.wgt.RDat' files for FUSION weights.
#'
#' @param weight_format a string or a vector, specifying format of each weight file, e.g. PredictDB, FUSION.
#'
#' @param filter_protein_coding_genes If TRUE, keep protein coding genes only. This option is only for PredictDB weights
#'
#' @param load_predictdb_LD If TRUE, load pre-computed LD among weight SNPs. This option is only for PredictDB weights
#'
#' @param fusion_method a string, specifying the method to choose in FUSION models.
#' "best.cv" option will use the best model (smallest p-value) under cross-validation.
#'
#' @param fusion_genome_version a string, specifying the genome version of FUSION models
#'
#' @param ncore integer, number of cores for parallel computing.
#'
#' @export
#'
load_weights <- function(weight_file,
                         weight_format = c("PredictDB", "FUSION"),
                         filter_protein_coding_genes = TRUE,
                         load_predictdb_LD = TRUE,
                         fusion_method = c("lasso","enet","top1","blup","bslmm","best.cv"),
                         fusion_genome_version = NA,
                         ncore = 1){

  weight_format <- match.arg(weight_format)
  fusion_method <- match.arg(fusion_method)

  if (weight_format == "PredictDB") {
    res <- load_predictdb_weights(weight_file,
                                  filter_protein_coding_genes = filter_protein_coding_genes,
                                  load_predictdb_LD = load_predictdb_LD)
  } else if (weight_format == "FUSION") {
    res <- load_fusion_weights(weight_file,
                               fusion_method = fusion_method,
                               fusion_genome_version = fusion_genome_version,
                               ncore = ncore)
  }

  return(res)
}

#' @title Loads weights in PredictDB format
#'
#' @param weight_file a string, pointing path to weights in PredictDB format.
#'
#' @param filter_protein_coding_genes TRUE/FALSE. If TRUE, keep protein coding
#' genes only.
#'
#' @param load_predictdb_LD TRUE/FALSE. If TRUE, load pre-computed LD
#' among weight variants.
#'
#' @importFrom tools file_path_sans_ext
#' @importFrom RSQLite dbDriver dbConnect dbGetQuery dbDisconnect
#' @importFrom utils read.table
#' @importFrom data.table fread
#' @importFrom logging loginfo logwarn
#'
#' @export
#'
load_predictdb_weights <- function(weight_file,
                                   filter_protein_coding_genes = TRUE,
                                   load_predictdb_LD = TRUE){

  # read the PredictDB weights
  stopifnot(file.exists(weight_file))

  loginfo("Load PredictDB weights")
  sqlite <- dbDriver("SQLite")
  db <- dbConnect(sqlite, weight_file)
  query <- function(...) dbGetQuery(db, ...)
  weight_table <- query("select * from weights")
  extra_table <- query("select * from extra")

  loginfo("Number of molecular traits in weights: %s", length(unique(weight_table$gene)))

  # subset to protein coding genes only
  if (filter_protein_coding_genes) {
    # filter only when protein coding genes exist
    if ("protein_coding" %in% extra_table$gene_type){
      loginfo("Limit to protein coding genes")
      extra_table <- extra_table[extra_table$gene_type=="protein_coding",,drop=FALSE]
      weight_table <- weight_table[weight_table$gene %in% extra_table$gene,]
      loginfo("Number of molecular traits in weights after filtering protein coding genes: %s", length(unique(weight_table$gene)))
    } else {
      loginfo("No 'protein_coding' in 'extra_table$gene_type'. Skipped filtering protein coding genes.")
    }
  }

  if (anyNA(weight_table)) {
    logwarn("weight_table contains NAs!")
  }

  # load pre-computed covariances from PredictDB LD file
  if (load_predictdb_LD) {
    loginfo("Load PredictDB LD")
    predictdb_LD_file <- paste0(file_path_sans_ext(weight_file), ".txt.gz")
    if (!file.exists(predictdb_LD_file)){
      stop(paste("PredictDB LD file", predictdb_LD_file, "does not exist!"))
    }
    cov_table <- as.data.frame(fread(predictdb_LD_file, header=TRUE))
    if (anyNA(cov_table)) {
      logwarn("cov_table contains NAs!")
    }
  }
  else{
    cov_table <- NULL
  }
  dbDisconnect(db)

  return(list("weight_table" = weight_table,
              "extra_table" = extra_table,
              "cov_table" = cov_table))
}

#' @title Loads weights in FUSION format
#'
#' @param weight_dir the directory containing FUSION weights ('.wgt.RDat' files).
#'
#' @param fusion_method a string, specifying the method to choose in FUSION models.
#' "best.cv" option will use the best model (smallest p-value) under cross-validation.
#'
#' @param fusion_genome_version a string, specifying the genome version of FUSION models
#'
#' @param make_extra_table If TRUE, make an extra table in predictDB format
#'
#' @param ncore integer, number of cores for parallel computing.
#'
#' @importFrom utils read.table
#' @importFrom stats ave
#' @importFrom tools file_path_sans_ext
#' @importFrom parallel mclapply
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise n ungroup
#' @importFrom rlang .data
#' @importFrom logging loginfo logwarn
#'
#' @export
#'
load_fusion_weights <- function(weight_dir,
                                fusion_method = c("lasso","enet","top1","blup","bslmm","best.cv"),
                                fusion_genome_version = NA,
                                make_extra_table = FALSE,
                                ncore = 1) {

  fusion_method <- match.arg(fusion_method)
  stopifnot(dir.exists(weight_dir))

  loginfo("Load FUSION weights")

  # list FUSION weight Rdata files
  wgt_dir <- dirname(weight_dir)
  wgt_pos_file <- file.path(wgt_dir, paste0(basename(weight_dir), ".pos"))
  if (file.exists(wgt_pos_file)) {
    wgt_pos <- read.table(wgt_pos_file, header = TRUE, stringsAsFactors = FALSE)
    wgt_pos$ID <-
      ifelse(duplicated(wgt_pos$ID) | duplicated(wgt_pos$ID,fromLast = TRUE),
             paste(wgt_pos$ID, ave(wgt_pos$ID,wgt_pos$ID,FUN = seq_along), sep="_ID"),
             wgt_pos$ID)
    wgt_pos <- wgt_pos[wgt_pos$ID!="NA_IDNA",] # filter NA genes
    wgt_rdata_files <- file.path(wgt_dir, wgt_pos[, "WGT"])
    wgt_IDs <- wgt_pos[, "ID"]
  } else {
    wgt_rdata_files <- list.files(weight_dir, full.names = TRUE)
    wgt_IDs <- gsub(".wgt.RDat", "", basename(wgt_rdata_files))
  }

  # Get list of all files in weight_dir
  if (length(wgt_rdata_files) == 0) {
    stop("No FUSION '.wgt.RDat' files found.")
  }

  loginfo("Load %d FUSION 'wgt' files", length(wgt_rdata_files))
  weight_table <- mclapply_check(1:length(wgt_rdata_files), function(i){
    loaded_wgt_res <- load_fusion_wgt_data(
      wgt_rdata_files[i],
      wgt_IDs[i],
      fusion_method = fusion_method,
      fusion_genome_version = fusion_genome_version)
    loaded_wgt_res$weight_table
  }, mc.cores = ncore)
  weight_table <- do.call(rbind, weight_table)

  if (anyNA(weight_table)) {
    logwarn("weight_table contains NAs!")
  }

  if (make_extra_table) {
    extra_table <- weight_table %>% group_by(.data$gene) %>%
      summarise(n.snps.in.model = n()) %>% ungroup()
    extra_table$genename <- NA
    extra_table$gene_type <- NA
    extra_table$pred.perf.R2 <- NA
    extra_table$pred.perf.pval <- NA
    extra_table$pred.perf.qval <- NA
    extra_table <- extra_table[, c("gene", "genename", "gene_type", "n.snps.in.model",
                                   "pred.perf.R2", "pred.perf.pval", "pred.perf.qval")]
  } else {
    extra_table <- NULL
  }

  cov_table <- NULL

  return(list("weight_table" = weight_table,
              "extra_table" = extra_table,
              "cov_table" = cov_table))
}


# load FUSION .wgt.RDat file
#' @importFrom dplyr left_join
load_fusion_wgt_data <- function(wgt_rdata_file,
                                 wgt_ID,
                                 fusion_method = c("lasso","enet","top1","blup","bslmm","best.cv"),
                                 fusion_genome_version = NA){

  fusion_method <- match.arg(fusion_method)

  # define the variables as NULL to binding the variables locally to the function.
  snps <- wgt.matrix <- cv.performance <- NULL

  # load FUSION wgt.RDat
  load(wgt_rdata_file)

  if (missing(wgt_ID)){
    wgt_ID <- gsub(".wgt.RDat", "", basename(wgt_rdata_file))
  }

  snps <- as.data.frame(snps)
  colnames(snps) <- c("chrom", "rsid", "cm", "pos", "alt", "ref")
  snps$varID <- sprintf("chr%s_%s_%s_%s_%s",
                        snps$chrom, snps$pos, snps$ref, snps$alt,
                        fusion_genome_version)
  # use varID for those missing rsIDs
  snps[is.na(snps$rsid),"rsid"] <- snps[is.na(snps$rsid),"varID"]
  snps <- snps[, c("chrom", "rsid", "varID", "pos", "ref", "alt")]

  rownames(wgt.matrix) <- snps$rsid

  # select FUSION method

  # which rows have rsq
  row.rsq <- grep("rsq" , rownames(cv.performance))
  # which rows have p-values
  row.pval <- grep("pval" , rownames(cv.performance))

  cv.rsq <- cv.performance[row.rsq,]
  if (fusion_method == "best.cv"){
    best.idx <- which.min(apply(cv.performance[row.pval,,drop=FALSE],2,min,na.rm=TRUE))
    g.method <- names(cv.performance)[best.idx]
  } else{
    g.method <- fusion_method
  }

  if (!g.method %in% colnames(wgt.matrix))
    stop(paste(g.method, "not found in the columns of wgt.matrix"))

  g.cv.rsq <- cv.rsq[g.method]

  # for top1 method, only top magnitude SNP weight have non-zero weights
  if(g.method == "top1"){
    wgt.matrix[, "top1"][-which.max(abs(wgt.matrix[, "top1"]))] <- 0
  }
  wgt.matrix <- wgt.matrix[wgt.matrix[, g.method] != 0, , drop = FALSE]
  wgt.matrix <- wgt.matrix[complete.cases(wgt.matrix), , drop = FALSE]
  if (nrow(wgt.matrix) > 0){
    weight_table <- data.frame(gene = wgt_ID,
                               rsid = rownames(wgt.matrix),
                               weight = wgt.matrix[,g.method])
    weight_table <- weight_table %>% left_join(snps, by = "rsid")
    weight_table <- weight_table[,c("gene","rsid","varID","ref","alt","weight")]
    colnames(weight_table) <- c("gene","rsid","varID","ref_allele","eff_allele","weight")
  } else {
    weight_table <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(weight_table) <- c("gene","rsid","varID","ref_allele","eff_allele","weight")
  }

  return(list("weight_table" = weight_table,
              "fusion_method" = g.method,
              "wgt_ID" = wgt_ID,
              "cv.rsq" = g.cv.rsq))
}

#' Creates weight files in PredictDB format from QTL data
#'
#' @param weight_table a data frame of the genes, QTLs and weights, with columns:
#' "gene", "rsid", "varID", "ref_allele", "eff_allele", "weight".
#' If you want to use multiple eQTLs per gene, you can set \code{use_top_QTL=FALSE}.
#' But we assume the weights of the eQTLs are learned from multiple regression
#' (instead of marginal effect sizes).
#'
#' @param gene_table a data frame (optional) with information of the genes
#' in \code{weight_table} ("gene","genename","gene_type", etc.).
#' If NULL, create a simple gene_table based on the weight_table
#'
#' @param cov_table a data frame of covariances between variants, with columns:
#' "GENE","RSID1","RSID2", "VALUE".
#' If NULL, do not create covariance files (.txg.gz), unless \code{use_top_QTL=TRUE}.
#'
#' @param use_top_QTL If TRUE, only keep the top QTL with
#' the largest abs(weight) for each gene (molecular trait), and
#' create a simple cov_table with covariance set to 1.
#'
#' @param outputdir output directory
#'
#' @param outname name of the output weight file
#'
#' @importFrom stats complete.cases
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise n ungroup
#' @importFrom rlang .data
#'
#' @export
#'
create_predictdb_from_QTLs <- function(weight_table,
                                       gene_table = NULL,
                                       cov_table = NULL,
                                       use_top_QTL = TRUE,
                                       outputdir = getwd(),
                                       outname){

  if (!dir.exists(outputdir))
    dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)

  loginfo("Makes PredictDB weights from QTL data")

  # check and clean the QTL data
  weight_table <- as.data.frame(weight_table)
  required_cols <- c("gene", "rsid", "varID", "ref_allele", "eff_allele", "weight")
  if (!all(required_cols %in% colnames(weight_table))){
    stop("QTL_data needs to contain the following columns: ",
         paste(required_cols, collapse = " "))
  }

  # if use_top_QTL, select the top SNP with the max abs(weight) for each gene
  if (use_top_QTL) {
    loginfo("select the top SNP with the max abs(weight) for each gene")
    weight_table <- weight_table[with(weight_table, order(gene, -abs(weight))),]
    weight_table <- weight_table[!duplicated(weight_table$gene), ]
  }

  weight_table <- weight_table[weight_table$weight != 0, ,drop = FALSE]
  weight_table <- weight_table[complete.cases(weight_table), ,drop = FALSE]

  # if NULL, create a simply extra_table based on weight_table
  if (is.null(gene_table)) {
    gene_table <- weight_table %>%
      group_by(.data$gene) %>%
      summarise(n.snps.in.model = n()) %>%
      ungroup() %>% as.data.frame()
    gene_table$genename <- NA
    gene_table$gene_type <- NA
    gene_table$pred.perf.R2 <- NA
    gene_table$pred.perf.pval <- NA
    gene_table$pred.perf.qval <- NA
    gene_table <- gene_table[, c("gene", "genename", "gene_type", "n.snps.in.model",
                                 "pred.perf.R2", "pred.perf.pval", "pred.perf.qval")]
  }

  if (use_top_QTL) {
    if (any(gene_table$n.snps.in.model > 1)){
      stop("each gene should have only one SNP when using top QTL only")
    }
    # set covariance to 1 as each gene only has one top QTL
    cov_table <- weight_table[,c("gene","varID","varID")]
    colnames(cov_table) <- c("GENE","RSID1","RSID2")
    cov_table$VALUE <- 1
  }

  # write PredictDB '.db' file
  write_predictdb(weight_table, gene_table, cov_table, outputdir, outname)

}

#' @title Converts fusion weights to predictDB format
#'
#' @param weight_dir the directory containing FUSION weights ('.wgt.RDat' files).
#'
#' @param fusion_method a string, specifying the method to choose in FUSION models.
#' "best.cv" option will use the best model (smallest p-value) under cross-validation.
#'
#' @param fusion_genome_version a string, specifying the genome version of FUSION models
#'
#' @param make_extra_table If TRUE, make an extra table in predictDB format
#'
#' @param cov_table a data frame of covariances between variants, with columns:
#' "GENE","RSID1","RSID2","VALUE".
#' If NULL, do not create covariance files (.txg.gz).
#'
#' @param outputdir output directory
#'
#' @param outname name of the output weight file
#'
#' @export
#'
convert_fusion_to_predictdb <- function(
    weight_dir,
    fusion_method = c("lasso","enet","top1","blup","bslmm","best.cv"),
    fusion_genome_version = NA,
    make_extra_table = TRUE,
    cov_table = NULL,
    outputdir = getwd(),
    outname){

  fusion_method <- match.arg(fusion_method)

  loaded_weights_res <- load_fusion_weights(weight_dir,
                                            fusion_method = fusion_method,
                                            fusion_genome_version = fusion_genome_version,
                                            make_extra_table = make_extra_table)
  weight_table <- loaded_weights_res$weight_table
  extra_table <- loaded_weights_res$extra_table

  if (missing(outname)){
    outname <- file_path_sans_ext(basename(weight_dir))
  }

  # write PredictDB weights
  write_predictdb(weight_table, extra_table, cov_table, outputdir, outname)

  return(list("weight_table" = weight_table,
              "extra_table" = extra_table,
              "cov_table" = cov_table))
}


# Makes PredictDB '.db' file from weight_table and extra_table
#' @importFrom utils write.table
#' @importFrom RSQLite dbDriver dbConnect dbWriteTable dbDisconnect
write_predictdb <- function(weight_table,
                            extra_table = NULL,
                            cov_table = NULL,
                            outputdir,
                            outname) {

  # check required columns
  required_cols <- c("gene", "rsid", "varID", "ref_allele", "eff_allele", "weight")
  if (!all(required_cols %in% colnames(weight_table))){
    stop("weight_table needs to contain the following columns: ",
         paste(required_cols, collapse = " "))
  }

  # Create a database connection
  driver <- dbDriver('SQLite')
  db <- dbConnect(drv = driver, file.path(outputdir, paste0(outname,".db")))
  # create weights table
  dbWriteTable(db, 'weights', weight_table, overwrite = TRUE)

  # create an empty extra table if NULL
  if (is.null(extra_table)) {
    extra_table <- data.frame(matrix(ncol = 7, nrow = 0))
    colnames(extra_table) <- c("gene", "genename", "gene_type", "n.snps.in.model",
                               "pred.perf.R2", "pred.perf.pval", "pred.perf.qval")
  }
  # check required columns
  required_cols <- c("gene", "genename", "gene_type")
  if (!all(required_cols %in% colnames(extra_table))) {
    stop("extra_table needs to contain the following columns: ",
         paste(required_cols, collapse = " "))
  }
  dbWriteTable(db, 'extra', extra_table, overwrite = TRUE)
  dbDisconnect(db)

  if (!is.null(cov_table)) {
    # check required columns
    required_cols <- c("GENE", "RSID1", "RSID2", "VALUE")
    if (!all(required_cols %in% colnames(cov_table))) {
      stop("cov_table needs to contain the following columns: ",
           paste(required_cols, collapse = " "))
    }
    write.table(cov_table,
                file = gzfile(file.path(outputdir,paste0(outname,".txt.gz"))),
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  }

}


# Gets LD for a gene from precomputed PredictDB covariances
#' @importFrom utils combn
get_weight_LD_from_predictdb <- function (g.cov_table,
                                          g.weight_table,
                                          convert_cov_to_cor = TRUE){

  # convert covariances to correlations
  if (convert_cov_to_cor){
    g.cov_table <- convert_predictdb_cov_to_cor(g.cov_table)
  }

  # use the "rsid" column in g.weight_table if cov_table uses rsIDs
  if (any(grepl("rs", g.cov_table$RSID1))) {
    g.weight_table$varID <- g.weight_table$rsid
  }

  if (!all(g.weight_table$varID %in% unique(c(g.cov_table$RSID1, g.cov_table$RSID2)))){
    stop("Not all variants in weight_table are in cov_table!")
  }

  # convert correlation table to LD matrix
  varIDs <- g.weight_table$varID
  n_wgt <- length(varIDs)

  if (n_wgt == 0) {
    return(NULL)
  }

  R_wgt <- diag(n_wgt)

  if (n_wgt > 1) {
    snp_pairs <- combn(length(varIDs), 2)
    R_snp_pairs <- apply(snp_pairs, 2, function(x){
      match.idx <- which(g.cov_table$RSID1 == varIDs[x[1]] & g.cov_table$RSID2 == varIDs[x[2]])
      if (length(match.idx) == 0){
        match.idx <- which(g.cov_table$RSID1 == varIDs[x[2]] & g.cov_table$RSID2 == varIDs[x[1]])
      }
      if (length(match.idx) > 0) {
        g.cov_table[match.idx, "VALUE"]
      } else {
        NA
      }
    })
    R_wgt[t(snp_pairs)] <- R_snp_pairs
    R_wgt[t(snp_pairs[c(2,1),])] <- R_snp_pairs
  }

  rownames(R_wgt) <- g.weight_table$rsid
  colnames(R_wgt) <- g.weight_table$rsid

  return(R_wgt)
}

# Converts predictDB covariance to correlation
#' @importFrom stats setNames
convert_predictdb_cov_to_cor <- function(cov_table){
  stdev_table <- cov_table[cov_table$RSID1==cov_table$RSID2,]
  stdev_table <- setNames(sqrt(stdev_table$VALUE), stdev_table$RSID1)
  cov_table$VALUE <- cov_table$VALUE/(stdev_table[cov_table$RSID1]*stdev_table[cov_table$RSID2])

  return(cov_table)
}

#' @title Computes LD for weight variants using reference LD
#'
#' @param weights a list of preprocessed weights.
#'
#' @param region_info a data frame of region definitions.
#'
#' @param LD_map a data frame with filenames of LD matrices and SNP information for the regions.
#' Required when \code{load_predictdb_LD = FALSE}.
#'
#' @param snp_map a list of SNP-to-region map for the reference.
#' If NUll, it will reads SNP info from the "SNP_file" column of LD_map.
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
#' @importFrom parallel mclapply
#' @importFrom Matrix bdiag
#' @importFrom logging loginfo
#'
#' @return a list of processed weights, with LD of weight variants included.
#'
#' @export
compute_weight_LD_from_ref <- function(weights,
                                       region_info,
                                       LD_map,
                                       snp_map = NULL,
                                       LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
                                       LD_loader_fun = NULL,
                                       snpinfo_loader_fun = NULL,
                                       ncore = 1) {

  LD_format <- match.arg(LD_format)

  if (!inherits(weights,"list"))
    stop("'weights' should be a list!")

  if (!inherits(LD_map,"data.frame"))
    stop("'LD_map' should be a data frame!")

  weight_info <- lapply(names(weights), function(x){
    as.data.frame(weights[[x]][c("chrom", "p0","p1", "molecular_id", "weight_name", "type","context")])})
  weight_info <- do.call(rbind, weight_info)
  weight_info$weight_id <- paste0(weight_info$molecular_id, "|", weight_info$weight_name)
  # get the regions overlapping with each gene
  for (k in 1:nrow(weight_info)) {
    chrom <- weight_info[k, "chrom"]
    p0 <- weight_info[k, "p0"]
    p1 <- weight_info[k, "p1"]
    idx <- which(region_info$chrom == chrom & region_info$start <= p1 & region_info$stop > p0)
    weight_info[k, "region_id"] <- paste(sort(region_info[idx, "region_id"]), collapse = ",")
  }

  # compute LD for weight variants on each chromosome
  chrs <- sort(unique(weight_info$chrom))
  for (b in chrs) {
    loginfo("Computing LD for variants in weights on chr%s", b)
    weightinfo <- weight_info[weight_info$chrom == b, ]
    if (nrow(weightinfo) > 0) {
      weight_region_ids <- names(sort(-table(weightinfo$region_id)))
      weight_LD_list <- mclapply_check(weight_region_ids, function(x){
        # load the R_snp and SNP info for the region
        # and extract LD for the weight variants
        curr_region_LD_list <- list()
        curr_region_ids <- unlist(strsplit(x, ","))
        curr_region_idx <- match(curr_region_ids, LD_map$region_id)

        LD_matrix_files <- unlist(strsplit(LD_map$LD_file[curr_region_idx], split = ","))
        stopifnot(all(file.exists(LD_matrix_files)))

        if (length(LD_matrix_files) > 1) {
          R_snp <- lapply(LD_matrix_files, load_LD, format = LD_format, LD_loader_fun = LD_loader_fun)
          R_snp <- suppressWarnings(as.matrix(bdiag(R_snp)))
        } else {
          R_snp <- load_LD(LD_matrix_files, format = LD_format, LD_loader_fun = LD_loader_fun)
        }

        # load SNP info of the region
        # if snp_map is available, reads SNP info from snp_map;
        # otherwise, reads SNP info from the "SNP_file" column of LD_map.
        if (!is.null(snp_map)){
          snpinfo <- as.data.frame(rbindlist(snp_map[curr_region_ids], idcol = "region_id"))
        } else {
          SNP_info_files <- LD_map$SNP_file[curr_region_idx]
          stopifnot(all(file.exists(SNP_info_files)))
          snpinfo <- read_snp_info_files(SNP_info_files, snpinfo_loader_fun = snpinfo_loader_fun)
        }

        rownames(R_snp) <- snpinfo$id
        colnames(R_snp) <- snpinfo$id
        weight_ids <- weightinfo[weightinfo$region_id == x, "weight_id"]

        for (weight_id in weight_ids) {
          wgt_snp_ids <- rownames(weights[[weight_id]]$wgt)
          R_wgt <- R_snp[wgt_snp_ids, wgt_snp_ids, drop=FALSE]
          curr_region_LD_list[[weight_id]] <- R_wgt
        }
        curr_region_LD_list
      }, mc.cores = ncore, stop_if_missing = TRUE)

      weight_LD_list <- unlist(weight_LD_list, recursive = FALSE)
      for(weight_id in names(weight_LD_list)){
        weights[[weight_id]][["R_wgt"]] <- weight_LD_list[[weight_id]]
      }
    }
  }
  return(weights)
}

#' @title Get genome build of PredictDB weight
#'
#' @param weight_file a string, pointing path to weights in PredictDB format.
#'
#' @importFrom RSQLite dbDriver dbConnect dbGetQuery dbDisconnect
#'
#' @export
#'
get_predictdb_genome_build <- function(weight_file){

  # read the PredictDB weights
  stopifnot(file.exists(weight_file))
  sqlite <- dbDriver("SQLite")
  db <- dbConnect(sqlite, weight_file)
  query <- function(...) dbGetQuery(db, ...)
  weight_table <- query("select * from weights")
  varID_genomebuild <- unique(sapply(strsplit(weight_table$varID, split = "_"), "[[", 5))
  dbDisconnect(db)

  return(varID_genomebuild)
}

