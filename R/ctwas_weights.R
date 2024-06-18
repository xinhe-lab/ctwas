
#' @title Loads weights in PredictDB or FUSION format
#'
#' @param weight_file a string or a vector, pointing path to one or multiple sets of weights in PredictDB or FUSION format.
#'
#' @param weight_format a string or a vector, specifying format of each weight file, e.g. PredictDB, FUSION.
#'
#' @param filter_protein_coding_genes TRUE/FALSE. If TRUE, keep protein coding genes only. This option is only for PredictDB weights
#'
#' @param load_predictdb_LD TRUE/FALSE. If TRUE, load pre-computed LD among weight SNPs. This option is only for PredictDB weights
#'
#' @param fusion_method a string, specifying the method to choose in FUSION models
#'
#' @param fusion_genome_version a string, specifying the genome version of FUSION models
#'
#' @param ncore integer, number of cores for parallel computing.
#'
#' @export
#'
load_weights <- function(weight_file,
                         weight_format = c("PredictDB", "FUSION"),
                         filter_protein_coding_genes = FALSE,
                         load_predictdb_LD = FALSE,
                         fusion_method = c("lasso","enet","top1","blup"),
                         fusion_genome_version = c("b38","b37"),
                         ncore = 1){

  weight_format <- match.arg(weight_format)
  fusion_method <- match.arg(fusion_method)
  fusion_genome_version <- match.arg(fusion_genome_version)
  stopifnot(file.exists(weight_file))

  if (weight_format == "PredictDB") {
    res <- load_predictdb_weights(
      weight_file,
      filter_protein_coding_genes = filter_protein_coding_genes,
      load_predictdb_LD = load_predictdb_LD)
  } else if (weight_format == "FUSION") {
    res <- load_fusion_weights(
      weight_file,
      fusion_method = fusion_method,
      fusion_genome_version = fusion_genome_version,
      ncore = ncore)
  }
  return(res)
}

#' @title Loads weights in PredictDB format
#'
#' @param weight_file a string or a vector, pointing path to one or multiple
#' sets of weights in PredictDB format.
#'
#' @param filter_protein_coding_genes TRUE/FALSE. If TRUE, keep protein coding
#' genes only.
#'
#' @param load_predictdb_LD TRUE/FALSE. If TRUE, load pre-computed LD
#' among weight variants.
#'
#' @importFrom tools file_path_sans_ext
#' @importFrom RSQLite dbDriver dbConnect dbGetQuery dbDisconnect
#'
#' @export
#'
load_predictdb_weights <- function(weight_file,
                                   filter_protein_coding_genes = FALSE,
                                   load_predictdb_LD = FALSE){

  # read the PredictDB weights
  stopifnot(file.exists(weight_file))
  weight_name <- file_path_sans_ext(basename(weight_file))
  sqlite <- dbDriver("SQLite")
  db <- dbConnect(sqlite, weight_file)
  query <- function(...) dbGetQuery(db, ...)
  weight_table <- query("select * from weights")
  weight_table <- weight_table[weight_table$weight!=0,]
  extra_table <- query("select * from extra")

  # subset to protein coding genes only
  if (filter_protein_coding_genes) {
    # filter only when protein coding genes exist
    if ("protein_coding" %in% extra_table$gene_type){
      loginfo("Keep protein coding genes only")
      extra_table <- extra_table[extra_table$gene_type=="protein_coding",,drop=F]
      weight_table <- weight_table[weight_table$gene %in% extra_table$gene,]
    }
  }

  # load pre-computed LD from PredictDB LD file
  if (isTRUE(load_predictdb_LD)){
    predictdb_LD_file <- paste0(file_path_sans_ext(weight_file), ".txt.gz")
    if (!file.exists(predictdb_LD_file)){
      stop(paste("PredictDB LD file", predictdb_LD_file, "does not exist!"))
    }
    R_wgt <- read.table(gzfile(predictdb_LD_file), header=T)
  }
  else{
    R_wgt <- NULL
  }
  dbDisconnect(db)

  return(list(weight_table=weight_table,
              extra_table=extra_table,
              weight_name=weight_name,
              R_wgt=R_wgt))
}

#' @title Loads weights in FUSION format
#'
#' @param weight_file a string or a vector, pointing path to one or multiple
#' sets of weights in FUSION format.
#'
#' @param fusion_method a string, specifying the method to choose in FUSION models.
#'
#' @param fusion_genome_version a string, specifying the genome version of FUSION models
#'
#' @param ncore integer, number of cores for parallel computing.
#'
#' @importFrom utils read.table
#' @importFrom stats complete.cases
#' @importFrom stats ave
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble
#' @importFrom dplyr left_join mutate select
#' @importFrom tools file_path_sans_ext
#' @importFrom parallel mclapply
#'
#' @export
#'
load_fusion_weights <- function(weight_file,
                                fusion_method = c("lasso","enet","top1","blup"),
                                fusion_genome_version = c("b38","b37"),
                                ncore = 1) {

  fusion_method <- match.arg(fusion_method)
  fusion_genome_version <- match.arg(fusion_genome_version)

  stopifnot(file.exists(weight_file))
  weight_name <- file_path_sans_ext(basename(weight_file))
  wgt_dir <- dirname(weight_file)
  wgt_pos_file <- file.path(wgt_dir, paste0(basename(weight_file), ".pos"))
  wgt_pos <- read.table(wgt_pos_file, header = T, stringsAsFactors = F)
  wgt_pos$ID <-
    ifelse(duplicated(wgt_pos$ID) | duplicated(wgt_pos$ID,fromLast = TRUE),
           paste(wgt_pos$ID,ave(wgt_pos$ID,wgt_pos$ID,FUN = seq_along),
                 sep = "_ID"),
           wgt_pos$ID)
  wgt_pos <- wgt_pos[wgt_pos$ID!="NA_IDNA",] # filter NA genes

  loginfo("Loading FUSION weights ...")
  weight_table <- NULL
  if (nrow(wgt_pos) > 0) {
    weight_table_list <- mclapply(1:nrow(wgt_pos), function(i){
      load(file.path(wgt_dir, wgt_pos[i, "WGT"]))
      gname <- wgt_pos[i, "ID"]
      colnames(snps) <- c("chrom", "rsid", "cm", "pos", "alt", "ref")
      snps[is.na(snps$rsid),"rsid"] <- paste0("chr",snps[is.na(snps$rsid),"chrom"],
                                              "_",snps[is.na(snps$rsid),"pos"],
                                              "_",snps[is.na(snps$rsid),"ref"],
                                              "_",snps[is.na(snps$rsid),"alt"],
                                              "_",fusion_genome_version)

      snps[,"varID"] <- paste0("chr",snps[,"chrom"],"_",snps[,"pos"],"_",
                               snps[,"ref"],"_",snps[,"alt"],
                               "_",fusion_genome_version)

      rownames(wgt.matrix) <- snps$rsid
      g.method <- fusion_method
      # Ensure only top magnitude snp weight in the top1 wgt.matrix column
      if (g.method == "top1"){
        wgt.matrix[,"top1"][-which.max(wgt.matrix[,"top1"]^2)] <- 0
      }
      wgt.matrix <- wgt.matrix[abs(wgt.matrix[, g.method]) > 0, , drop = FALSE]
      wgt.matrix <- wgt.matrix[complete.cases(wgt.matrix), , drop = FALSE]
      if (nrow(wgt.matrix) > 0){
        out_table <- as_tibble(wgt.matrix[,g.method]) %>%
          mutate(rsid = rownames(wgt.matrix)) %>%
          left_join(as_tibble(snps) %>% select(-"cm"), by = "rsid")
        out_table$gene <- gname
        out_table <- out_table[,c("gene","rsid","varID","ref","alt","value")]
        colnames(out_table) <- c("gene","rsid","varID","ref_allele","eff_allele","weight")
        out_table
      }
    })
    if (length(weight_table_list) != nrow(wgt_pos)) {
      stop("Not all cores returned results. Try rerun with bigger memory or fewer cores")
    }
    weight_table <- do.call(rbind, weight_table_list)
  }
  extra_table <- NULL
  R_wgt <- NULL

  return(list(weight_table=weight_table,
              extra_table=extra_table,
              weight_name=weight_name,
              R_wgt=R_wgt))
}

# Gets pre-computed LD matrix from predictedDB weights
#' @importFrom stats setNames
get_weight_LD <- function (R_wgt_all, gname, rsid_varID){
  R_wgt <- R_wgt_all[R_wgt_all$GENE == gname,]
  #convert covariance to correlation
  R_wgt_stdev <- R_wgt[R_wgt$RSID1==R_wgt$RSID2,]
  R_wgt_stdev <- setNames(sqrt(R_wgt_stdev$VALUE), R_wgt_stdev$RSID1)
  R_wgt$VALUE <- R_wgt$VALUE/(R_wgt_stdev[R_wgt$RSID1]*R_wgt_stdev[R_wgt$RSID2])

  unique_id <- unique(c(R_wgt$RSID1, R_wgt$RSID2))

  # Create an empty correlation matrix
  n <- length(unique_id)
  cor_matrix <- matrix(NA, nrow = n, ncol = n)

  # Fill in the correlation values
  for (i in 1:n) {
    for (j in i:n) {  # Only iterate over half of the matrix
      if (i == j) {
        cor_matrix[i, j] <- 1  # Diagonal elements are 1
      } else {
        # Check if there are any matches for the RSID combination
        matches <- R_wgt[R_wgt$RSID1 == unique_id[i] & R_wgt$RSID2 == unique_id[j], "VALUE"]
        if (length(matches) > 0) {
          cor_matrix[i, j] <- matches
          cor_matrix[j, i] <- matches  # Set symmetric value
        } else {
          cor_matrix[i, j] <- NA  # No correlation value found
          cor_matrix[j, i] <- NA  # No correlation value found
        }
      }
    }
  }

  rownames(cor_matrix) <- rsid_varID$rsid[match(unique_id, rsid_varID$varID)]
  colnames(cor_matrix) <- rsid_varID$rsid[match(unique_id, rsid_varID$varID)]

  return(cor_matrix)
}

# Computes LD for weight variants using reference LD
#' @importFrom parallel mclapply
#' @importFrom Matrix bdiag
#' @importFrom logging loginfo
compute_weight_LD_from_ref <- function(weights,
                                       weight_name,
                                       region_info,
                                       LD_info,
                                       snp_info,
                                       LD_format = c("rds", "rdata", "csv", "txt", "custom"),
                                       LD_loader = NULL,
                                       ncore = 1) {

  if (is.null(LD_info) || is.null(snp_info)) {
    stop("LD_info and snp_info are required for computing LD")
  }

  LD_format <- match.arg(LD_format)

  weight_info <- lapply(names(weights), function(x){
    as.data.frame(weights[[x]][c("chrom", "p0","p1", "gene_name", "weight_name", "type","context")])})
  weight_info <- do.call(rbind, weight_info)
  weight_info$weight_id <- paste0(weight_info$gene_name, "|", weight_name)
  # get the regions overlapping with each gene
  for (k in 1:nrow(weight_info)) {
    chrom <- weight_info[k, "chrom"]
    p0 <- weight_info[k, "p0"]
    p1 <- weight_info[k, "p1"]
    idx <- which(region_info$chrom == chrom & region_info$start <= p1 & region_info$stop > p0)
    weight_info[k, "region_id"] <- paste(sort(region_info[idx, "region_id"]), collapse = ";")
  }

  # compute LD for weight variants on each chromosome
  chrs <- sort(unique(weight_info$chrom))
  for (b in chrs) {
    loginfo("Computing LD for weight variants on chr%s", b)
    weightinfo <- weight_info[weight_info$chrom == b, ]
    if (nrow(weightinfo) > 0) {
      weight_region_ids <- names(sort(-table(weightinfo$region_id)))
      weight_LD_list <- mclapply(weight_region_ids, function(x){
        # load the R_snp and SNP info for the region
        # and extract LD for the weight variants
        curr_region_LD_list <- list()
        curr_region_ids <- unlist(strsplit(x, ";"))
        curr_region_idx <- match(curr_region_ids, LD_info$region_id)
        LD_matrix_files <- LD_info$LD_matrix[curr_region_idx]
        if (length(LD_matrix_files) > 1) {
          R_snp <- lapply(LD_matrix_files, load_LD, format = LD_format, LD_loader = LD_loader)
          R_snp <- suppressWarnings(as.matrix(bdiag(R_snp)))
        } else {
          R_snp <- load_LD(LD_matrix_files, format = LD_format, LD_loader = LD_loader)
        }

        snpinfo <- do.call(rbind, snp_info[curr_region_ids])
        rownames(R_snp) <- snpinfo$id
        colnames(R_snp) <- snpinfo$id

        weight_ids <- weightinfo[weightinfo$region_id == x, "weight_id"]

        for (weight_id in weight_ids) {
          snpnames <- rownames(weights[[weight_id]]$wgt)
          R_wgt <- R_snp[snpnames, snpnames, drop=F]
          curr_region_LD_list[[weight_id]] <- R_wgt
        }
        curr_region_LD_list
      })
      if (length(weight_LD_list) != length(weight_region_ids)) {
        stop("Not all cores returned results. Try rerun with bigger memory or fewer cores")
      }
      weight_LD_list <- unlist(weight_LD_list, recursive = FALSE)
      for(weight_id in names(weight_LD_list)){
        weights[[weight_id]][["R_wgt"]] <- weight_LD_list[[weight_id]]
      }
    }
  }
  return(weights)
}

#' @title Makes PredictDB weights from top QTLs
#'
#' @description Make PredictDB weights from only top QTLs.
#' Each gene (or molecular trait) has only one top QTL, defined by users,
#' based on min p-value, max abs(weight), etc.
#'
#' @param QTL_data a data frame with required columns:
#' "gene", "varID", "weight".
#'
#' @param gene_info a data frame (optional) with information of the genes
#' ("gene","genename","gene_type", etc.)
#'
#' @param output_dir output directory
#'
#' @param weight_name name of the output weight files
#'
#' @importFrom utils read.table write.table
#' @importFrom stats complete.cases
#' @importFrom RSQLite dbDriver dbConnect dbWriteTable dbDisconnect
#'
#' @export
#'
make_predictdb_weights_from_top_QTLs <- function(QTL_data,
                                                 gene_info=NULL,
                                                 output_dir = getwd(),
                                                 weight_name){

  if (!dir.exists(output_dir))
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # check and clean the QTL data
  required_cols <- c("gene", "varID", "weight")
  if (!all(required_cols %in% colnames(QTL_data))){
    stop("QTL_data needs to contain the following columns: ",
         paste(required_cols, collapse = " "))
  }
  QTL_data <- QTL_data[abs(QTL_data[, "weight"]) > 0, ,drop = FALSE]
  # remove missing values
  QTL_data <- QTL_data[complete.cases(QTL_data), ,drop = FALSE]

  # Create a database connection
  driver <- dbDriver('SQLite')
  db <- dbConnect(drv = driver, file.path(output_dir, paste0(weight_name,".db")))
  # create weights table
  dbWriteTable(db, 'weights', QTL_data, overwrite = TRUE)
  # create extra table with gene info
  if (is.null(gene_info)) {
    gene_info <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(gene_info) <- c("gene","genename","gene_type")
  }
  dbWriteTable(db, 'extra', gene_info, overwrite = TRUE)
  dbDisconnect(db)

  # set covariance to 1 as each gene only has one top QTL
  covar_data <- QTL_data[,c("gene","varID","varID")]
  colnames(covar_data) <- c("GENE","RSID1","RSID2")
  covar_data$VALUE <- 1
  write.table(covar_data,
              file = gzfile(file.path(output_dir,paste0(weight_name,".txt.gz"))),
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
}
