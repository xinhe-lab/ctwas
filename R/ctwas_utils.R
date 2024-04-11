#' Function to identify consecutive regions and label them
label_merged_regions <- function(df) {
  # Identify consecutive regions
  df$merged <- c(FALSE, df$start[-1] == df$stop[-length(df$stop)])
  df$label <- cumsum(!df$merged)

  # Remove the consecutive column as it's no longer needed
  df$merged <- NULL

  return(df)
}

#' read all LD SNP info files as a data frame
read_LD_SNP_files <- function(files){
  if (length(files)>0){
    LD_SNP <- do.call(rbind, lapply(files, read_LD_SNP_file))
  } else {
    LD_SNP <- data.table::data.table(chrom=as.integer(), id=as.character(),
                                     pos=as.integer(), alt=as.character(),
                                     ref=as.character(), variance=as.numeric())
  }
  LD_SNP
}

#' read a single LD SNP info file as a data frame
read_LD_SNP_file <- function(file){
  LD_SNP <- data.table::fread(file, header = T)
  target_header <- c("chrom", "id", "pos", "alt", "ref")

  if (!all(target_header %in% colnames(LD_SNP))){
    stop("The .Rvar file needs to contain the following columns: ",
         paste(target_header, collapse = " "))
  }
  if (length(unique(LD_SNP$chrom)) != 1){
    stop("LD region needs to be on only one chromosome.")
  }

  return(LD_SNP)
}

#' Load PredictDB or Fusion weights
#'
#' @param weight_file a string or a vector, pointing path to one or multiple sets of weights in PredictDB or Fusion format.
#'
#' @param weight_format a string or a vector, specifying format of each weight file, e.g. PredictDB, Fusion.
#'
#' @param filter_protein_coding_genes TRUE/FALSE. If TRUE, keep protein coding genes only. This option is only for PredictDB weights
#'
#' @param load_predictdb_LD TRUE/FALSE. If TRUE, load pre-computed LD among weight SNPs. This option is only for PredictDB weights
#'
#' @param method_Fusion a string, specifying the method to choose in Fusion models
#'
#' @param ncore integer, number of cores for parallel computing.
#'
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble
#' @importFrom dplyr left_join mutate select
#'
#' @export
#'
load_weights <- function(weight_file,
                         weight_format = c("PredictDB", "Fusion"),
                         filter_protein_coding_genes = FALSE,
                         load_predictdb_LD = FALSE,
                         method_Fusion = c("lasso","enet","top1","blup"),
                         ncore = 1){
  #library(tidyverse)
  weight_format <- match.arg(weight_format)
  method_Fusion <- match.arg(method_Fusion)
  if(weight_format == "PredictDB"){
    weight_name <- tools::file_path_sans_ext(basename(weight_file))
    # read the PredictDB weights
    sqlite <- RSQLite::dbDriver("SQLite")
    db = RSQLite::dbConnect(sqlite, weight_file)
    query <- function(...) RSQLite::dbGetQuery(db, ...)
    weight_table <- query("select * from weights")
    weight_table <- weight_table[weight_table$weight!=0,]
    extra_table <- query("select * from extra")

    # subset to protein coding genes only
    if (isTRUE(filter_protein_coding_genes)){
      if ("protein_coding" %in% extra_table$gene_type){ #filter only there exist protein coding genes
        loginfo("Keep protein coding genes only")
        extra_table <- extra_table[extra_table$gene_type=="protein_coding",,drop=F]
        weight_table <- weight_table[weight_table$gene %in% extra_table$gene,]
      }
    }

    if (isTRUE(load_predictdb_LD)){
      R_wgt <- read.table(gzfile(paste0(tools::file_path_sans_ext(weight_file), ".txt.gz")), header=T)
    }
    else{
      R_wgt <- NULL
    }
    RSQLite::dbDisconnect(db)
  }
  else if (weight_format == "Fusion"){
    weight_name <- tools::file_path_sans_ext(basename(weight_file))
    wgtdir <- dirname(weight_file)
    wgtposfile <- file.path(wgtdir, paste0(basename(weight_file), ".pos"))
    wgtpos <- read.table(wgtposfile, header = T, stringsAsFactors = F)
    wgtpos <- transform(wgtpos, ID = ifelse(duplicated(ID) | duplicated(ID, fromLast = TRUE),
                                            paste(ID, ave(ID, ID, FUN = seq_along), sep = "_ID"), ID))
    wgtpos <- wgtpos[wgtpos$ID!="NA_IDNA",] #filter NA genes
    loginfo("Loading FUSION weights ...")
    cl <- parallel::makeCluster(ncore, outfile = "", type = "FORK")
    doParallel::registerDoParallel(cl)
    weight_table <- NULL
    if (nrow(wgtpos) > 0) {
      weight_table <- foreach(i = 1:nrow(wgtpos), .combine = "rbind") %dopar% {
        wf <- file.path(wgtdir, wgtpos[i, "WGT"])
        load(wf)
        gname <- wgtpos[i, "ID"]
        colnames(snps) <- c("chrom", "rsid", "cm", "pos", "alt", "ref")
        snps[is.na(snps$rsid),"rsid"] <- paste0("chr",snps[is.na(snps$rsid),"chrom"],"_",snps[is.na(snps$rsid),"pos"],
                                                "_",snps[is.na(snps$rsid),"ref"], "_",snps[is.na(snps$rsid),"alt"],"_","b38")

        snps[,"varID"] <- paste0("chr",snps[,"chrom"],"_",snps[,"pos"],"_",snps[,"ref"],"_",snps[,"alt"],"_","b38")

        rownames(wgt.matrix) <- snps$rsid
        g.method = method_Fusion
        # Ensure only top magnitude snp weight in the top1 wgt.matrix column
        if (g.method == "top1"){
          wgt.matrix[,"top1"][-which.max(wgt.matrix[,"top1"]^2)] <- 0
        }

        wgt.matrix <- wgt.matrix[abs(wgt.matrix[, g.method]) > 0, , drop = F]
        wgt.matrix <- wgt.matrix[complete.cases(wgt.matrix), , drop = F]
        if (nrow(wgt.matrix) > 0){
          out_table <- as_tibble(wgt.matrix[,g.method]) %>%
            mutate(rsid = rownames(wgt.matrix)) %>%
            left_join(tibble::as_tibble(snps) %>% select(-cm), by = "rsid")
          out_table$gene <- gname
          out_table <- out_table[,c("gene","rsid","varID","ref","alt","value")]
          colnames(out_table) <- c("gene","rsid","varID","ref_allele","eff_allele","weight")
          out_table
        }
      }
      parallel::stopCluster(cl)
    }
    extra_table <- NULL
    R_wgt <- NULL
  }
  return(list(weight_table=weight_table,extra_table=extra_table,weight_name=weight_name,R_wgt=R_wgt))
}


get_weight_LD <- function(R_wgt_all, gname, rsid_varID){
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


#' Load LD matrix
load_LD <- function(file, format = c("rds", "rdata", "csv", "txt", "tsv")) {
  format <- match.arg(format)

  # if format is missing, try to guess format by file extension
  if (missing(format)) {
    file_ext_lower <- tolower(tools::file_ext(file))

    if (file_ext_lower == "rds"){
      format <- "rds"
    } else if (file_ext_lower %in% c("rdata", "rd", "rda", "rdat")){
      format <- "rdata"
    } else if (file_ext_upper %in% c("csv", "csv.gz")) {
      format <- "csv"
    } else if (file_ext_lower %in% c("txt", "txt.gz")){
      format <- "txt"
    } else if (file_ext_lower %in% c("tsv", "tsv.gz")){
      format <- "tsv"
    } else {
      stop("Unknown LD file format!")
    }
  }

  if (format == "rds"){
    res <- readRDS(file)
  } else if (format == "rdata"){
    res <- get(load(file))
  } else if (format == "csv"){
    res <- as.matrix(read.csv(file, sep=",", row.names=1))
  } else if (format %in% c("txt", "tsv")){
    res <- as.matrix(data.table::fread(file))
  } else {
    stop("Unknown file format!")
  }

  return(res)
}

