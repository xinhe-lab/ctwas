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

#' Load weight
load_weights <- function(weight_file, weight_format = c("PredictDB","Fusion"), filter_protein_coding_genes = FALSE, load_predictdb_LD = FALSE){
  weight_format <- match.arg(weight_format)
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
      loginfo("Keep protein coding genes only")
      extra_table <- extra_table[extra_table$gene_type=="protein_coding",,drop=F]
      weight_table <- weight_table[weight_table$gene %in% extra_table$gene,]
    }

    if (isTRUE(load_predictdb_LD)){
      R_wgt <- read.table(gzfile(paste0(tools::file_path_sans_ext(weight_file), ".txt.gz")), header=T)
    }
    else{
      R_wgt <- NULL
    }
    RSQLite::dbDisconnect(db)
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

