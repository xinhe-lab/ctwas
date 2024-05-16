
#' Read a single SNP info file as a data frame
read_snp_info_file <- function(file){
  snp_info <- as.data.frame(data.table::fread(file, header = TRUE))
  target_header <- c("chrom", "id", "pos", "alt", "ref")
  if (!all(target_header %in% colnames(snp_info))){
    stop("The SNP info file needs to contain the following columns: ",
         paste(target_header, collapse = " "))
  }
  return(snp_info)
}

#' Read all SNP info files as a data frame
read_snp_info_files <- function(files){
  snp_info <- do.call(rbind, lapply(files, read_snp_info_file))
  snp_info <- unique(as.data.frame(snp_info))
  return(snp_info)
}

#' Load PredictDB or FUSION weights
#'
#' @param weight_file a string or a vector, pointing path to one or multiple sets of weights in PredictDB or FUSION format.
#'
#' @param weight_format a string or a vector, specifying format of each weight file, e.g. PredictDB, FUSION.
#'
#' @param filter_protein_coding_genes TRUE/FALSE. If TRUE, keep protein coding genes only. This option is only for PredictDB weights
#'
#' @param load_predictdb_LD TRUE/FALSE. If TRUE, load pre-computed LD among weight SNPs. This option is only for PredictDB weights
#'
#' @param method_FUSION a string, specifying the method to choose in FUSION models
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
                         weight_format = c("PredictDB", "FUSION"),
                         filter_protein_coding_genes = FALSE,
                         load_predictdb_LD = FALSE,
                         method_FUSION = c("lasso","enet","top1","blup"),
                         ncore = 1){

  weight_format <- match.arg(weight_format)
  method_FUSION <- match.arg(method_FUSION)
  stopifnot(file.exists(weight_file))

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
    if (filter_protein_coding_genes) {
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
  else if (weight_format == "FUSION"){
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
        g.method <- method_FUSION
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


#' Prepare .pvar file
#'
#' @param pgenf pgen file
#' .pvar file format: https://www.cog-genomics.org/plink/2.0/formats#pvar
#'
#' @param outputdir a string, the directory to store output
#'
#' @return corresponding pvar file
#'
#' @importFrom tools file_ext file_path_sans_ext
#'
prep_pvar <- function(pgenf, outputdir = getwd()){

  if (file_ext(pgenf) == "pgen"){
    pvarf <- paste0(file_path_sans_ext(pgenf), ".pvar")
    pvarf2 <-  paste0(outputdir, basename(file_path_sans_ext(pgenf)), ".hpvar")

    # pgenlib can't read pvar without header, check if header present
    firstl <- read.table(file = pvarf, header = F, comment.char = '',
                         nrows = 1, stringsAsFactors = F)

    if (substr(firstl[1,1],1,1) == '#') {
      pvarfout <- pvarf
    } else {
      pvarfout <- pvarf2

      if (!file.exists(pvarf2)) {
        pvar <- data.table::fread(pvarf, header = F)

        if (ncol(pvar) == 6) {
          colnames(pvar) <- c('#CHROM', 'ID', 'CM', 'POS', 'ALT', 'REF')
        } else if (ncol(pvar) == 5){
          colnames(pvar) <- c('#CHROM', 'ID', 'POS', 'ALT', 'REF')
        } else {
          stop(".pvar file has incorrect format")
        }

        data.table::fwrite(pvar, file = pvarf2 , sep="\t", quote = F)
      }
    }

  } else if (file_ext(pgenf) == "bed"){
    # .bim file has no header
    pvarf <- paste0(file_path_sans_ext(pgenf), ".bim")
    pvarf2 <-  file.path(outputdir, paste0(basename(file_path_sans_ext(pgenf)), ".hbim"))

    if (!file.exists(pvarf2)){
      pvar <- data.table::fread(pvarf, header = F)
      colnames(pvar) <- c('#CHROM', 'ID', 'CM', 'POS', 'ALT', 'REF')
      data.table::fwrite(pvar, file = pvarf2 , sep="\t", quote = F)
    }
    pvarfout <- pvarf2
  } else {
    stop("Unrecognized genotype input format")
  }

  return(pvarfout)
}

#' Read .pgen file into R
#'
#' @param pgenf .pgen file or .bed file
#'
#' @param pvarf .pvar file or .bim file with have proper
#'  header.  Matching `pgenf`.
#'
#' @return  A matrix of allele count for each variant (columns) in each sample
#'  (rows). ALT allele in pvar file is counted (A1 allele in .bim file is the ALT
#'   allele).
#'
#' @importFrom pgenlibr NewPvar
#' @importFrom pgenlibr NewPgen
#' @importFrom tools file_ext file_path_sans_ext
#'
prep_pgen <- function(pgenf, pvarf){

  pvar <- pgenlibr::NewPvar(pvarf)

  if (file_ext(pgenf) == "pgen"){
    pgen <- pgenlibr::NewPgen(pgenf, pvar = pvar)

  } else if (file_ext(pgenf) == "bed"){
    famf <- paste0(file_path_sans_ext(pgenf), ".fam")
    fam <- data.table::fread(famf, header = F)
    raw_s_ct <- nrow(fam)
    pgen <- pgenlibr::NewPgen(pgenf, pvar = pvar, raw_sample_ct = raw_s_ct)

  } else{
    stop("unrecognized input")
  }

  return(pgen)
}

#' Read pgen file into R
#'
#' @param pgen .pgen file or .bed file
#'
#' @param variantidx variant index. If NULL, all variants will be extracted.
#'
#' @return A matrix, columns are allele count for each SNP, rows are
#'  for each sample.
#'
#' @importFrom pgenlibr GetVariantCt
#' @importFrom pgenlibr ReadList
#'
read_pgen <- function(pgen, variantidx = NULL, meanimpute = F ){
  if (is.null(variantidx)){
    variantidx <- 1: pgenlibr::GetVariantCt(pgen)}

  pgenlibr::ReadList(pgen,
                     variant_subset = variantidx,
                     meanimpute = meanimpute)
}


#' Read .pvar file into R
#' @param pvarf .pvar file with proper format: https://www.cog-genomics.org/plink/2.0/formats#pvar
#'
#' @return A data.table. variant info
#'
read_pvar <- function(pvarf){
  pvar <- data.table::fread(pvarf, skip = "#CHROM")
  pvar <- dplyr::rename(pvar, "chrom" = "#CHROM", "pos" = "POS",
                        "alt" = "ALT", "ref" = "REF", "id" = "ID")
  pvar <- pvar[, c("chrom", "id", "pos", "alt", "ref")]
  return(pvar)
}

#' Read .bim file into R
#' @param bimf .bim file with proper format: https://www.cog-genomics.org/plink/2.0/formats#bim
#'
#' @return A data.table. variant info
#'
read_bim <- function(bimf) {
  bim <- data.table::fread(bimf)
  colnames(bim) <- c("chr", "id", "cm", "pos", "alt", "ref")
  return(bim)
}

#' Read variant information from .pvar or .bim file into R
#' @param var_info_file .pvar or .bim file with proper format:
#' .pvar: https://www.cog-genomics.org/plink/2.0/formats#pvar
#' .bim: https://www.cog-genomics.org/plink/2.0/formats#bim
#'
#' @return A data.table. variant info
#' @importFrom tools file_ext
#'
read_var_info <- function(var_info_file){
  if (file_ext(var_info_file) == "pvar"){
    var_info <- read_pvar(var_info_file)
  } else if (file_ext(var_info_file) == "bim"){
    var_info <- read_bim(var_info_file)
  } else{
    stop("unrecognized input")
  }
  return(var_info)
}

#' assign regions to cores
region2core <- function(region_data, ncore = 1){
  region_ids <- names(region_data)
  if (ncore > 1) {
    d <- cut(1:length(region_ids), ncore, labels = FALSE)
    corelist <- split(region_ids,d)
  } else {
    corelist <- list("1" = region_ids)
  }
  return(corelist)
}

