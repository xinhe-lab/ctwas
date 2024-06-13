
# Read a single SNP info file as a data frame
#' @importFrom data.table fread
read_snp_info_file <- function (file){
  snp_info <- as.data.frame(fread(file, header = TRUE))
  target_header <- c("chrom", "id", "pos", "alt", "ref")
  if (!all(target_header %in% colnames(snp_info))){
    stop("The SNP info file needs to contain the following columns: ",
         paste(target_header, collapse = " "))
  }
  return(snp_info)
}

# Read all SNP info files as a data frame
read_snp_info_files <- function (files){
  snp_info <- do.call(rbind, lapply(files, read_snp_info_file))
  snp_info <- unique(as.data.frame(snp_info))
  return(snp_info)
}

#' @title Loads PredictDB or FUSION weights
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
#' @importFrom RSQLite dbDriver dbConnect dbGetQuery dbDisconnect
#' @importFrom parallel mclapply
#'
#' @export
#'
load_weights <- function(weight_file,
                         weight_format = c("PredictDB", "FUSION"),
                         filter_protein_coding_genes = FALSE,
                         load_predictdb_LD = FALSE,
                         method_FUSION = c("lasso","enet","top1","blup"),
                         fusion_genome_version = c("b38","b37"),
                         ncore = 1){

  weight_format <- match.arg(weight_format)
  method_FUSION <- match.arg(method_FUSION)
  fusion_genome_version <- match.arg(fusion_genome_version)
  stopifnot(file.exists(weight_file))

  if(weight_format == "PredictDB"){

    weight_name <- file_path_sans_ext(basename(weight_file))
    # read the PredictDB weights
    sqlite <- dbDriver("SQLite")
    db <- dbConnect(sqlite, weight_file)
    query <- function(...) dbGetQuery(db, ...)
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
  }
  else if (weight_format == "FUSION"){
    weight_name <- file_path_sans_ext(basename(weight_file))
    wgtdir <- dirname(weight_file)
    wgtposfile <- file.path(wgtdir, paste0(basename(weight_file), ".pos"))
    wgtpos <- read.table(wgtposfile, header = T, stringsAsFactors = F)
    wgtpos$ID <-
      ifelse(duplicated(wgtpos$ID) | duplicated(wgtpos$ID,fromLast = TRUE),
             paste(wgtpos$ID,ave(wgtpos$ID,wgtpos$ID,FUN = seq_along),
                   sep = "_ID"),
             ID)
    wgtpos <- wgtpos[wgtpos$ID!="NA_IDNA",] #filter NA genes
    loginfo("Loading FUSION weights ...")
    weight_table <- NULL
    if (nrow(wgtpos) > 0) {
      weight_table_list <- mclapply(1:nrow(wgtpos), function(i){
        load(file.path(wgtdir, wgtpos[i, "WGT"]))
        gname <- wgtpos[i, "ID"]
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
        g.method <- method_FUSION
        # Ensure only top magnitude snp weight in the top1 wgt.matrix column
        if (g.method == "top1"){
          wgt.matrix[,"top1"][-which.max(wgt.matrix[,"top1"]^2)] <- 0
        }
        wgt.matrix <- wgt.matrix[abs(wgt.matrix[, g.method]) > 0, , drop = FALSE]
        wgt.matrix <- wgt.matrix[complete.cases(wgt.matrix), , drop = FALSE]
        if (nrow(wgt.matrix) > 0){
          out_table <- as_tibble(wgt.matrix[,g.method]) %>%
            mutate(rsid = rownames(wgt.matrix)) %>%
            left_join(as_tibble(snps) %>% select(-cm), by = "rsid")
          out_table$gene <- gname
          out_table <- out_table[,c("gene","rsid","varID","ref","alt","value")]
          colnames(out_table) <- c("gene","rsid","varID","ref_allele","eff_allele","weight")
          out_table
        }
      })
      if (length(weight_table_list) != nrow(wgtpos)) {
        stop("Not all cores returned results. Try rerun with bigger memory or fewer cores")
      }
      weight_table <- do.call(rbind, weight_table_list)
    }
    extra_table <- NULL
    R_wgt <- NULL
  }
  return(list(weight_table=weight_table,
              extra_table=extra_table,
              weight_name=weight_name,
              R_wgt=R_wgt))
}

# Load LD matrix
# @param file path to LD matrix
# @param format file format for LD matrix
# @param LD_loader a user defined function to load LD matrix
#' @importFrom utils read.csv
#' @importFrom data.table fread
#' @importFrom tools file_ext
load_LD <- function (file, format = c("rds", "rdata", "csv", "txt", "custom"), LD_loader = NULL) {
  format <- match.arg(format)

  # if format is missing, try to guess format by file extension
  if (missing(format)) {
    file_ext_lower <- tolower(file_ext(file))
    if (file_ext_lower == "rds"){
      format <- "rds"
    } else if (file_ext_lower %in% c("rdata", "rd", "rda", "rdat")){
      format <- "rdata"
    } else if (file_ext_lower %in% c("csv", "csv.gz")) {
      format <- "csv"
    } else if (file_ext_lower %in% c("txt", "txt.gz")){
      format <- "txt"
    } else {
      # set format to "custom" if not matched with known formats
      format <- "custom"
    }
  }

  if (format == "rds"){
    R <- readRDS(file)
  } else if (format == "rdata"){
    R <- get(load(file))
  } else if (format == "csv"){
    R <- as.matrix(read.csv(file, sep=",", row.names=1))
  } else if (format == "txt"){
    R <- as.matrix(fread(file))
  } else if (format == "custom") {
    # use LD_loader() function to load LD matrix
    R <- LD_loader(file)
  }

  return(R)
}


# Prepare .pvar file
#
# @param pgenf pgen file
# .pvar file format: https://www.cog-genomics.org/plink/2.0/formats#pvar
#
# @param outputdir a string, the directory to store output
#
# @return corresponding pvar file
#
#' @importFrom utils read.table
#' @importFrom data.table fread fwrite
#' @importFrom tools file_ext file_path_sans_ext
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
        pvar <- fread(pvarf, header = F)

        if (ncol(pvar) == 6) {
          colnames(pvar) <- c('#CHROM', 'ID', 'CM', 'POS', 'ALT', 'REF')
        } else if (ncol(pvar) == 5){
          colnames(pvar) <- c('#CHROM', 'ID', 'POS', 'ALT', 'REF')
        } else {
          stop(".pvar file has incorrect format")
        }

        fwrite(pvar, file = pvarf2 , sep="\t", quote = F)
      }
    }

  } else if (file_ext(pgenf) == "bed"){
    # .bim file has no header
    pvarf <- paste0(file_path_sans_ext(pgenf), ".bim")
    pvarf2 <-  file.path(outputdir, paste0(basename(file_path_sans_ext(pgenf)), ".hbim"))

    if (!file.exists(pvarf2)){
      pvar <- fread(pvarf, header = F)
      colnames(pvar) <- c('#CHROM', 'ID', 'CM', 'POS', 'ALT', 'REF')
      fwrite(pvar, file = pvarf2 , sep="\t", quote = F)
    }
    pvarfout <- pvarf2
  } else {
    stop("Unrecognized genotype input format")
  }

  return(pvarfout)
}

# Read .pgen file into R
#
# @param pgenf .pgen file or .bed file
#
# @param pvarf .pvar file or .bim file with have proper
#  header.  Matching `pgenf`.
#
# @return  A matrix of allele count for each variant (columns) in each sample
#  (rows). ALT allele in pvar file is counted (A1 allele in .bim file is the ALT
#   allele).
#
#' @importFrom data.table fread
#' @importFrom pgenlibr NewPgen NewPvar
#' @importFrom tools file_ext file_path_sans_ext
prep_pgen <- function(pgenf, pvarf){

  pvar <- NewPvar(pvarf)

  if (file_ext(pgenf) == "pgen"){
    pgen <- NewPgen(pgenf, pvar = pvar)

  } else if (file_ext(pgenf) == "bed"){
    famf <- paste0(file_path_sans_ext(pgenf), ".fam")
    fam <- fread(famf, header = F)
    raw_s_ct <- nrow(fam)
    pgen <- NewPgen(pgenf, pvar = pvar, raw_sample_ct = raw_s_ct)

  } else{
    stop("unrecognized input")
  }

  return(pgen)
}

# Read pgen file into R
#
# @param pgen .pgen file or .bed file
#
# @param variantidx variant index. If NULL, all variants will be extracted.
#
# @return A matrix, columns are allele count for each SNP, rows are
#  for each sample.
#
#' @importFrom pgenlibr GetVariantCt
#' @importFrom pgenlibr ReadList
read_pgen <- function(pgen, variantidx = NULL, meanimpute = F ){
  if (is.null(variantidx)){
    variantidx <- 1: GetVariantCt(pgen)}

  ReadList(pgen,
           variant_subset = variantidx,
           meanimpute = meanimpute)
}


# Read .pvar file into R
# @param pvarf .pvar file with proper format: https://www.cog-genomics.org/plink/2.0/formats#pvar
#
# @return A data.table. variant info
#' @importFrom data.table fread
#' @importFrom dplyr rename
read_pvar <- function(pvarf){
  pvar <- fread(pvarf, skip = "#CHROM")
  pvar <- rename(pvar, "chrom" = "#CHROM", "pos" = "POS",
                 "alt" = "ALT", "ref" = "REF", "id" = "ID")
  pvar <- pvar[, c("chrom", "id", "pos", "alt", "ref")]
  return(pvar)
}

# Read .bim file into R
# @param bimf .bim file with proper format: https://www.cog-genomics.org/plink/2.0/formats#bim
#
# @return A data.table. variant info
#' @importFrom data.table fread
read_bim <- function(bimf) {
  bim <- fread(bimf)
  colnames(bim) <- c("chr", "id", "cm", "pos", "alt", "ref")
  return(bim)
}

# Read variant information from .pvar or .bim file into R
# @param var_info_file .pvar or .bim file with proper format:
# .pvar: https://www.cog-genomics.org/plink/2.0/formats#pvar
# .bim: https://www.cog-genomics.org/plink/2.0/formats#bim
#
# @return A data.table. variant info
#' @importFrom tools file_ext
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

# assign regions to cores
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

