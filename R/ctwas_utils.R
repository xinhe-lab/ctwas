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

# check SNP_info files from the region_info table, and return a list of region_info tables
# adapted from old write_ld_Rf() function
get_region_info_list <- function(region_info){

  # read and check SNP info (Rvar) data and extract min and max positions in each region
  for (i in 1:nrow(region_info)){
    region_snp_info <- read_LD_SNP_file(region_info$SNP_info[i])
    region_info$SNP_start[i] <- min(region_snp_info$pos)
    region_info$SNP_stop[i] <- max(region_snp_info$pos) + 1
  }

  region_info <- data.frame(region_info, stringsAsFactors = F)
  region_info <- transform(region_info,
                           chr = as.integer(chr),
                           SNP_start = as.integer(SNP_start),
                           SNP_stop = as.integer(SNP_stop))

  # save a data frame for each chromosome with region_name
  region_info_list <- vector("list", length = 22)

  for (b in sort(unique(region_info$chr))) {
    region_info_chr <- region_info[region_info$chr == b, , drop = F]
    region_info_chr <- region_info_chr[order(region_info_chr$start), ]
    region_info_chr$region_name <- 1:nrow(region_info_chr)
    region_info_chr$region_tag <- paste0(region_info_chr$start, "-", region_info_chr$stop)
    region_info_list[[b]] <- region_info_chr
  }

  return(region_info_list)
}


#' read LD by file format
read_LD <- function(file, format = c("rds", "rdata", "csv", "txt", "tsv")) {
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

#' read weights in PredictDB format, do not perform harmonization on weights
read_weights <- function (weight_files,
                          chr = 1:22,
                          ld_snpinfo,
                          z_snp = NULL,
                          scale_by_ld_variance=TRUE,
                          ncore=1){

  sqlite <- RSQLite::dbDriver("SQLite")

  # read gene names in each weight file
  gnames_all <- list()
  for (i in 1:length(weight_files)){
    weight_file <- weight_files[i]
    db <- RSQLite::dbConnect(sqlite, weight_file)
    query <- function(...) RSQLite::dbGetQuery(db, ...)
    gnames <- unique(query("select gene from weights")[, 1])
    gnames_all[[i]] <- cbind(gnames, weight_file)
    RSQLite::dbDisconnect(db)
  }
  gnames_all <- as.data.frame(do.call(rbind, gnames_all))
  colnames(gnames_all) <- c("gname", "weight")

  loginfo("Number of genes with weights provided: %s", nrow(gnames_all))
  loginfo("Collecting gene weight information ...")

  # parallelize genes into cores
  corelist <- lapply(1:ncore, function(core){
    njobs <- ceiling(nrow(gnames_all)/ncore);
    jobs <- ((core-1)*njobs+1):(core*njobs);
    jobs[jobs<=nrow(gnames_all)]})
  names(corelist) <- 1:ncore

  cl <- parallel::makeCluster(ncore, outfile = "")
  doParallel::registerDoParallel(cl)

  outlist <- foreach(core = 1:ncore, .combine = "c", .packages = "ctwas") %dopar% {
    gnames_core <- gnames_all[corelist[[core]],,drop=F]
    weights_core <- unique(gnames_core$weight)

    outlist_core <- list()

    for (weight in weights_core){

      cat("weight", weight, "\n")
      weight_name <- tools::file_path_sans_ext(basename(weight))
      gnames_core_weight <- gnames_core$gname[gnames_core$weight==weight]

      db <- RSQLite::dbConnect(sqlite, weight)
      query <- function(...) RSQLite::dbGetQuery(db, ...)

      for (gname in gnames_core_weight) {

        # if more than one weight file, append gene name with weight name
        if (length(weights)>1){
          gname_weight <- paste0(gname, "|", weight_name)
        } else {
          gname_weight <- gname
        }

        wgt <- query("select * from weights where gene = ?", params = list(gname))
        wgt.matrix <- as.matrix(wgt[, "weight", drop = F])
        rownames(wgt.matrix) <- wgt$rsid

        chrpos <- do.call(rbind, strsplit(wgt$varID, "_"))

        chrom <- as.integer(gsub("chr", "", chrpos[1, 1]))

        snps <- data.frame(chrom = as.integer(chrom),
                           id = wgt$rsid,
                           cm = "0",
                           pos = as.integer(chrpos[, 2]),
                           alt = wgt$eff_allele,
                           ref = wgt$ref_allele,
                           stringsAsFactors = F)

        if (!chrom %in% chr)
          next

        # select weight matrix by the prediction method used for this gene
        g.method <- "weight"
        wgt.matrix <- wgt.matrix[abs(wgt.matrix[, g.method]) > 0, , drop = F]
        wgt.matrix <- wgt.matrix[complete.cases(wgt.matrix),, drop = F]
        if (nrow(wgt.matrix) == 0)
          next

        if (is.null(z_snp)) {
          # select SNPs in both weight and LD reference
          snpnames <- intersect(rownames(wgt.matrix), ld_snpinfo$id)
        } else {
          # if z_snp is available, select SNPs in weight and LD reference and SNP zscores
          snpnames <- Reduce(intersect, list(rownames(wgt.matrix), ld_snpinfo$id, z_snp$id))
        }

        if (length(snpnames) == 0)
          next

        wgt.idx <- match(snpnames, rownames(wgt.matrix))
        wgt <- wgt.matrix[wgt.idx, g.method, drop = F]

        # scale weights by standard deviation of variant in LD reference
        if (isTRUE(scale_by_ld_variance)){
          ld_snpinfo.idx <- match(snpnames, ld_snpinfo$id)
          wgt <- wgt*sqrt(ld_snpinfo$variance[ld_snpinfo.idx])
        }

        p0 <- min(snps[snps[, "id"] %in% snpnames, "pos"])
        p1 <- max(snps[snps[, "id"] %in% snpnames, "pos"])
        nwgt <- nrow(wgt.matrix)
        nmiss <- nrow(wgt.matrix) - length(snpnames)
        # TODO: check missrate calculation
        outlist_core[[gname_weight]] <- list(chrom = chrom, p0 = p0, p1 = p1,
                                             wgt = wgt, gname=gname, weight_name=weight_name,
                                             n = nwgt, nmiss = nmiss, missrate = nwgt/nmiss)
      }
      RSQLite::dbDisconnect(db)
    }
    outlist_core
  }

  parallel::stopCluster(cl)

  exprlist <- lapply(names(outlist), function(x){
    outlist[[x]][c("chrom","p0","p1","wgt","gname","weight_name")]})
  names(exprlist) <- names(outlist)

  qclist <- lapply(names(outlist), function(x){
    outlist[[x]][c("n","nmiss","missrate")]})
  names(qclist) <- names(outlist)

  wgtlist <- lapply(outlist, "[[", "wgt")
  names(wgtlist) <- names(outlist)

  weight_info <- lapply(names(outlist), function(x){
    as.data.frame(outlist[[x]][c("chrom", "p0","p1", "gname", "weight_name", "n", "nmiss", "missrate")])})
  weight_info <- do.call(rbind, weight_info)
  rownames(weight_info) <- names(outlist)

  return(list(exprlist = exprlist, qclist = qclist,
              wgtlist = wgtlist, weight_info = weight_info))
}

#' Get region info with filenames of LD matrices and SNP information
#'
#' @param regions A data frame of LD regions
#'
#' @param LD_R_dir Directory of LD reference files
#'
#' @param pattern pattern of the LD reference filenames
#'
#' @param prefix prefix of the LD reference filenames
#'
#' @param LD_matrix_ext File extension of LD matrix files
#'
#' @param snp_info_ext File extension of SNP information files
#'
#' @return A data frame with information of the variants in the LD matrix.
#' @export
get_UKB_LD_region_info <- function(regions,
                                   LD_R_dir,
                                   pattern = "%s_chr%d.R_snp.%d_%d",
                                   prefix = "ukb_b37_0.1",
                                   LD_matrix_ext = "RDS",
                                   snp_info_ext = "Rvar") {

  LD_file <- sprintf(pattern, prefix, regions$chr, regions$start, regions$stop)

  if (is.null(regions$region_tag)){
    regions$region_tag <- paste0(regions$chr, ":", regions$start, "-", regions$stop)
  }

  region_info <- data.frame(regions,
                            LD_matrix = file.path(LD_R_dir, paste0(LD_file, ".", LD_matrix_ext)),
                            SNP_info = file.path(LD_R_dir, paste0(LD_file, ".", snp_info_ext)))

  return(region_info)
}
