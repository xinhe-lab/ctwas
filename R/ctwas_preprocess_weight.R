
#' Preprocess PredictDB weights and harmonize with LD reference
#' (adapted from preharmonize_wgt_ld)
#'
#' @param weight a string, pointing to a directory with the fusion/twas format of weights, or a .db file in predictdb format.
#' A vector of multiple sets of weights in PredictDB format can also be specified; genes will have their filename appended
#' to their gene name to ensure IDs are unique.
#'
#' @param region_info a data frame of region definition and associated file names.
#'
#' @param drop_strand_ambig TRUE/FALSE, if TRUE remove strand ambiguous variants (A/T, G/C).
#'
#' @param write_db TRUE/FALSE, if TRUE, write processed weights as .db file in predictdb format
#'
#' @param outputdir a string, the directory to store output
#'
#' @param outname a string, the output name.
#'
#' @param logfile the log file, if NULL will print log info on screen.
#'
#' @importFrom logging addHandler loginfo writeToFile
#'
#' @return a list of processed weight table and extra table
#'
#' @export
#'
preprocess_weight <- function(weight,
                              region_info,
                              drop_strand_ambig = TRUE,
                              write_db = FALSE,
                              outputdir = getwd(),
                              outname = NULL,
                              logfile = NULL){

  # read the PredictDB weights
  sqlite <- RSQLite::dbDriver("SQLite")
  db = RSQLite::dbConnect(sqlite, weight)
  query <- function(...) RSQLite::dbGetQuery(db, ...)
  weight_table <- query("select * from weights")
  extra_table <- query("select * from extra")
  RSQLite::dbDisconnect(db)

  # subset to protein coding genes only
  loginfo("Subset to protein coding genes only")
  extra_table <- extra_table[extra_table$gene_type=="protein_coding",,drop=F]
  weight_table <- weight_table[weight_table$gene %in% extra_table$gene,]

  # read and subset the covariances
  weight_info <- read.table(gzfile(paste0(tools::file_path_sans_ext(weight), ".txt.gz")), header = T)
  weight_info <- weight_info[weight_info$GENE %in% extra_table$gene,]

  gnames <- unique(weight_table$gene)
  loginfo("Number of genes with weights provided: %d", length(gnames))

  # load all variants in LD reference
  ld_snpinfo <- do.call(rbind, lapply(region_info$SNP_info, read_LD_SNP_file))

  # remove variants in weight table, but not in LD reference
  loginfo("Number of variants in weights: %s", length(unique(weight_table$rsid)))
  loginfo("Remove %d variants in weights but not in LD reference", length(setdiff(weight_table$rsid, ld_snpinfo$id)))
  weight_table <- weight_table[weight_table$rsid %in% ld_snpinfo$id, ]

  # remove genes with no variants in LD reference
  loginfo("Remove %d genes with no variants in LD reference", length(setdiff(gnames, weight_table$gene)))
  gnames <- unique(weight_table$gene)
  loginfo("Number of genes left after filtering by LD reference: %d", length(gnames))

  # subset to variants in weight table
  ld_snpinfo <- ld_snpinfo[ld_snpinfo$id %in% weight_table$rsid,]

  loginfo("Harmonizing weights for %d genes...", length(gnames))

  weight_table_harmonized <- list()

  for (i in 1:length(gnames)){

    if (i %% 1000 == 0){
      loginfo("Current gene: %s", i)
    }

    gname <- gnames[i]
    wgt <- weight_table[weight_table$gene==gname,]
    wgt.matrix <- as.matrix(wgt[, "weight", drop = F])
    rsid_varID <- wgt[,c("rsid", "varID")]
    rownames(wgt.matrix) <- wgt$rsid
    chrpos <- do.call(rbind, strsplit(wgt$varID, "_"))

    snps <- data.frame(gsub("chr", "", chrpos[, 1]), wgt$rsid,
                       "0", chrpos[, 2], wgt$eff_allele, wgt$ref_allele,
                       stringsAsFactors = F)
    colnames(snps) <- c("chrom", "id", "cm", "pos", "alt", "ref")
    snps$chrom <- as.integer(snps$chrom)
    snps$pos <- as.integer(snps$pos)

    w <- harmonize_wgt_ld(wgt.matrix, snps, ld_snpinfo, drop_strand_ambig)
    wgt.matrix <- w[["wgt"]]
    snps <- w[["snps"]]

    wgt.matrix <- wgt.matrix[abs(wgt.matrix[, "weight"]) > 0, , drop = F]
    wgt.matrix <- wgt.matrix[complete.cases(wgt.matrix),, drop = F]

    snpnames <- intersect(rownames(wgt.matrix), ld_snpinfo$id)
    wgt.idx <- match(snpnames, rownames(wgt.matrix))
    wgt <- wgt.matrix[wgt.idx, "weight", drop = F]
    snps.idx <- match(snpnames, snps$id)
    snps <- snps[snps.idx,]

    if (length(snpnames)>0){
      weight_table_harmonized[[gname]] <- data.frame(gene=gname,
                                                     rsid=snps$id,
                                                     varID=rsid_varID$varID[match(snps$id, rsid_varID$rsid)],
                                                     ref_allele=snps$ref,
                                                     eff_allele=snps$alt,
                                                     weight=wgt[,"weight"])
    }
  }
  weight_table_harmonized <- do.call(rbind, weight_table_harmonized)
  loginfo("Number of genes with weights harmonized with LD reference: %s", length(unique(weight_table_harmonized$gene)))

  extra_table <- extra_table[extra_table$gene %in% weight_table_harmonized$gene,]

  if (isTRUE(write_db)) {
    if (!dir.exists(outputdir)){
      dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)
    }

    if (file.exists(file.path(outputdir, paste0(outname, ".db")))){
      invisible(file.remove(file.path(outputdir, paste0(outname, ".db"))))
      invisible(file.remove(file.path(outputdir, paste0(outname, ".txt.gz"))))
    }

    loginfo("Write preprocessed weights to %s", file.path(outputdir, paste0(outname, ".db")))

    db <- RSQLite::dbConnect(sqlite, file.path(outputdir, paste0(outname, ".db")))
    RSQLite::dbWriteTable(db, "extra", extra_table)
    RSQLite::dbWriteTable(db, "weights", weight_table_harmonized)
    RSQLite::dbDisconnect(db)

    weight_info_gz <- gzfile(file.path(outputdir, paste0(outname, ".txt.gz")), "w")
    write.table(weight_info, weight_info_gz, sep=" ", quote=F, row.names=F, col.names=T)
    close(weight_info_gz)
  }

  return(list(weight_table = weight_table_harmonized, extra_table = extra_table))
}

#' read (preprocessed) weights in PredictDB format
read_weights <- function (weight_files,
                          chrom,
                          ld_snpinfo,
                          z_snp = NULL,
                          scale_by_ld_variance=TRUE,
                          ncore=1){

  sqlite <- RSQLite::dbDriver("SQLite")

  # browser()
  # read gene names in each weight file
  gnames_all <- list()
  for (i in 1:length(weight_files)){
    weight_file <- weight_files[i]
    stopifnot(file.exists(weight_file))
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
          g_wgt_id <- paste0(gname, "|", weight_name)
        } else {
          g_wgt_id <- gname
        }

        wgt <- query("select * from weights where gene = ?", params = list(gname))
        wgt.matrix <- as.matrix(wgt[, "weight", drop = F])
        rownames(wgt.matrix) <- wgt$rsid

        chrpos <- do.call(rbind, strsplit(wgt$varID, "_"))

        snps <- data.frame(chrom = as.integer(gsub("chr", "", chrpos[, 1])),
                           id = wgt$rsid,
                           cm = "0",
                           pos = as.integer(chrpos[, 2]),
                           alt = wgt$eff_allele,
                           ref = wgt$ref_allele,
                           stringsAsFactors = F)

        if (!any(snps$chrom==chrom)){
          next
        }

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
        outlist_core[[g_wgt_id]] <- list(chrom = chrom, p0 = p0, p1 = p1,
                                         wgt = wgt, gene_name=gname, weight_name=weight_name,
                                         n = nwgt, nmiss = nmiss, missrate = nwgt/nmiss)
      }
      RSQLite::dbDisconnect(db)
    }
    outlist_core
  }

  parallel::stopCluster(cl)

  # exprlist <- lapply(names(outlist), function(x){
  #   outlist[[x]][c("chrom","p0","p1","wgt","gname","weight_name")]})
  # names(exprlist) <- names(outlist)
  #
  # qclist <- lapply(names(outlist), function(x){
  #   outlist[[x]][c("n","nmiss","missrate")]})
  # names(qclist) <- names(outlist)

  weight_list <- lapply(outlist, "[[", "wgt")
  names(weight_list) <- names(outlist)

  weight_info <- lapply(names(outlist), function(x){
    as.data.frame(outlist[[x]][c("chrom", "p0","p1", "gene_name", "weight_name", "n", "nmiss", "missrate")])})
  weight_info <- do.call(rbind, weight_info)
  weight_info$id <- names(outlist)
  rownames(weight_info) <- names(outlist)

  return(list(weight_list = weight_list, weight_info = weight_info))
}
