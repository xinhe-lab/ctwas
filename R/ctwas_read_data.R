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

  pvarfout
}

#' Read .pvar file into R
#' @param pvarf .pvar file or .bim file with have proper
#'  .pvar file format: https://www.cog-genomics.org/plink/2.0/formats#pvar
#'
#' @return A data.table. variant info
#'
read_pvar <- function(pvarf){

  pvardt <- data.table::fread(pvarf, skip = "#CHROM")
  pvardt <- dplyr::rename(pvardt, "chrom" = "#CHROM", "pos" = "POS",
                "alt" = "ALT", "ref" = "REF", "id" = "ID")
  pvardt <- pvardt[, c("chrom", "id", "pos", "alt", "ref")]
  pvardt
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

  pgen
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


#' Prepare .exprvar file
#'
#' @param exprf expression variable info files, the output of \code{impute_expr}
#'
#' @return corresponding exprvar file
#'
#' @importFrom tools file_ext file_path_sans_ext
#'
prep_exprvar <- function(exprf){
  if (file_ext(exprf) == "gz"){
    exprf <- file_path_sans_ext(exprf)
  }
  exprvarf <- paste0(exprf, "var")
  exprvarf
}

#' Read .exprvar file into R
#'
#' @param exprvarf expression variable info files, prepared by the \code{prep_exprvar} function
#'
#' @return A data.table. variant info
#'
read_exprvar <- function(exprvarf){

  exprvar <- try(data.table::fread(exprvarf, header = T))

  if (inherits(exprvar, "try-error")){
    exprvar <-  setNames(data.table(matrix(nrow = 0, ncol = 4)),
                         c("chrom", "id", "p0", "p1"))
  }
  exprvar
}

#' Read .expr file into R
#'
#' @param exprf expression variable info files, the output of \code{impute_expr}
#'
#' @param variantidx variant index. If NULL, all variants will be extracted.
#'
#' @return A matrix, columns are imputed expression for each gene, rows are
#'  for each sample.
#'
read_expr <- function(exprf, variantidx = NULL){
  if (!is.null(variantidx) & length(variantidx)==0){
    return(NULL)
  } else{
    return(as.matrix(data.table::fread(exprf, header = F,
                                       select = variantidx)))
  }
}


#' read variant information associated with a LD R matrix .RDS file.
#'
#' @param ld_RDSf files containing the variant information for the LD matrices
#'
#' @return a data frame with columns: "chrom", "id", "pos", "alt", "ref". "alt" is
#' the coded allele
#'
#' @importFrom tools file_ext file_path_sans_ext
#'
read_ld_Rvar_RDS <- function(ld_RDSf){
  ld_Rvarf <- paste0(file_path_sans_ext(ld_RDSf), ".Rvar")
  ld_Rvar <- data.table::fread(ld_Rvarf, header = T)
  target_header <- c("chrom", "id", "pos", "alt", "ref")
  if (all(target_header %in% colnames(ld_Rvar))){
      return(ld_Rvar)
  } else {
    stop("The .Rvar file needs to contain the following columns: ",
         paste(target_header, collapse = " "))
  }
}

#' combine variant information associated with a LD R matrix .RDS file.
#'
#' @param ld_R_dir The directory that contains all ld R matrices.
#' the ld R matrices should not have overlapping positions.
#'
#' @param outputdir a string, the directory to store output
#'
#' @param outname a string, the output name
#'
#' @return A vector of the `ld_Rf` file names. The function will write one `ld_Rf` file
#' for each chromosome, so the vector has length 22. The `ld_Rf` file has the following
#' columns: chr region_name start stop RDS_file.
#'
write_ld_Rf <- function(ld_R_dir, outname = outname , outputdir = getwd()){
  ld_RDSfs <- list.files(path = ld_R_dir, pattern = "\\.RDS$", full.names = T)
  ldinfolist <- list()
  for (ld_RDSf in ld_RDSfs){
    Rvar <- read_ld_Rvar_RDS(ld_RDSf)
    chrom <- unique(Rvar$chrom)
    if (length(chrom) != 1){
      stop("R matrix on multiple chromosomes,
           can't handle this. Need to be on one chromosome:", ld_RDSf)
    }
    start <- min(Rvar$pos)
    stop <- max(Rvar$pos) + 1
    ldinfolist[[ld_RDSf]] <- c(chrom, start, stop, ld_RDSf)
  }
  ldinfo <- do.call(rbind, ldinfolist)
  colnames(ldinfo) <- c("chrom", "start", "stop", "RDS_file")
  rownames(ldinfo) <- NULL
  ldinfo <- data.frame(ldinfo, stringsAsFactors = F)
  ldinfo <- transform(ldinfo, chrom = as.numeric(chrom),
                      start = as.numeric(start),
                      stop = as.numeric(stop))

  ld_Rfs <- vector()
  for (b in 1:22) {
    ldinfo.b <- ldinfo[ldinfo$chrom == b, , drop = F]
    ldinfo.b <- ldinfo.b[order(ldinfo.b$start), ]

    if (nrow(ldinfo.b) == 0) {
      loginfo(paste0("no region on chromosome ", b))
      ldinfo.b <- cbind(ldinfo.b, data.frame(region_name=as.character()))
    } else {
      ldinfo.b$region_name <- 1:nrow(ldinfo.b)
    }

    ld_Rf <- file.path(outputdir, paste0(outname, "_ld_R_chr",
                                         b, ".txt"))
    write.table(ldinfo.b, file = ld_Rf, row.names = F, col.names = T,
                sep = "\t", quote = F)
    ld_Rfs[b] <- ld_Rf
  }
  ld_Rfs
}

#' read variant information for all ld matrices in `ld_Rf`.
#'
#' @param ld_Rf a vector of paths to the LD matrices
#'
#' @return a data frame with columns: "chrom", "id", "pos", "alt", "ref"
#'
read_ld_Rvar <- function(ld_Rf){
  Rinfo <- data.table::fread(ld_Rf, header = T)
  if (nrow(Rinfo)>0){
    ld_Rvar <- do.call(rbind, lapply(Rinfo$RDS_file, read_ld_Rvar_RDS))
  } else {
    ld_Rvar <- data.table::data.table(chrom=as.integer(), id=as.character(), pos=as.integer(), alt=as.character(), ref=as.character(), variance=as.numeric())
  }
  ld_Rvar
}


read_weight_fusion <- function (weight, chrom, ld_snpinfo, z_snp = NULL, method = "lasso", harmonize_wgt = T){
  exprlist <- list()
  qclist <- list()
  wgtdir <- dirname(weight)
  wgtposfile <- file.path(wgtdir, paste0(basename(weight),
                                         ".pos"))
  wgtpos <- read.table(wgtposfile, header = T, stringsAsFactors = F)
  wgtpos <- transform(wgtpos, ID = ifelse(duplicated(ID) | duplicated(ID, fromLast = TRUE),
                                          paste(ID, ave(ID, ID, FUN = seq_along), sep = "_ID"), ID))
  loginfo("number of genes with weights provided: %s", nrow(wgtpos))
  wgtpos <- wgtpos[wgtpos$CHR == chrom, ]
  loginfo("number of genes on chromosome %s: %s", chrom, nrow(wgtpos))
  loginfo("collecting gene weight information ...")
  if (nrow(wgtpos) > 0) {
    for (i in 1:nrow(wgtpos)) {
      wf <- file.path(wgtdir, wgtpos[i, "WGT"])
      load(wf)
      gname <- wgtpos[i, "ID"]
      if (isTRUE(harmonize_wgt)) {
        w <- harmonize_wgt_ld(wgt.matrix, snps, ld_snpinfo, recover_strand_ambig=F)
        wgt.matrix <- w[["wgt"]]
        snps <- w[["snps"]]
      }
      g.method = method
      if (g.method == "best") {
        g.method = names(which.max(cv.performance["rsq", ]))
      }
      if (exists("cv.performance")){
        if (!(g.method %in% names(cv.performance[1,]))){
          next
        }
      }

      # Ensure only top magnitude snp weight in the top1 wgt.matrix column
      if (g.method == "top1"){
        wgt.matrix[,"top1"][-which.max(wgt.matrix[,"top1"]^2)] <- 0
      }

      wgt.matrix <- wgt.matrix[abs(wgt.matrix[, g.method]) >
                                 0, , drop = F]
      wgt.matrix <- wgt.matrix[complete.cases(wgt.matrix),
                               , drop = F]
      if (nrow(wgt.matrix) == 0)
        next
      if (is.null(z_snp)) {
        snpnames <- intersect(rownames(wgt.matrix), ld_snpinfo$id)
      } else {
        snpnames <- Reduce(intersect, list(rownames(wgt.matrix),
                                           ld_snpinfo$id, z_snp$id))
      }
      if (length(snpnames) == 0)
        next
      wgt.idx <- match(snpnames, rownames(wgt.matrix))
      wgt <- wgt.matrix[wgt.idx, g.method, drop = F]
      p0 <- min(snps[snps[, "id"] %in% snpnames, "pos"])
      p1 <- max(snps[snps[, "id"] %in% snpnames, "pos"])
      exprlist[[gname]] <- list(chrom = chrom, p0 = p0,
                                p1 = p1, wgt = wgt)
      nwgt <- nrow(wgt.matrix)
      nmiss <- nrow(wgt.matrix) - length(snpnames)
      qclist[[gname]] <- list(n = nwgt, nmiss = nmiss,
                              missrate = nwgt/nmiss)
    }
  }
  return(list(exprlist = exprlist, qclist = qclist))
}

read_weight_predictdb <- function (weight, chrom, ld_snpinfo, z_snp = NULL, harmonize_wgt = T,
                                   recover_strand_ambig=T, ld_pgenfs=NULL, ld_Rinfo=NULL){
  exprlist <- list()
  qclist <- list()
  sqlite <- RSQLite::dbDriver("SQLite")
  db = RSQLite::dbConnect(sqlite, weight)
  query <- function(...) RSQLite::dbGetQuery(db, ...)
  gnames <- unique(query("select gene from weights")[, 1])
  loginfo("Number of genes with weights provided: %s", length(gnames))
  loginfo("Collecting gene weight information ...")
  if (harmonize_wgt){
    loginfo("Flipping weights to match LD reference")
    if (recover_strand_ambig){
      loginfo("Harmonizing strand ambiguous weights using correlations with unambiguous variants")
      R_wgt_all = read.table(gzfile(paste0(file_path_sans_ext(weight), ".txt.gz")), header=T) #load correlations for variants in each gene (accompanies .db file)
    }
  }
  for (gname in gnames) {
    wgt <- query("select * from weights where gene = ?",
                 params = list(gname))
    wgt.matrix <- as.matrix(wgt[, "weight", drop = F])
    rownames(wgt.matrix) <- wgt$rsid
    chrpos <- do.call(rbind, strsplit(wgt$varID, "_"))
    snps <- data.frame(gsub("chr", "", chrpos[, 1]), wgt$rsid,
                       "0", chrpos[, 2], wgt$eff_allele, wgt$ref_allele,
                       stringsAsFactors = F)
    colnames(snps) <- c("chrom", "id", "cm", "pos", "alt",
                        "ref")
    snps$chrom <- as.integer(snps$chrom)
    snps$pos <- as.integer(snps$pos)

    if (!any(snps$chrom==chrom)){
      next
    }

    if (isTRUE(harmonize_wgt)) {
      if (recover_strand_ambig){
        #subset R_wgt_all to current weight
        R_wgt <- R_wgt_all[R_wgt_all$GENE == gname,]

        #convert covariance to correlation
        R_wgt_stdev <- R_wgt[R_wgt$RSID1==R_wgt$RSID2,]
        R_wgt_stdev <- setNames(sqrt(R_wgt_stdev$VALUE), R_wgt_stdev$RSID1)
        R_wgt$VALUE <- R_wgt$VALUE/(R_wgt_stdev[R_wgt$RSID1]*R_wgt_stdev[R_wgt$RSID2])

        #discard variances
        R_wgt <- R_wgt[R_wgt$RSID1!=R_wgt$RSID2,]

        #fix edge case where variance=0; treat correlations with these variants as uninformative (=0) for harmonization
        R_wgt$VALUE[is.nan(R_wgt$VALUE)] <- 0
      } else {
        R_wgt <- NULL
      }
      w <- harmonize_wgt_ld(wgt.matrix, snps, ld_snpinfo,
                            recover_strand_ambig=recover_strand_ambig,
                            ld_Rinfo=ld_Rinfo, R_wgt=R_wgt, wgt=wgt, )
      wgt.matrix <- w[["wgt"]]
      snps <- w[["snps"]]
    }
    g.method = "weight"
    wgt.matrix <- wgt.matrix[abs(wgt.matrix[, g.method]) > 0, , drop = F]
    wgt.matrix <- wgt.matrix[complete.cases(wgt.matrix),, drop = F]
    if (nrow(wgt.matrix) == 0)
      next
    if (is.null(z_snp)) {
      snpnames <- intersect(rownames(wgt.matrix), ld_snpinfo$id)
    } else {
      snpnames <- Reduce(intersect, list(rownames(wgt.matrix), ld_snpinfo$id, z_snp$id))
    }
    if (length(snpnames) == 0)
      next
    wgt.idx <- match(snpnames, rownames(wgt.matrix))
    wgt <- wgt.matrix[wgt.idx, g.method, drop = F]

    #scale weights by standard deviation of variant in LD reference
    ld_snpinfo.idx <- match(snpnames, ld_snpinfo$id)
    wgt <- wgt*sqrt(ld_snpinfo$variance[ld_snpinfo.idx])

    p0 <- min(snps[snps[, "id"] %in% snpnames, "pos"])
    p1 <- max(snps[snps[, "id"] %in% snpnames, "pos"])
    exprlist[[gname]] <- list(chrom = chrom, p0 = p0, p1 = p1, wgt = wgt)
    nwgt <- nrow(wgt.matrix)
    nmiss <- nrow(wgt.matrix) - length(snpnames)
    qclist[[gname]] <- list(n = nwgt, nmiss = nmiss, missrate = nwgt/nmiss)
  }
  RSQLite::dbDisconnect(db)
  return(list(exprlist = exprlist, qclist = qclist))
}

