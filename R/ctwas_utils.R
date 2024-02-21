
#' read LD Rvar file as a data frame
read_ld_Rvar_file <- function(ld_Rvar_file){
  ld_Rvar <- data.table::fread(ld_Rvar_file, header = T)
  header <- c("chrom", "id", "pos", "alt", "ref")
  if (!all(header %in% colnames(ld_Rvar))){
    stop("The .Rvar file needs to contain the following columns: ",
         paste(header, collapse = " "), "\nCheck: ", ld_Rvar_file)
  }
  if (length(unique(ld_Rvar$chrom)) != 1){
    stop("LD region needs to be on only one chromosome. \nCheck: ", ld_Rvar_file)
  }

  return(ld_Rvar)
}

#' read LD Rvar files from all regions and combine SNP info a data frame
#' adapted from read_ld_Rvar()
read_ld_Rvar_snp_info <- function(Rvar_files){
  if (!is.null(Rvar_files)){
    ld_Rvar_snpinfo <- do.call(rbind, lapply(Rvar_files, read_ld_Rvar_file))
  } else {
    ld_Rvar_snpinfo <- data.table::data.table(chrom=as.integer(), id=as.character(), pos=as.integer(), alt=as.character(), ref=as.character(), variance=as.numeric())
  }
  return(ld_Rvar_snpinfo)
}


#' read and check Rvar from the region_info table, adapted from write_ld_Rf()
read_region_ld_Rinfo <- function(region_info){

  region_info$region_start <- region_info$start
  region_info$region_stop <- region_info$stop

  # read and check Rvar data and extract Rvar start and stop positions in each region
  for (i in 1:nrow(region_info)){
    Rvar <- read_ld_Rvar_file(region_info$Rvar_file[i])
    region_info$chrom[i] <- unique(Rvar$chrom)
    region_info$start[i] <- min(Rvar$pos)
    region_info$stop[i] <- max(Rvar$pos) + 1
  }

  region_info <- data.frame(region_info, stringsAsFactors = F)
  region_info <- transform(region_info,
                           chrom = as.numeric(chrom),
                           start = as.numeric(start),
                           stop = as.numeric(stop))

  chroms <- sort(unique(region_info$chrom))
  loginfo("chromosomes found in region_info: %s", chroms)

  # save a data frame for each chromosome with region_name and region_id
  ld_Rinfo_list <- vector("list", length = 22)

  for (b in chroms) {
    ld_Rinfo_chr <- region_info[region_info$chrom == b, , drop = F]
    ld_Rinfo_chr <- ld_Rinfo_chr[order(ld_Rinfo_chr$start), ]
    ld_Rinfo_chr$region_name <- 1:nrow(ld_Rinfo_chr)
    ld_Rinfo_chr$region_id <- paste0("chr", ld_Rinfo_chr$chrom, ":", ld_Rinfo_chr$start, "-",ld_Rinfo_chr$stop)

    ld_Rinfo_list[[b]] <- ld_Rinfo_chr
  }

  return(ld_Rinfo_list)
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

read_weight_fusion <- function (weight, chrom,
                                ld_snpinfo,
                                z_snp = NULL,
                                method = "lasso",
                                harmonize_wgt = T,
                                strand_ambig_action=c("drop", "none")){

  strand_ambig_action <- match.arg(strand_ambig_action)

  weight_name <- tools::file_path_sans_ext(basename(weight))

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
        w <- harmonize_wgt_ld(wgt.matrix, snps, ld_snpinfo, strand_ambig_action=strand_ambig_action)
        wgt.matrix <- w[["wgt"]]
        snps <- w[["snps"]]
      } else {
        colnames(snps) <- c("chrom", "id", "cm", "pos", "alt", "ref")
      }
      g.method = method
      if (g.method == "best") {
        g.method = names(which.max(cv.performance["rsq",]))
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
      p0 <- min(snps$pos[snps$id %in% snpnames])
      p1 <- max(snps$pos[snps$id %in% snpnames])
      exprlist[[gname]] <- list(chrom = chrom, p0 = p0, p1 = p1, wgt = wgt, gname=gname, weight_name=weight_name)
      nwgt <- nrow(wgt.matrix)
      nmiss <- nrow(wgt.matrix) - length(snpnames)
      qclist[[gname]] <- list(n = nwgt, nmiss = nmiss, missrate = nwgt/nmiss)
    }
  }
  return(list(exprlist = exprlist, qclist = qclist))
}

read_weight_predictdb <- function (weight,
                                   chrom,
                                   ld_snpinfo,
                                   z_snp = NULL,
                                   harmonize_wgt = T,
                                   strand_ambig_action = c("drop", "none", "recover"),
                                   ld_pgenfs=NULL,
                                   ld_Rinfo=NULL,
                                   scale_by_ld_variance=T,
                                   ncore=1){

  strand_ambig_action <- match.arg(strand_ambig_action)

  exprlist <- list()
  qclist <- list()
  weights <- weight

  sqlite <- RSQLite::dbDriver("SQLite")

  gnames_all <- list()

  for (i in 1:length(weights)){
    weight <- weights[i]

    db = RSQLite::dbConnect(sqlite, weight)
    query <- function(...) RSQLite::dbGetQuery(db, ...)
    gnames <- unique(query("select gene from weights")[, 1])

    gnames_all[[i]] <- cbind(gnames,weight)

    RSQLite::dbDisconnect(db)
  }

  gnames_all <- as.data.frame(do.call(rbind, gnames_all))
  colnames(gnames_all) <- c("gname", "weight")

  loginfo("Number of genes with weights provided: %s", nrow(gnames_all))
  loginfo("Collecting gene weight information ...")

  if (harmonize_wgt){
    loginfo("Flipping weights to match LD reference")
    if (strand_ambig_action=="recover"){
      loginfo("Harmonizing strand ambiguous weights using correlations with unambiguous variants")
    }
  }

  corelist <- lapply(1:ncore, function(core){njobs <- ceiling(nrow(gnames_all)/ncore); jobs <- ((core-1)*njobs+1):(core*njobs); jobs[jobs<=nrow(gnames_all)]})
  names(corelist) <- 1:ncore

  cl <- parallel::makeCluster(ncore, outfile = "")
  doParallel::registerDoParallel(cl)

  outlist <- foreach(core = 1:ncore, .combine = "c", .packages = "ctwas") %dopar% {
    gnames_core <- gnames_all[corelist[[core]],,drop=F]
    weights_core <- unique(gnames_core$weight)

    outlist_core <- list()

    for (weight in weights_core){
      loginfo("Current weight: %s (core %s)", weight, core)

      weight_name <- tools::file_path_sans_ext(basename(weight))
      gnames_core_weight <- gnames_core$gname[gnames_core$weight==weight]

      if (harmonize_wgt & strand_ambig_action=="recover"){
        R_wgt_all <- read.table(gzfile(paste0(file_path_sans_ext(weight), ".txt.gz")), header=T) #load covariances for variants in each gene (accompanies .db file)
        R_wgt_all <- R_wgt_all[R_wgt_all$GENE %in% gnames_core_weight,]
      }

      db = RSQLite::dbConnect(sqlite, weight)
      query <- function(...) RSQLite::dbGetQuery(db, ...)

      for (gname in gnames_core_weight) {

        if (length(weights)>1){
          gname_weight <- paste0(gname, "|", weight_name)
        } else {
          gname_weight <- gname
        }

        wgt <- query("select * from weights where gene = ?", params = list(gname))
        wgt.matrix <- as.matrix(wgt[, "weight", drop = F])

        rownames(wgt.matrix) <- wgt$rsid
        chrpos <- do.call(rbind, strsplit(wgt$varID, "_"))


        snps <- data.frame(gsub("chr", "", chrpos[, 1]), wgt$rsid,
                           "0", chrpos[, 2], wgt$eff_allele, wgt$ref_allele,
                           stringsAsFactors = F)
        colnames(snps) <- c("chrom", "id", "cm", "pos", "alt", "ref")
        snps$chrom <- as.integer(snps$chrom)
        snps$pos <- as.integer(snps$pos)

        if (!any(snps$chrom==chrom)){
          next
        }

        if (isTRUE(harmonize_wgt)) {
          if (strand_ambig_action=="recover"){
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
          w <- harmonize_wgt_ld(wgt.matrix,
                                snps,
                                ld_snpinfo,
                                strand_ambig_action=strand_ambig_action,
                                ld_Rinfo=ld_Rinfo,
                                R_wgt=R_wgt,
                                wgt=wgt)
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
        if (scale_by_ld_variance){
          ld_snpinfo.idx <- match(snpnames, ld_snpinfo$id)
          wgt <- wgt*sqrt(ld_snpinfo$variance[ld_snpinfo.idx])
        }

        p0 <- min(snps[snps[, "id"] %in% snpnames, "pos"])
        p1 <- max(snps[snps[, "id"] %in% snpnames, "pos"])
        nwgt <- nrow(wgt.matrix)
        nmiss <- nrow(wgt.matrix) - length(snpnames)
        outlist_core[[gname_weight]] <- list(chrom = chrom, p0 = p0, p1 = p1, wgt = wgt, gname=gname, weight_name=weight_name,
                                             n = nwgt, nmiss = nmiss, missrate = nwgt/nmiss)
      }

      RSQLite::dbDisconnect(db)
    }

    outlist_core
  }

  parallel::stopCluster(cl)

  exprlist_weight <- lapply(names(outlist), function(x){outlist[[x]][c("chrom","p0","p1","wgt","gname","weight_name")]})
  names(exprlist_weight) <- names(outlist)

  qclist_weight <- lapply(names(outlist), function(x){outlist[[x]][c("n","nmiss","missrate")]})
  names(qclist_weight) <- names(outlist)

  exprlist <- c(exprlist, exprlist_weight)
  qclist <- c(qclist, qclist_weight)

  rm(outlist, exprlist_weight, qclist_weight)

  return(list(exprlist = exprlist, qclist = qclist))
}
