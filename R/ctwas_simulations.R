#' Impute expression
#' 
#' @param pgenfs A character vector of .pgen or .bed files. One file for one
#' chromosome, in the order of 1 to 22. Therefore, the length of this vector
#' needs to be 22. If .pgen files are given, then .pvar and .psam are assumed
#' to present in the same directory. If .bed files are given, then .bim and
#' .fam files are assumed to present in the same directory.
#' 
#' @param weight a string, pointing to a directory with the fusion/twas format of weights, 
#' or a vector of one or more .db files in PredictDB format.
#' 
#' @param method a string, blup/bslmm/lasso/top1/enet/best. This option is only used for fusion weights. 
#' "best" means the method giving the best cross #' validation R^2. Note that top1 uses only the weight 
#' with largest effect. 
#'   
#' @param outputdir a string, the directory to store output
#' 
#' @param outname a string, the output name
#' 
#' @param logfile the log file, if NULL will print log info on screen
#' 
#' @param compress TRUE/FALSE. If TRUE, the imputed expression files are compressed
#' 
#' @param harmonize_wgt TRUE/FALSE. If TRUE, GWAS and eQTL genotype alleles are harmonized
#' 
#' @param strand_ambig_action_wgt the action to take to harmonize strand ambiguous variants (A/T, G/C) between 
#' the weights and LD reference. "drop" removes the ambiguous variant from the prediction models. "none" treats the variant 
#' as unambiguous, flipping the weights to match the LD reference and then taking no additional action. "recover" uses a procedure
#' to recover strand ambiguous variants. This procedure compares correlations between variants in the 
#' LD reference and prediction models, and it can only be used with PredictDB format prediction models, which include this
#' information.
#' 
#' @param ncore The number of cores used to parallelize imputation over weights
#' 
#' @param scale_by_ld_variance TRUE/FALSE. If TRUE, PredictDB weights are scaled by genotype variance, which is the default 
#' behavior for PredictDB
#' 
#' @importFrom logging addHandler loginfo
#'
#' @export
impute_expr <- function(pgenfs,
                        weight,
                        method = "lasso",
                        outputdir = getwd(),
                        outname = NULL,
                        logfile = NULL,
                        compress = T,
                        harmonize_wgt = T,
                        ncore=1,
                        scale_by_ld_variance=F,
                        strand_ambig_action_wgt = c("drop", "none", "recover")){

  if (!is.null(logfile)){
    addHandler(writeToFile, file= logfile, level='DEBUG')
  }

  if (isTRUE(harmonize_wgt)){
    loginfo("Harmonize set to True: will flip weights to match LD reference")
  }

  outname <- file.path(outputdir, outname)

  exprfs <- vector()

  for (b in 1:22){
    pgenf <- pgenfs[b]
    pvarf <- prep_pvar(pgenf, outputdir = outputdir)
    snpinfo <- read_pvar(pvarf)

    pgen <- prep_pgen(pgenf, pvarf)

    b <- unique(snpinfo$chrom)
    if (length(b) !=1){
      stop("Input genotype not splited by chromosome")
    }

    exprf <- paste0(outname, "_chr", b, ".expr")
    exprvarf <- paste0(outname, "_chr", b, ".exprvar")
    exprqcf <-  paste0(outname, "_chr", b, ".exprqc.Rd")

    loginfo("Reading weights for chromosome %s", b)
    
    if (isTRUE(dir.exists(weight))) {
      weight_type <- "fusion"
      weightall <- read_weight_fusion(weight, 
                                      b, 
                                      snpinfo, 
                                      z_snp=NULL, 
                                      method = method, 
                                      harmonize_wgt=harmonize_wgt)
    } else if (all(file_ext(weight) == "db")) {
      weight_type <- "predictdb"
      weightall <- read_weight_predictdb(weight, 
                                         b, 
                                         snpinfo, 
                                         z_snp=NULL, 
                                         harmonize_wgt=harmonize_wgt, 
                                         strand_ambig_action=strand_ambig_action_wgt,
                                         ncore=ncore,
                                         scale_by_ld_variance=F)
    } else {
      stop("Unrecognized weight format, need to use either FUSION format or predict.db format")
    }

    exprlist <- weightall[["exprlist"]]
    qclist <- weightall[["qclist"]]
    
    if (length(exprlist) > 0) {
      loginfo("Start gene experssion imputation ...")
      for (gname in names(exprlist)){
        wgt <- exprlist[[gname]][["wgt"]]
        snpnames <- rownames(wgt)
        gwas.idx <-  match(snpnames, snpinfo$id)
        g <- read_pgen(pgen, variantidx = gwas.idx)
        
        if (weight_type=="fusion" | scale_by_ld_variance==F){
          g <- scale(g)  # genotypes are standardized for FUSION
        }
        
        gexpr <- g %*% wgt
        exprlist[[gname]][["expr"]] <- gexpr
      }
    }
    
    loginfo ("Imputation done, writing results to output...")
    expr <- do.call(cbind, lapply(exprlist, '[[', "expr"))
    gnames <- names(exprlist)
    chrom <- unlist(lapply(exprlist,'[[', "chrom"))
    p0 <- unlist(lapply(exprlist,'[[', "p0"))
    p1 <- unlist(lapply(exprlist,'[[', "p1"))
    wgtlist <- lapply(exprlist, '[[', "wgt")

    if (length(exprlist) == 0){
      expr <- data.table::data.table(NULL)
      geneinfo <- data.table::data.table(NULL)
    }

    loginfo("Number of genes with imputed expression: %s for chr %s",
            ncol(expr), b)

    data.table::fwrite(expr, file = exprf,
                       row.names = F, col.names = F,
                       sep = "\t", quote = F)

    if (isTRUE(compress)){
      system(paste0("gzip -f ", exprf))
      exprf <- paste0(exprf, '.gz')
    }

    geneinfo <- data.frame("chrom" = chrom,
                           "id" = gnames,
                           "p0" = p0,
                           "p1" = p1)
    data.table::fwrite(geneinfo, file = exprvarf, sep = "\t", quote = F)

    save(wgtlist, qclist, file = exprqcf)

    loginfo('expression inmputation done for chr %s.', b)

    exprfs[b] <- exprf
  }

  return(exprfs)
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
          w <- harmonize_wgt_ld_old(wgt.matrix, 
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


harmonize_wgt_ld_old <- function (wgt.matrix, snps, ld_snpinfo, strand_ambig_action = c("drop", "none", "recover"), 
                              ld_pgenfs=NULL, ld_Rinfo=NULL, R_wgt=NULL, wgt=NULL){
  
  strand_ambig_action <- match.arg(strand_ambig_action)
  
  colnames(snps) <- c("chrom", "id", "cm", "pos", "alt", "ref")
  snps <- snps[match(rownames(wgt.matrix), snps$id), ]
  snpnames <- intersect(snps$id, ld_snpinfo$id)
  
  if (length(snpnames) != 0) {
    snps.idx <- match(snpnames, snps$id)
    ld.idx <- match(snpnames, ld_snpinfo$id)
    qc <- allele.qc(snps[snps.idx,]$alt, snps[snps.idx,]$ref, ld_snpinfo[ld.idx,]$alt, ld_snpinfo[ld.idx,]$ref)
    ifflip <- qc[["flip"]]
    ifremove <- !qc[["keep"]]
    flip.idx <- snps.idx[ifflip]
    snps[flip.idx, c("alt", "ref")] <- snps[flip.idx, c("ref", "alt")]
    wgt.matrix[flip.idx, ] <- -wgt.matrix[flip.idx, ]
    
    if (strand_ambig_action == "recover" & 
        any(ifremove) & 
        sum(!ifremove)>0){
      if (is.null(ld_pgenfs)){
        #load correlation matrix(es) for LD reference(s) containing current weight
        wgt_pos <- ld_snpinfo$pos[ld_snpinfo$id %in% snpnames]
        regnames <- unique(sapply(wgt_pos, function(x){which(x >= ld_Rinfo$start & x <= ld_Rinfo$stop)}))
        regRDS <- ld_Rinfo$RDS_file[match(regnames, ld_Rinfo$region_name)]
        R_snp <- lapply(regRDS, readRDS)
        R_snp <- suppressWarnings({as.matrix(Matrix::bdiag(R_snp))})
        R_snp_anno <- do.call(rbind, lapply(regRDS, read_ld_Rvar_RDS))
        
        #index the variant positions in LD reference
        R_snp.idx <- match(snpnames, R_snp_anno$id)
        R_snp.idx.unambig <- R_snp.idx[!ifremove]
        
        #drop R_wgt correlations between ambiguous variants
        R_wgt <- R_wgt[R_wgt$RSID1 %in% wgt$varID[!ifremove] | R_wgt$RSID2 %in% wgt$varID[!ifremove],]
        
        #flip correlations if weights were flipped in previous step
        for (i in flip.idx){
          ifflip_rwgt <- R_wgt$RSID1 == wgt$varID[i] | R_wgt$RSID2 == wgt$varID[i]
          R_wgt$VALUE[ifflip_rwgt] <- -R_wgt$VALUE[ifflip_rwgt]
        }
        
        unrecoverable.idx <- c()
        
        #iterate over ambiguous snps
        for (i in snpnames[ifremove]){
          #index current ambiguous snp
          snpnames.idx <- match(i, snpnames)
          
          #sum of correlations in LD reference between current ambiguous variant and unambiguous variants
          sumcor_R_snp <- sum(R_snp[R_snp.idx[snpnames.idx], R_snp.idx.unambig])
          
          #sum of correlations in weights between current ambiguous variant and unambiguous variants
          sumcor_R_wgt <- sum(R_wgt$VALUE[R_wgt$RSID1==wgt$varID[snpnames.idx] | R_wgt$RSID2==wgt$varID[snpnames.idx]])
          
          if (sumcor_R_snp==0 | sumcor_R_wgt==0){
            #collect ambiguous variants that do not have an unambiguous variant in the same LD region: all off-diagonal correlations = 0
            #also collect ambiguous variants independent of unambiguous variants in weights (trivial, correlations must = exactly zero)
            unrecoverable.idx <- c(unrecoverable.idx, snpnames.idx)
          } else {
            #flip weight if sign of correlations is not the same
            if (sign(sumcor_R_snp)!=sign(sumcor_R_wgt)){
              wgt.matrix[snpnames.idx,] <- -wgt.matrix[snpnames.idx,]
            }
          }
        }
        
        #drop ambiguous variants that cannot be recovered
        if (length(unrecoverable.idx)>0){
          snps <- snps[-unrecoverable.idx, , drop = F]
          wgt.matrix <- wgt.matrix[-unrecoverable.idx, , drop = F]
        }
      } else {
        stop("Recovering strand ambiguous variants is not currently suppported when using .pgen files")
        #TO-DO: mirror following section but compute R_snp for each region using individual level data
      }
    } else if (strand_ambig_action == "recover" & 
               any(ifremove) & 
               sum(ifremove)==1){
      #take no action if single variant. wrote this as separate if-statement for clarity, but it could be rolled into the following if-statement
    } else if (strand_ambig_action != "none" & 
               any(ifremove)){
      #if dropping ambiguous variants, or >2 ambiguous variants and 0 unambiguous variants, discard the ambiguous variants
      remove.idx <- snps.idx[ifremove]
      snps <- snps[-remove.idx, , drop = F]
      wgt.matrix <- wgt.matrix[-remove.idx, , drop = F]
    }
  }
  return(list(wgt = wgt.matrix, snps = snps))
}

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
