
# flip alleles
# Copied from allele.qc function
# in https://github.com/gusevlab/fusion_twas/blob/master/FUSION.assoc_test.R,
allele.qc = function(a1,a2,ref1,ref2) {
  a1 = toupper(a1)
  a2 = toupper(a2)
  ref1 = toupper(ref1)
  ref2 = toupper(ref2)

  ref = ref1
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip1 = flip

  ref = ref2
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip2 = flip

  snp = list()
  snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
  snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)

  return(snp)
}

#' Harmonize z scores from GWAS to match ld reference genotypes.
#' Flip signs when reverse complement matches.
#' 
#' @param z_snp a data frame, with columns "id", "A1", "A2" and "z".
#'     Z scores for every SNP. "A1" is the effect allele.
#'     
#' @param ld_snpinfo a data frame, snp info for LD reference,
#'  with columns "chrom", "id", "pos", "alt", "ref".
#'  
#' @param strand_ambig_action the action to take to harmonize strand ambiguous variants (A/T, G/C) between 
#' the z scores and LD reference. "drop" removes the ambiguous variant from the z scores. "none" treats the variant 
#' as unambiguous, flipping the z score to match the LD reference and then taking no additional action. "recover" 
#' imputes the sign of ambiguous z scores using unambiguous z scores and the LD reference and flips the z scores 
#' if there is a mismatch between the imputed sign and the observed sign of the z score. This option is computationally intensive
#' 
#' @param ld_pgenfs a character vector of .pgen or .bed files. One file for one
#'  chromosome, in the order of 1 to 22. Therefore, the length of this vector
#'  needs to be 22. If .pgen files are given, then .pvar and .psam are assumed
#'  to present in the same directory. If .bed files are given, then .bim and
#'  .fam files are assumed to present in the same directory.
#'  
#' @param ld_Rinfo a vector of paths to the variant information for all LD matrices 
#' 
#' @return a data frame, z_snp with the "z" columns flipped to match LD ref.
#' 
harmonize_z_ld <- function(z_snp, ld_snpinfo, strand_ambig_action = c("drop", "none", "recover"), ld_pgenfs = NULL, ld_Rinfo = NULL){
  strand_ambig_action <- match.arg(strand_ambig_action)
  snpnames <- intersect(z_snp$id, ld_snpinfo$id)
    if (length(snpnames) != 0) {
    z.idx <- match(snpnames, z_snp$id)
    ld.idx <- match(snpnames, ld_snpinfo$id)
    qc <- allele.qc(z_snp[z.idx,]$A1, z_snp[z.idx,]$A2, ld_snpinfo[ld.idx,]$alt, ld_snpinfo[ld.idx,]$ref)
    ifflip <- qc[["flip"]]
    ifremove <- !qc[["keep"]]
    flip.idx <- z.idx[ifflip]
    z_snp[flip.idx, c("A1", "A2")] <- z_snp[flip.idx, c("A2", "A1")]
    z_snp[flip.idx, "z"] <- -z_snp[flip.idx, "z"]
    if (strand_ambig_action=="recover" & any(ifremove)){
      #compare sign of imputed z score with observed z score for strand ambiguous variants 
      #following imputation strategy in https://dx.doi.org/10.1093%2Fbioinformatics%2Fbtu416
      if (is.null(ld_pgenfs)){
        loginfo("Harmonizing strand-ambiguous z scores using imputation by region")
        for (i in 1:nrow(ld_Rinfo)){
          R_snp <- readRDS(ld_Rinfo$RDS_file[i])
          R_snp_anno <- read_ld_Rvar_RDS(ld_Rinfo$RDS_file[i])
          
          #indicator if z.idx (on chromsome, in LD ref) is in current region
          z.idx.ifreg <- z_snp$id[z.idx] %in% R_snp_anno$id
          
          #index if z.idx (on chromsome, in LD ref) is in current region and unambiguous/ambiguous
          z.idx.unambig <- z.idx[z.idx.ifreg & qc$keep]
          z.idx.ambig <- z.idx[z.idx.ifreg & !qc$keep]
          
          #skip region if there are no ambiguous variants
          if (length(z.idx.ambig)>0){
            #z scores for unambiguous and ambiguous variants in current region
            z_t <- z_snp$z[z.idx.unambig]
            z_i_obs <- z_snp$z[z.idx.ambig]
            
            #impute z scores for ambiguous variants using unambiguous variants in current region
            R_t.idx <- match(z_snp$id[z.idx.unambig], R_snp_anno$id)
            R_i.idx <- match(z_snp$id[z.idx.ambig], R_snp_anno$id)
            lambda <- 0.001
            sigma_tt_inv <- solve(R_snp[R_t.idx, R_t.idx, drop=F] + lambda*diag(length(R_t.idx)))
            sigma_it <- R_snp[R_i.idx, R_t.idx, drop=F]
            z_i <- sigma_it%*%sigma_tt_inv%*%z_t
            
            #flip z scores that do not match the sign of the imputation
            #NOTE: consider replacing this with a test of significance - see "Fine-mapping from summary data with the Sum of Single Effects model" - Yuxin
            if_sign_neq <- sign(z_i_obs) != sign(z_i)
            z_snp[z.idx.ambig[if_sign_neq], "z"] <- -z_snp[z.idx.ambig[if_sign_neq], "z"]
          }
        }
      } else {
        #TO-DO: mirror previous section but compute R for each region using X
      }
    } else if (strand_ambig_action=="drop" & any(ifremove)) {
      remove.idx <- z.idx[ifremove]
      z_snp <- z_snp[-remove.idx, , drop = F]
    }
  }
  return(z_snp)
}

#' Harmonize z scores from GWAS to match ld reference genotypes.
#' Flip signs when reverse complement matches, remove/recover strand ambiguous SNPs
#' 
#' @param wgt.matrix from FUSION weight .Rdat file
#' 
#' @param snps from FUSION weight .Rdat file
#'  with columns "chrom", "id", "pos", "alt", "ref". The effect allele
#'  for FUSION is alt.
#'  
#' @param ld_snpinfo a data frame, snp info for LD reference,
#'  with columns "chrom", "id", "pos", "alt", "ref".
#' 
#' @param strand_ambig_action the action to take to harmonize strand ambiguous variants (A/T, G/C) between 
#' the weights and LD reference. "drop" removes the ambiguous variant from the prediction models. "none" treats the variant 
#' as unambiguous, flipping the weights to match the LD reference and then taking no additional action. "recover" uses a procedure
#' to recover strand ambiguous variants. This procedure compares correlations between variants in the 
#' LD reference and prediction models, and it can only be used with PredictDB format prediction models, which include this
#' information.
#'  
#' @param ld_pgenfs a character vector of .pgen or .bed files. One file for one
#'  chromosome, in the order of 1 to 22. Therefore, the length of this vector
#'  needs to be 22. If .pgen files are given, then .pvar and .psam are assumed
#'  to present in the same directory. If .bed files are given, then .bim and
#'  .fam files are assumed to present in the same directory.
#'  
#' @param ld_Rinfo a vector of paths to the variant information for all LD matrices 
#'  
#' @param R_wgt the LD matrix for the variants in \code{wgt.matrix} 
#' 
#' @param wgt allele information from the weights
#'  
#' @return wgt.matrix and snps with alleles flipped to match
#' 
harmonize_wgt_ld <- function (wgt.matrix, snps, ld_snpinfo, strand_ambig_action = c("drop", "none", "recover"), 
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

#' Harmonize GWAS summary statistics and LD reference
#' 
#' @param z_snp A data frame with two columns: "id", "A1", "A2", "z". giving the z scores for
#' snps. "A1" is effect allele. "A2" is the other allele. If `harmonize= False`, A1 and A2 are not required.
#' 
#' @param LD_R_dir a string, pointing to a directory containing all LD matrix files and variant information. Expects .RDS files which contain LD correlation matrices for a region/block.
#' For each RDS file, a file with same base name but ended with .Rvar needs to be present in the same folder. the .Rvar file has 5 required columns: "chrom", "id", "pos", "alt", "ref". 
#' If using PredictDB format weights and \code{scale_by_ld_variance=T}, a 6th column is also required: "variance", which is the variance of the each SNP.
#' The order of rows needs to match the order of rows in .RDS file.
#'   
#' @param outputdir a string, the directory to store output
#' 
#' @param outname a string, the output name
#' 
#' @param logfile the log file, if NULL will print log info on screen
#' 
#' @param harmonize_z TRUE/FALSE. If TRUE, GWAS and eQTL genotype alleles are harmonized
#' 
#' @param strand_ambig_action_z the action to take to harmonize strand ambiguous variants (A/T, G/C) between 
#' the z scores and LD reference. "drop" removes the ambiguous variant from the z scores. "none" treats the variant 
#' as unambiguous, flipping the z score to match the LD reference and then taking no additional action. "recover" 
#' imputes the sign of ambiguous z scores using unambiguous z scores and the LD reference and flips the z scores 
#' if there is a mismatch between the imputed sign and the observed sign of the z score. This option is computationally intensive
#' 
#' @param drop_multiallelic TRUE/FALSE. If TRUE, multiallelic variants will be dropped from the summary statistics
#' 
#' @importFrom logging addHandler loginfo
#' @importFrom tools file_ext
#'
#' @export
#' 
preharmonize_z_ld <- function (z_snp, ld_R_dir, outputdir = getwd(), outname = NULL, logfile = NULL,
                           harmonize_z = T, strand_ambig_action_z = c("drop", "none", "recover"),
                           drop_multiallelic=T){
  dir.create(outputdir, showWarnings = F, recursive=T)
  
  if (!is.null(logfile)) {
    addHandler(writeToFile, file = logfile, level = "DEBUG")
  }
  
  ld_Rfs <- write_ld_Rf(ld_R_dir, outname = outname, outputdir = outputdir)
  
  ld_snplist <- c()
  
  if (drop_multiallelic){
    z_snp <- z_snp[!(z_snp$id %in% z_snp$id[duplicated(z_snp$id)]),] #drop multiallelic variants (id not unique)
  }
  
  for (b in 1:22){
    loginfo("Harmonizing summary statistics for chromosome %s", b)
    
    ld_Rf <- ld_Rfs[b]
    ld_Rinfo <- data.table::fread(ld_Rf, header = T)
    ld_snpinfo <- read_ld_Rvar(ld_Rf)
    
    chrom <- unique(ld_snpinfo$chrom)
    if (length(chrom) > 1) {
      stop("Input LD reference not split by chromosome")
    }
    ld_snplist <- c(ld_snplist, ld_snpinfo$id) #store names of snps in ld reference
    
    if (isTRUE(harmonize_z)) {
      loginfo("Flipping z scores to match LD reference")
      z_snp <- harmonize_z_ld(z_snp, ld_snpinfo,
                              strand_ambig_action = strand_ambig_action_z, 
                              ld_pgenfs = NULL, 
                              ld_Rinfo = ld_Rinfo)
    }
  }
  
  z_snp <- z_snp[z_snp$id %in% ld_snplist,]

  return(list(z_snp = z_snp))
}

#' Harmonize PredictDB weights and LD reference
#' 
#' @param weight a string, pointing to a directory with the fusion/twas format of weights, or a .db file in predictdb format.
#' A vector of multiple sets of weights in PredictDB format can also be specified; genes will have their filename appended
#' to their gene name to ensure IDs are unique.
#' 
#' @param LD_R_dir a string, pointing to a directory containing all LD matrix files and variant information. Expects .RDS files which contain LD correlation matrices for a region/block.
#' For each RDS file, a file with same base name but ended with .Rvar needs to be present in the same folder. the .Rvar file has 5 required columns: "chrom", "id", "pos", "alt", "ref". 
#' If using PredictDB format weights and \code{scale_by_ld_variance=T}, a 6th column is also required: "variance", which is the variance of the each SNP.
#' The order of rows needs to match the order of rows in .RDS file.
#'   
#' @param outputdir a string, the directory to store output
#' 
#' @param outname a string, the output name
#' 
#' @param strand_ambig_action_wgt the action to take to harmonize strand ambiguous variants (A/T, G/C) between 
#' the weights and LD reference. "drop" removes the ambiguous variant from the prediction models. "none" treats the variant 
#' as unambiguous, flipping the weights to match the LD reference and then taking no additional action. "recover" uses a procedure
#' to recover strand ambiguous variants. This procedure compares correlations between variants in the 
#' LD reference and prediction models, and it can only be used with PredictDB format prediction models, which include this
#' information.
#' 
#' @importFrom logging addHandler loginfo
#' @importFrom tools file_ext
#'
#' @export
#' 
preharmonize_wgt_ld <- function (weight, 
                                 ld_R_dir, 
                                 outputdir = getwd(), 
                                 outname, 
                                 strand_ambig_action_wgt = c("drop", "none", "recover")){
  
  strand_ambig_action_wgt <- match.arg(strand_ambig_action_wgt)
  
  dir.create(outputdir, showWarnings = F)
  
  #read the PredictDB weights
  sqlite <- RSQLite::dbDriver("SQLite")
  db = RSQLite::dbConnect(sqlite, weight)
  query <- function(...) RSQLite::dbGetQuery(db, ...)
  weight_table <- query("select * from weights")
  extra_table <- query("select * from extra")
  RSQLite::dbDisconnect(db)
  
  gnames <- unique(weight_table$gene)
  loginfo("Number of genes with weights provided: %s", length(gnames))
  
  #load ld information and subset to variants in weights
  ld_Rfs <- ctwas:::write_ld_Rf(ld_R_dir, outname = outname, outputdir = outputdir)
  
  ld_Rinfo <- lapply(ld_Rfs, data.table::fread, header=T)
  ld_Rinfo <- as.data.frame(do.call(rbind, ld_Rinfo))
  
  ld_snpinfo <- lapply(ld_Rfs, ctwas:::read_ld_Rvar)
  ld_snpinfo <- as.data.frame(do.call(rbind, ld_snpinfo))
  ld_snpinfo <- ld_snpinfo[ld_snpinfo$id %in% weight_table$rsid,]
  
  loginfo("Flipping weights to match LD reference")
  
  if (strand_ambig_action_wgt=="recover"){
    R_wgt_all <- read.table(gzfile(paste0(file_path_sans_ext(weight), ".txt.gz")), header=T) #load covariances for variants in each gene (accompanies .db file)
    loginfo("Harmonizing strand ambiguous weights using correlations with unambiguous variants")
  }
  
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
    
    chrom <- unique(snps$chrom)
    ld_Rinfo_chrom <- ld_Rinfo[ld_Rinfo$chrom==chrom,]
    
    if (strand_ambig_action_wgt=="recover"){
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
    
    w <- ctwas:::harmonize_wgt_ld(wgt.matrix,
                                  snps,
                                  ld_snpinfo,
                                  strand_ambig_action=strand_ambig_action_wgt,
                                  ld_Rinfo=ld_Rinfo_chrom,
                                  R_wgt=R_wgt,
                                  wgt=wgt)
    
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
      weight_table_current <- data.frame(gene=gname,
                                         rsid=snps$id,
                                         varID=rsid_varID$varID[match(snps$id, rsid_varID$rsid)],
                                         ref_allele=snps$ref,
                                         eff_allele=snps$alt,
                                         weight=wgt[,"weight"])
      
      weight_table_harmonized[[gname]] <- weight_table_current
    }
  }
  
  weight_table_harmonized <- do.call(rbind, weight_table_harmonized)
  extra_table <- extra_table[extra_table$gene %in% weight_table_harmonized$gene,]
  
  if (file.exists(paste0(outputdir, outname, ".db"))){
    invisible(file.remove(paste0(outputdir, outname, ".db")))
  }
  
  db <- RSQLite::dbConnect(sqlite, paste0(outputdir, outname, ".db"))
  RSQLite::dbWriteTable(db, "extra", extra_table)
  RSQLite::dbWriteTable(db, "weights", weight_table_harmonized)
  RSQLite::dbDisconnect(db)
  
  invisible(file.remove(paste0(outputdir, outname, "_ld_R_chr", 1:22, ".txt")))
}
