
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
#' @param z_snp a data frame, with columns "id", "A1", "A2" and "z".
#'     Z scores for every SNP. "A1" is the effect allele.
#' @param ld_snpinfo A data frame, snp info for LD reference,
#'  with columns "chrom", "id", "pos", "alt", "ref".
#' @return A data frame, z_snp with the "z" columns flipped to match LD ref.
#'
#' @export
harmonize_z_ld <- function(z_snp, ld_snpinfo, recover_strand_ambig = T, ld_pgenfs = NULL, ld_Rinfo = NULL){
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
    if (recover_strand_ambig & any(ifremove)){
      #compare sign of imputed z score with observed z score for strand ambiguous variants 
      #following imputation strategy in https://dx.doi.org/10.1093%2Fbioinformatics%2Fbtu416
      if (is.null(ld_pgenfs)){
        loginfo("harmonizing strand-ambiguous z scores using imputation by region")
        for (i in 1:nrow(ld_Rinfo)){
          R_snp <- readRDS(ld_Rinfo$RDS_file[i])
          R_snp_anno <- read_ld_Rvar_RDS(ld_Rinfo$RDS_file[i])
          
          #indicator if z.idx (on chromsome, in LD ref) is in current region
          z.idx.ifreg <- z_snp$id[z.idx] %in% R_snp_anno$id
          
          #index if z.idx (on chromsome, in LD ref) is in current region and unambiguous/ambiguous
          z.idx.unambig <- z.idx[z.idx.ifreg & qc$keep]
          z.idx.ambig <- z.idx[z.idx.ifreg & !qc$keep]
          
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
          if_sign_neq <- sign(z_i_obs) != sign(z_i)
          z_snp[z.idx.ambig[if_sign_neq], "z"] <- -z_snp[z.idx.ambig[if_sign_neq], "z"]
        }
      } else {
        #TO-DO: mirror previous section but compute R for each region using X
      }
    } else if (any(ifremove)) {
      remove.idx <- z.idx[ifremove]
      z_snp <- z_snp[-remove.idx, , drop = F]
    }
  }
  return(z_snp)
}

#' Harmonize z scores from GWAS to match ld reference genotypes.
#' Flip signs when reverse complement matches, remove strand ambiguous SNPs
#' @param wgt.matrix from FUSION weight .Rdat file
#' @param snps from FUSION weight .Rdat file
#'  with columns "chrom", "id", "pos", "alt", "ref". The effect allele
#'  for FUSION is alt.
#' @return wgt.matrix and snps with alleles flipped to match
#'
#' @export
harmonize_wgt_ld <- function (wgt.matrix, snps, ld_snpinfo, recover_strand_ambig=T, 
                              ld_pgenfs=NULL, ld_Rinfo=NULL, R_wgt_all=NULL, gname=NULL, wgt=NULL){
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
    
    if (recover_strand_ambig & 
        any(ifremove) & 
        sum(!ifremove)>0){
      if (is.null(ld_pgenfs)){
        #load correlation matrix(es) for LD reference(s) containing current weight
        wgt_pos <- ld_snpinfo$pos[ld_snpinfo$id %in% snpnames]
        regnames <- unique(sapply(wgt_pos, function(x){which(x > ld_Rinfo$start & x < ld_Rinfo$stop)}))
        regRDS <- ld_Rinfo[match(regnames, ld_Rinfo$region_name), "RDS_file"]
        R_snp <- lapply(regRDS, readRDS)
        R_snp <- as.matrix(Matrix::bdiag(R_snp))
        R_snp_anno <- do.call(rbind, lapply(regRDS, read_ld_Rvar_RDS))
        
        #index the variant positions in LD reference
        R_snp.idx <- match(snpnames, R_snp_anno$id)
        R_snp.idx.unambig <- R_snp.idx[!ifremove]
        
        #subset R_wgt_all to current weight
        R_wgt <- R_wgt_all[R_wgt_all$GENE == gname,]
        R_wgt <- R_wgt[R_wgt$RSID1!=R_wgt$RSID2,] #discard variant correlations with itself (NOTE: not equal to one? not scaled?)
        
        #drop R_wgt correlations between ambiguous variants
        R_wgt <- R_wgt[R_wgt$RSID1 %in% wgt$varID[!ifremove] | R_wgt$RSID2 %in% wgt$varID[!ifremove],]
        
        #flip correlations if weights were flipped in previous step
        for (i in flip.idx){
          ifflip_rwgt <- R_wgt$RSID1 == wgt$varID[i] | R_wgt$RSID2 == wgt$varID[i]
          R_wgt$VALUE[ifflip_rwgt] <- -R_wgt$VALUE[ifflip_rwgt]
        }
        
        #iterate over ambiguous snps
        for (i in snpnames[ifremove]){
          #index current ambiguous snp
          snpnames.idx <- match(i, snpnames)
          
          #sum of correlations in LD reference between current ambiguous variant and unambiguous variants
          sumcor_R_snp <- sum(R_snp[R_snp.idx[snpnames.idx], R_snp.idx.unambig])
          
          #sum of correlations in weights between current ambiguous variant and unambiguous variants
          sumcor_R_wgt <- sum(R_wgt$VALUE[R_wgt$RSID1==wgt$varID[snpnames.idx] | R_wgt$RSID2==wgt$varID[snpnames.idx]])
          
          #flip weight if sign of correlations is not the same
          if (sign(sumcor_R_snp)!=sign(sumcor_R_wgt)){
            wgt.matrix[snpnames.idx,] <- -wgt.matrix[snpnames.idx,]
          }
        }
      } else {
        #TO-DO: mirror following section but compute R_snp each region using Xs
      }
    } else if (recover_strand_ambig & 
               any(ifremove) & 
               sum(ifremove)==1){
      #take no action if single variant. wrote this as separate if-statement for clarity, but it could be rolled into the following if-statement
    } else if (any(ifremove)){
      #if recover_strand_ambig=F, or >2 ambiguous variants and 0 unambiguous variants, discard the ambiguous variants
      remove.idx <- snps.idx[ifremove]
      snps <- snps[-remove.idx, , drop = F]
      wgt.matrix <- wgt.matrix[-remove.idx, , drop = F]
    }
  }
  return(list(wgt = wgt.matrix, snps = snps))
}

