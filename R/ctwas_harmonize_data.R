
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
  flip2 = flip;

  snp = list()
  snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
  snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F
  snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F
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
harmonize_z_ld <- function(z_snp, ld_snpinfo){
    snpnames <- intersect(z_snp$id, ld_snpinfo$id)
    if (length(snpnames) != 0){
      z.idx <- match(snpnames, z_snp$id)
      ld.idx <- match(snpnames, ld_snpinfo$id)
      qc <- allele.qc(z_snp[z.idx, ]$A1, z_snp[z.idx, ]$A2,
                            ld_snpinfo[ld.idx, ]$alt, ld_snpinfo[ld.idx, ]$ref)

      ifflip <- qc[["flip"]]
      ifremove <- !qc[["keep"]]

      flip.idx <- z.idx[ifflip]
      z_snp[flip.idx, c("A1", "A2")] <- z_snp[flip.idx, c("A2", "A1")]
      z_snp[flip.idx, "z"] <- - z_snp[flip.idx, "z"]

      remove.idx <- z.idx[ifremove]
      if (length(remove.idx) != 0) {
        z_snp <- z_snp[-remove.idx, ,drop = F]
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
harmonize_wgt_ld <- function(wgt.matrix, snps, ld_snpinfo){
  # `snps` from FUSION follows .bim format
  colnames(snps) <- c("chrom", "id", "cm", "pos", "alt", "ref")
  snps <- snps[match(rownames(wgt.matrix), snps$id),]
  snpnames <- intersect(snps$id, ld_snpinfo$id)
  if (length(snpnames) != 0) {
    snps.idx <- match(snpnames, snps$id)
    ld.idx <- match(snpnames, ld_snpinfo$id)
    qc <- allele.qc(snps[snps.idx, ]$alt, snps[snps.idx, ]$ref,
                          ld_snpinfo[ld.idx, ]$alt, ld_snpinfo[ld.idx, ]$ref)
    ifflip <- qc[["flip"]]
    ifremove <- !qc[["keep"]]

    flip.idx <- snps.idx[ifflip]
    snps[flip.idx, c("alt", "ref")] <-  snps[flip.idx, c("ref", "alt")]
    wgt.matrix[flip.idx, ] <- - wgt.matrix[flip.idx, ]

    remove.idx <- snps.idx[ifremove]
    if (length(remove.idx) != 0) {
      snps <- snps[-remove.idx, ,drop = F]
      wgt.matrix <- wgt.matrix[-remove.idx, , drop = F]
    }
  }
  return(list("wgt" = wgt.matrix, "snps" = snps))
}

