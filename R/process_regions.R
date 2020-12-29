#' Get gene and SNP index for each region
#' @description For each region, get the index for snp and gene (index
#' is location/column number in .pgen file or .expr file) located within
#' this region.
#'
#' @param regionfile Has three columns, chr, start, end
#' @param down_sample_ratio  A scalar in (0,1]. The proportion of SNPs
#'  left after down sampling.
#'
#' @return A list. Items correspond to each pvarf/exprvarf. Each Item is
#'  also a list, the items in this list are for each region.
#'
#' @importFrom logging loginfo
#'
index_regions <- function(pvarfs,
                          exprvarfs,
                          regionfile,
                          down_sample_ratio = 1) {

  reg <- read.table(regionfile, header = T, stringsAsFactors = F)

  loginfo("No. LD regions: %s", nrow(reg))

  if (down_sample_ratio > 1 | down_sample_ratio <= 0){
    stop("down_sample_ratio needs to be in (0,1]")
  }

  regionlist <- list()
  for (b in 1:22){
    # get snp info (from pvarf file)
    pvarf <- pvarfs[b]
    snpinfo <- read_pvar(pvarf)

    snpinfo$down_sample_tag <- 0
    nkept <- round(nrow(snpinfo) * down_sample_ratio)
    seed(99)
    snpinfo$down_sample_tag[sample(1:nrow(snpinfo), nkept)] <- 1


    if (unique(snpinfo$chrom) != b){
      stop("Input genotype not by chromosome or not in correct order")
    }

    # get gene info (from exprf file)
    exprvarf <- exprvarfs[b]
    geneinfo <- read_exprvar(exprvarf)

    if (unique(geneinfo$chrom) != b){
      stop("Imputed expression not by chromosome or not in correct order")
    }

    regions <- reg[reg$chr == b | reg$chr == paste0("chr", b), ]

    regionlist[[b]] <- list()

    for (rn in 1:nrow(regions)){
      p0 <- regions[rn, "start"]
      p1 <- regions[rn, "stop"]

      gidx <- which(geneinfo$chrom == b & geneinfo$p0 > p0 & geneinfo$p0 < p1)
      sidx <- which(snpinfo$chrom == b & snpinfo$pos > p0 & snpinfo$pos < p1
                    & snpinfo$down_sample_tag == 1)

      if (length(gidx) == 0 & length(sidx) == 0) {next}

      regionlist[[b]][[rn]][["gidx"]] <-  gidx
      regionlist[[b]][[rn]][["sidx"]] <- sidx

    }
    loginfo("No. regions with at least one SNP/gene for chr%s: %s",
            b, length(regionlist[[b]]))
  }

  regionlist
}

#' filter regions based on probality of at most 1 causal effect
filter_regions <- function(regionlist, group_prior, prob_single = 0.8){
  prior.gene <- group_prior[1]
  prior.SNP <- group_prior[2]

  for (b in 1:22){

    for (rn in names(regionlist[[b]])) {

      gidx <- regionlist[[b]][[rn]][["gidx"]]
      sidx <- regionlist[[b]][[rn]][["sidx"]]
      p.g <- length(gidx)
      p.s <- length(sidx)
      P2 <- 1 - (1- prior.gene)**p.g * (1- prior.SNP)**p.s -
        p.g * prior.gene * (1 - prior.gene) ** (p.g - 1) * (1- prior.SNP) ** p.s -
        p.s * (1 - prior.gene) ** p.g * prior.SNP * (1- prior.SNP) ** (p.s - 1)

      if (P2 >= 1-prob_single){
        regionlist[[b]][[rn]] <- NULL
      }
    }
  }
  regionlist
}
