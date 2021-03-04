#' Get gene and SNP index for each region
#' @description For each region, get the index for snp and gene (index
#' is location/column number in .pgen file or .expr file) located within
#' this region.
#'
#' @param regionfile regions file. Has three columns: chr, start, end. The regions file
#' should provide non overlaping regions defining LD blocks. currently does not support
#' chromsome X/Y etc.
#'
#' @param select variant ids to include. If NULL, all variants will be selected
#' @param thin  A scalar in (0,1]. The proportion of SNPs
#'  left after down sampling. Only applied on SNPs after selecting variants.
#' @param minvar minimum number of variatns in a region
#'
#' @param merge True/False. If merge regions when a gene belong to multiple regions.
#'
#' @return A list. Items correspond to each pvarf/exprvarf. Each Item is
#'  also a list, the items in this list are for each region.
#'
#' @importFrom logging loginfo
#'
index_regions <- function(pvarfs,
                          exprvarfs,
                          regionfile,
                          select = NULL,
                          thin = 1,
                          minvar = 1,
                          merge = T) {

  reg <- read.table(regionfile, header = T, stringsAsFactors = F)
  if (is.character(reg$chr)){
    reg$chr <- readr::parse_number(reg$chr)
  }
  # sort regions
  reg <- reg[order(reg$chr, reg$start),]

  loginfo("No. LD regions: %s", nrow(reg))

  if (thin > 1 | thin <= 0){
    stop("thin needs to be in (0,1]")
  }

  regionlist <- list()
  for (b in 1:length(pvarfs)){
    # get snp info (from pvarf file)
    pvarf <- pvarfs[b]
    snpinfo <- read_pvar(pvarf)

    if (unique(snpinfo$chrom) != b){
      stop("Input genotype file not splitted by chromosome or not in correct order")
    }

    # select variant
    snpinfo$keep <- 1
    if (!is.null(select)){
      snpinfo[!(snpinfo$id %in% select), "keep"] <- 0
    }

    # downsampling for SNPs
    snpinfo$thin_tag <- 0
    nkept <- round(nrow(snpinfo) * thin)
    set.seed(99)
    snpinfo$thin_tag[sample(1:nrow(snpinfo), nkept)] <- 1

    # get gene info (from exprf file)
    exprvarf <- exprvarfs[b]
    geneinfo <- read_exprvar(exprvarf)

    if (nrow(geneinfo)!=0){
      if (unique(geneinfo$chrom) != b){
        stop("Imputed expression not by chromosome or not in correct order")
      }

      # select variant
      geneinfo$keep <- 1
      if (!is.null(select)){
        geneinfo[!(geneinfo$id %in% select), "keep"] <- 0
      }
    }


    regions <- reg[reg$chr == b | reg$chr == paste0("chr", b), ]

    regionlist[[b]] <- list()

    for (rn in 1:nrow(regions)){
      p0 <- regions[rn, "start"]
      p1 <- regions[rn, "stop"]

      if (isTRUE(merge)){
        gidx <- which(geneinfo$chrom == b & geneinfo$p1 > p0 & geneinfo$p0 < p1
                    & geneinfo$keep == 1) # allow overlap
      } else {
        gidx <- which(geneinfo$chrom == b & geneinfo$p0 > p0 & geneinfo$p0 < p1
                      & geneinfo$keep == 1) # unique assignment to regions
      }

      sidx <- which(snpinfo$chrom == b & snpinfo$pos > p0 & snpinfo$pos < p1
                    & snpinfo$keep == 1 & snpinfo$thin_tag == 1)

      gid <- geneinfo[gidx, "id"]
      sid <- snpinfo[sidx, "id"]

      if (length(gidx) + length(sidx) < minvar) {next}

      regionlist[[b]][[as.character(rn)]] <- list("gidx" = gidx,
                                    "gid"  = gid,
                                    "sidx" = sidx,
                                    "sid"  = sid,
                                    "start" = p0,
                                    "stop" = p1)

    }
    loginfo("No. regions with at least one SNP/gene for chr%s: %s",
            b, length(regionlist[[b]]))

    if (isTRUE(merge) & nrow(regions) >= 2){
      loginfo("Merge regions for chr%s", b)
      for (rn in 2:nrow(regions)){
        current <- regionlist[[b]][[as.character(rn)]]
        previous <- regionlist[[b]][[as.character(rn - 1)]]

          if (length(intersect(current[["gid"]], previous[["gid"]]))> 0){

            merged <- lapply(names(current), function(x) unique(c(previous[[x]],
                                                                current[[x]])))

            names(merged) <- names(current)
            merged[["start"]] <- regions[rn - 1, "start"]
            merged[["stop"]] <- regions[rn, "stop"]

            regionlist[[b]][[as.character(rn)]] <- merged

            regionlist[[b]][[as.character(rn -1)]] <- NULL
          }
      }
    }

    loginfo("No. regions with at least one SNP/gene for chr%s after merging: %s",
            b, length(regionlist[[b]]))
  }

  regionlist
}



#' filter regions based on probality of at most 1 causal effect
filter_regions <- function(regionlist, group_prior, prob_single = 0.8){
  prior.gene <- group_prior[1]
  prior.SNP <- group_prior[2]

  regionlist2 <- regionlist
  for (b in 1: length(regionlist)){

    for (rn in names(regionlist[[b]])) {

      gidx <- regionlist[[b]][[rn]][["gidx"]]
      sidx <- regionlist[[b]][[rn]][["sidx"]]
      p.g <- length(gidx)
      p.s <- length(sidx)
      P2 <- 1 - (1- prior.gene)**p.g * (1- prior.SNP)**p.s -
        p.g * prior.gene * (1 - prior.gene) ** (p.g - 1) * (1- prior.SNP) ** p.s -
        p.s * (1 - prior.gene) ** p.g * prior.SNP * (1- prior.SNP) ** (p.s - 1)

      if (P2 >= 1-prob_single){
        regionlist2[[b]][[rn]] <- NULL
      }
    }
  }
  regionlist2
}

#' parallel regions
#' @param ncore integer, numeber of cores, at least 1
#' regions allocated to given number of cores
#' regionlist need to contain at least 1 non-empty
region2core <- function(regionlist, ncore = 1){
  dflist <- list()
  for (b in 1:length(regionlist)){
    if (length(regionlist[[b]]) > 0){
      dflist[[b]] <- data.frame("b" = b, "rn"= names(regionlist[[b]]), stringsAsFactors = FALSE)
    }
  }
  df <- do.call(rbind, dflist)
  if (ncore > 1) {
    d <- cut(1:nrow(df), ncore, labels = FALSE)
    corelist <- split(df,d)
  } else {
    corelist <- list("1" = df)
  }
  return(corelist)
}

