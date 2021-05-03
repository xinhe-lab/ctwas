#' Get gene and SNP index for each region
#' @description For each region, get the index for snp and gene (index
#' is location/column number in .pgen file or .expr file) located within
#' this region.
#'
#' @param regionfile regions file. Has three columns: chr, start, end. The regions file
#' should provide non overlaping regions defining LD blocks. currently does not support
#' chromsome X/Y etc.
#'
#' @param select Default is NULL, all variants will be selected. Or a vector of variant IDs,
#' or a data frame with columns id and z (id is for gene or SNP id, z is for z scores).
#' z will be used for remove SNPs if the total number of SNPs exceeds limit. See
#' parameter `maxSNP` for more information.
#'
#' @param maxSNP Default is Inf, no limit for the maximum number of SNPs in a region. Or an
#' integer indicating the maximum number of SNPs allowed in a region. This
#' parameter is useful when a region contains many SNPs and you don't have enough memory to
#' run the program. In this case, you can put a limit on the number of SNPs in the region.
#' If z scores are given in the parameter `select`, i.e. a data frame with columns id and z is
#' provided, SNPs are ranked based on |z| from high to low and only the top `maxSNP` SNPs
#' are kept. If only variant ids are provided, then `maxSNP` number of SNPs will be chosen
#' randomly.
#'
#' @param minvar minimum number of variants in a region.
#'
#' @param merge True/False. If merge regions when a gene's eQTL locate in multiple regions.
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
                          maxSNP = Inf,
                          minvar = 1,
                          merge = T) {

  reg <- read.table(regionfile, header = T, stringsAsFactors = F)
  if (is.character(reg$chr)){
    reg$chr <- readr::parse_number(reg$chr)
  }
  # sort regions
  reg <- reg[order(reg$chr, reg$start),]

  loginfo("No. LD regions: %s", nrow(reg))

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

    regions <- reg[reg$chr == b, ]

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
                    & snpinfo$keep == 1)

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

          gidx <-  unique(c(previous[["gidx"]], current[["gidx"]]))
          sidx <-  unique(c(previous[["sidx"]], current[["sidx"]]))

          gid <- geneinfo[gidx, "id"]
          sid <- snpinfo[sidx, "id"]

          regionlist[[b]][[as.character(rn)]] <- list("gidx" = gidx,
                                                      "gid"  = gid,
                                                      "sidx" = sidx,
                                                      "sid"  = sid,
                                                      "start" =  regionlist[[b]][[as.character(rn - 1)]]$start,
                                                      "stop" = regions[rn, "stop"])

          regionlist[[b]][[as.character(rn -1)]] <- NULL
        }
      }
    }

    loginfo("No. regions with at least one SNP/gene for chr%s after merging: %s",
            b, length(regionlist[[b]]))
  }


  loginfo("Trim regions with SNPs more than %s", maxSNP)

  if ("z" %in% colnames(select)) {
    # z score is given, trim snps with lower |z|
    for (b in 1: length(regionlist)){
      for (rn in names(regionlist[[b]])) {
        if (length(regionlist[[b]][[rn]][["sid"]]) > maxSNP){
            z.abs <- abs(select[match(regionlist[[b]][[rn]][["sid"]], select[, "id"]), "z"])
            ifkeep <- rank(-z.abs, ties.method = "first") <= maxSNP
            regionlist[[b]][[rn]][["sidx"]] <-  regionlist[[b]][[rn]][["sidx"]][ifkeep]
            regionlist[[b]][[rn]][["sid"]] <-  regionlist[[b]][[rn]][["sid"]][ifkeep]
            regionlist[[b]][[rn]][["prop"]] <- maxSNP/length(ifkeep)
        }
      }
    }
  } else{
    # if no z score information, randomly select snps
    for (b in 1: length(regionlist)){
      for (rn in names(regionlist[[b]])) {
        if (length(regionlist[[b]][[rn]][["sid"]]) > maxSNP){
          n.ori <- length(regionlist[[b]][[rn]][["sid"]])
          ifkeep <- rep(F, n.ori)
          set.seed <- 99
          ifkeep[sample.int(n.ori, size = maxSNP)] <- T
          regionlist[[b]][[rn]][["sidx"]] <-  regionlist[[b]][[rn]][["sidx"]][ifkeep]
          regionlist[[b]][[rn]][["sid"]] <-  regionlist[[b]][[rn]][["sid"]][ifkeep]
          regionlist[[b]][[rn]][["prop"]] <- maxSNP/length(ifkeep)
        }
      }
    }
  }


  return(regionlist)



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
