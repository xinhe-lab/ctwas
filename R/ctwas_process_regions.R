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
index_regions <- function(regionfile,
                          exprvarfs,
                          pvarfs = NULL,
                          ld_Rfs = NULL,
                          select = NULL,
                          thin = 1,
                          maxSNP = Inf,
                          minvar = 1,
                          merge = T,
                          outname = NULL,
                          outputdir = getwd()) {

  if (is.null(pvarfs) & is.null(ld_Rfs)){
    stop("Stopped: missing LD/genotype information.
         LD/genotype information needs to be provided either in genotype form
         (see parameter description for pvarfs) or R matrix form
         (see parameter description for ld_Rf) ")
  }

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

  if (is.null(dim(select))) {
    selectid <- select
  } else {
    selectid <- select$id
  }

  regionlist <- list()
  for (b in 1:length(exprvarfs)){

    if (!is.null(pvarfs)){
      # get snp info (from pvarf file)
      pvarf <- pvarfs[b]
      snpinfo <- read_pvar(pvarf)
    } else {
      # get snp info (from LD R matrix file)
      ld_Rf <- ld_Rfs[b]
      snpinfo <- read_ld_Rvar(ld_Rf)
    }

    if (unique(snpinfo$chrom) != b){
      stop("Input genotype file not splitted by chromosome or not in correct order")
    }

    # select variant
    snpinfo$keep <- 1
    if (!is.null(selectid)){
      snpinfo[!(snpinfo$id %in% selectid), "keep"] <- 0
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
      if (!is.null(selectid)){
        geneinfo[!(geneinfo$id %in% selectid), "keep"] <- 0
      }
    }

    regions <- reg[reg$chr == b, ]

    regionlist[[b]] <- list()

    for (rn in 1:nrow(regions)){
      rn.start <- regions[rn, "start"]
      rn.stop <- regions[rn, "stop"]

      if (isTRUE(merge)){
        gidx <- which(geneinfo$chrom == b & geneinfo$p1 >= rn.start & geneinfo$p0 < rn.stop
                    & geneinfo$keep == 1) # allow overlap
      } else {
        gidx <- which(geneinfo$chrom == b & geneinfo$p0 >= rn.start & geneinfo$p0 < rn.stop
                      & geneinfo$keep == 1) # unique assignment to regions
      }

      sidx <- which(snpinfo$chrom == b & snpinfo$pos >= rn.start & snpinfo$pos < rn.stop
                    & snpinfo$keep == 1 & snpinfo$thin_tag == 1)

      if (length(gidx) + length(sidx) < minvar) {next}

      gid <- geneinfo[gidx, "id"]
      sid <- snpinfo[sidx, "id"]

      minpos <- min(c(geneinfo[gidx, "p0"], snpinfo[sidx, "pos"]))
      maxpos <- max(c(geneinfo[gidx, "p1"], snpinfo[sidx, "pos"]))

      regionlist[[b]][[as.character(rn)]] <- list("gidx" = gidx,
                                    "gid"  = gid,
                                    "sidx" = sidx,
                                    "sid"  = sid,
                                    "start" = rn.start,
                                    "stop" = rn.stop,
                                    "minpos" = minpos,
                                    "maxpos" = maxpos)
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
                                                        "start" = previous$start,
                                                        "stop" = current$stop,
                                                        "minpos" = previous$minpos,
                                                        "maxpos" = current$maxpos)

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
          idx <- match(regionlist[[b]][[rn]][["sid"]], select[, "id"])
          z.abs <- abs(select[idx, "z"])
          ifkeep <- rank(-z.abs) <= maxSNP
          regionlist[[b]][[rn]][["sidx"]] <-  regionlist[[b]][[rn]][["sidx"]][ifkeep]
          regionlist[[b]][[rn]][["sid"]] <-  regionlist[[b]][[rn]][["sid"]][ifkeep]
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
        }
      }
    }
  }

  if (is.null(pvarfs)){

    loginfo("Adding R matrix info, as genotype is not given")

    dir.create(file.path(outputdir, paste0(outname, "_LDR")), showWarnings = F)

    wgtall <- lapply(exprvarfs, function(x){
      load(paste0(strsplit(x, ".exprvar")[[1]], ".exprqc.Rd")); wgtlist})
    wgtlistall <- do.call(c, wgtall)
    names(wgtlistall) <- do.call(c, lapply(wgtall, names))

    for (b in 1: length(ld_Rfs)){
      loginfo("Adding R matrix info for chrom %s", b)
      ld_Rf <- ld_Rfs[b]
      ld_Rinfo <- data.table::fread(ld_Rf, header = T)
      for (rn in names(regionlist[[b]])){
        ifreg <- ifelse(regionlist[[b]][[rn]][["minpos"]] < ld_Rinfo[, "stop"]
                        & regionlist[[b]][[rn]][["maxpos"]] >= ld_Rinfo[, "start"], T, F)
        regRDS <- ld_Rinfo[ifreg, "RDS_file"]
        R_snp <- lapply(regRDS, readRDS)
        R_snp <- as.matrix(Matrix::bdiag(R_snp))
        R_snp_anno <- do.call(rbind, lapply(regRDS, read_ld_Rvar_RDS))

        #update sidx to match R matrix info
        sidx <-  match(regionlist[[b]][[rn]][["sid"]], R_snp_anno$id)
        regionlist[[b]][[rn]][["sidx"]] <- sidx

        gnames <- regionlist[[b]][[rn]][["gid"]]
        R_snp_gene <- matrix( , nrow(R_snp), length(gnames))
        R_gene <- diag(length(gnames))
        if (length(gnames) > 0) {

          ldr <- list()
          for (i in 1:length(gnames)){
            gname <- gnames[i]
            wgt <- wgtlistall[[gname]]
            snpnames <- rownames(wgt)
            ld.idx <- match(snpnames, R_snp_anno$id)
            ldr[[gname]] <- ld.idx
            R.s <- R_snp[ld.idx, ld.idx]
            R_snp_gene[,i] <- sapply(1:nrow(R_snp),
                                     function(x){crossprod(wgt,R_snp[ld.idx,x])/sqrt(crossprod(wgt,R.s)%*%wgt*R_snp[x,x])})
          }

          if (length(gnames) > 1){
            gene_pairs <- combn(length(gnames), 2)
            wgtr <- wgtlistall[gnames]
            gene_corrs <- apply(gene_pairs, 2, function(x){t(wgtr[[x[1]]])%*%R_snp[ldr[[x[1]]], ldr[[x[2]]]]%*%wgtr[[x[2]]]/(
              sqrt(t(wgtr[[x[1]]])%*%R_snp[ldr[[x[1]]], ldr[[x[1]]]]%*%wgtr[[x[1]]]) *
                sqrt(t(wgtr[[x[2]]])%*%R_snp[ldr[[x[2]]], ldr[[x[2]]]]%*%wgtr[[x[2]]]))})
            R_gene[t(gene_pairs)] <- gene_corrs
            R_gene[t(gene_pairs[c(2,1),])] <- gene_corrs
          }
        }

        regionlist[[b]][[rn]][["regRDS"]] <- regRDS

        R_sg_file <- file.path(outputdir, paste0(outname, "_LDR"), paste0("chr", b, "_reg", rn, ".R_snp_gene.RDS"))
        R_g_file <- file.path(outputdir, paste0(outname, "_LDR"), paste0("chr", b, "_reg", rn, ".R_gene.RDS"))
        saveRDS(R_snp_gene, file=R_sg_file)
        saveRDS(R_gene, file=R_g_file)
        regionlist[[b]][[rn]][["R_sg_file"]] <- R_sg_file
        regionlist[[b]][[rn]][["R_g_file"]] <- R_g_file

        if (thin < 1){
          R_s_file <- file.path(outputdir, paste0(outname, "_LDR"), paste0("chr", b, "_reg", rn, ".R_snp.RDS"))
          R_snp <- R_snp[sidx, sidx, drop = F]
          saveRDS(R_snp, file=R_s_file)
          regionlist[[b]][[rn]][["R_s_file"]] <- R_s_file
        }
      }
    }
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

