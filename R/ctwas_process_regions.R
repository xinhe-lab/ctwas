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
                          outputdir = getwd(),
                          ncore = 1,
                          reuse_R_gene = F) {

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

    if (isTRUE(unique(snpinfo$chrom) != b)){
      stop("Input genotype file not split by chromosome or not in correct order")
    }

    # select variant
    snpinfo$keep <- rep(1, nrow(snpinfo))
    if (!is.null(selectid)){
      snpinfo$keep[!(snpinfo$id %in% selectid)] <- 0
    }

    # downsampling for SNPs
    snpinfo$thin_tag <- rep(0, nrow(snpinfo))
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
      
      #loginfo(rn)
      
      rn.start <- regions$start[rn]
      rn.stop <- regions$stop[rn]

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

      gid <- geneinfo$id[gidx]
      sid <- snpinfo$id[sidx]

      minpos <- min(c(geneinfo$p0[gidx], snpinfo$pos[sidx]))
      maxpos <- max(c(geneinfo$p1[gidx], snpinfo$pos[sidx]))

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

            gid <- geneinfo$id[gidx]
            sid <- snpinfo$id[sidx]

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
    
    wgtall <- lapply(exprvarfs, function(x){load(paste0(strsplit(x, ".exprvar")[[1]], ".exprqc.Rd")); wgtlist})
    wgtlistall <- do.call(c, wgtall)
    names(wgtlistall) <- do.call(c, lapply(wgtall, names))
    rm(wgtall)
    
    regionlist_all <- list()
    
    for (b in 1:22){
      regionlist_all[[b]] <- cbind(b, names(regionlist[[b]]))
    }
    
    regionlist_all <- as.data.frame(do.call(rbind, regionlist_all))
    colnames(regionlist_all) <- c("b", "rn")
    regionlist_all$b <- as.integer(regionlist_all$b)
    
    corelist <- lapply(1:ncore, function(core){njobs <- ceiling(nrow(regionlist_all)/ncore); jobs <- ((core-1)*njobs+1):(core*njobs); jobs[jobs<=nrow(regionlist_all)]})
    names(corelist) <- 1:ncore
    
    
    cl <- parallel::makeCluster(ncore, outfile = "")
    doParallel::registerDoParallel(cl)
    
    outlist <- foreach(core = 1:ncore, .combine = "c", .packages = c("ctwas", "tools")) %dopar% {
      
      regionlist_core <- regionlist_all[corelist[[core]],]
      outlist_core <- list()
      
      b_core <- unique(regionlist_core$b)
      
      for (b in b_core){
        logging::loginfo("Adding R matrix info for chrom %s (core %s)", b, core)
        
        ld_Rf <- ld_Rfs[b]
        ld_Rinfo <- as.data.frame(data.table::fread(ld_Rf, header = T))
        
        regionlist_core_b <- regionlist_core$rn[regionlist_core$b==b]
        
        for (rn in regionlist_core_b){
          outlist_core_region <- list(b=b, rn=rn)
          
          ifreg <- ifelse(regionlist[[b]][[rn]][["start"]] < ld_Rinfo[, "stop"] & regionlist[[b]][[rn]][["stop"]] >= ld_Rinfo[, "start"], T, F)
          
          regRDS <- ld_Rinfo[ifreg, "RDS_file"]
          R_snp <- lapply(regRDS, readRDS)
          
          if (length(R_snp)==1){
            R_snp <- unname(R_snp[[1]])
          } else {
            R_snp <- suppressWarnings(as.matrix(Matrix::bdiag(R_snp)))
          }
          
          R_snp_anno <- as.data.frame(do.call(rbind, lapply(regRDS, read_ld_Rvar_RDS)))
          sidx <-  match(regionlist[[b]][[rn]][["sid"]], R_snp_anno$id)
          outlist_core_region[["sidx"]] <- sidx
          
          gnames <- regionlist[[b]][[rn]][["gid"]]
          R_snp_gene <- matrix(NA, nrow(R_snp), length(gnames))
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
              R_snp_gene[,i] <- sapply(1:nrow(R_snp), function(x){t(wgt)%*%R_snp[ld.idx,x]/sqrt(t(wgt)%*%R.s%*%wgt*R_snp[x,x])})
            }
            
            if (length(gnames) > 1 & !reuse_R_gene){
              gene_pairs <- combn(length(gnames), 2)
              wgtr <- wgtlistall[gnames]
              
              gene_corrs <- apply(gene_pairs, 2, function(x){t(wgtr[[x[1]]])%*%R_snp[ldr[[x[1]]], ldr[[x[2]]]]%*%wgtr[[x[2]]]/(
                sqrt(t(wgtr[[x[1]]])%*%R_snp[ldr[[x[1]]], ldr[[x[1]]]]%*%wgtr[[x[1]]]) *
                  sqrt(t(wgtr[[x[2]]])%*%R_snp[ldr[[x[2]]], ldr[[x[2]]]]%*%wgtr[[x[2]]]))})
              R_gene[t(gene_pairs)] <- gene_corrs
              R_gene[t(gene_pairs[c(2,1),])] <- gene_corrs
            }
          }
          
          #remove genes with NA in R_snp_gene (e.g. boundary spanning with merge=F) from analysis
          gene_not_NA <- which(apply(R_snp_gene, 2, function(x){!any(is.na(x))}))
          outlist_core_region[["gidx"]] <- regionlist[[b]][[rn]][["gidx"]][gene_not_NA]
          outlist_core_region[["gid"]] <- regionlist[[b]][[rn]][["gid"]][gene_not_NA]
          R_snp_gene <- R_snp_gene[,gene_not_NA,drop=F]
          R_gene <- R_gene[gene_not_NA,gene_not_NA,drop=F]
          
          #save R_snp_gene and R_gene
          outlist_core_region[["regRDS"]] <- regRDS
          R_sg_file <- file.path(outputdir, paste0(outname, "_LDR"), paste0("chr", b, "_reg", rn, ".R_snp_gene.RDS"))
          R_g_file <- file.path(outputdir, paste0(outname, "_LDR"), paste0("chr", b, "_reg", rn, ".R_gene.RDS"))
          saveRDS(R_snp_gene, file=R_sg_file)
          
          if (!reuse_R_gene){
            saveRDS(R_gene, file=R_g_file)
          }
          
          outlist_core_region[["R_sg_file"]] <- R_sg_file
          outlist_core_region[["R_g_file"]] <- R_g_file
          
          if (thin < 1){
            R_s_file <- file.path(outputdir, paste0(outname, "_LDR"), paste0("chr", b, "_reg", rn, ".R_snp.RDS"))
            R_snp <- R_snp[sidx, sidx, drop = F]
            saveRDS(R_snp, file=R_s_file)
            outlist_core_region[["R_s_file"]] <- R_s_file
          }
          
          outlist_core[[length(outlist_core)+1]] <- outlist_core_region
        }
      }
      outlist_core
    }

    parallel::stopCluster(cl)
    
    for (i in 1:length(outlist)){
      b <- outlist[[i]][["b"]]
      rn <- outlist[[i]][["rn"]]
      regionlist[[b]][[rn]][["sidx"]] <- outlist[[i]][["sidx"]]
      regionlist[[b]][[rn]][["gidx"]] <- outlist[[i]][["gidx"]]
      regionlist[[b]][[rn]][["gid"]] <- outlist[[i]][["gid"]]
      regionlist[[b]][[rn]][["regRDS"]] <- outlist[[i]][["regRDS"]]
      regionlist[[b]][[rn]][["R_sg_file"]] <- outlist[[i]][["R_sg_file"]]
      regionlist[[b]][[rn]][["R_g_file"]] <- outlist[[i]][["R_g_file"]]
      regionlist[[b]][[rn]][["R_s_file"]] <- outlist[[i]][["R_s_file"]]
    }
  }

  regionlist
}


#' filter regions based on probability of at most 1 causal effect
filter_regions <- function(regionlist, group_prior, prob_single = 0.8, zdf){
  regionlist2 <- regionlist
  for (b in 1: length(regionlist)){
    for (rn in names(regionlist[[b]])) {
      gid <- regionlist[[b]][[rn]][["gid"]]
      sid <- regionlist[[b]][[rn]][["sid"]]
      gs_type <- zdf$type[match(c(gid,sid), zdf$id)]
      
      #pi_prior <- unname(group_prior[gs_type])
      #P1 <- prod(1-pi_prior) * (1 + sum(pi_prior/(1-pi_prior)))
      
      group_size <- table(gs_type)[names(group_prior)]
      group_size[is.na(group_size)] <- 0
      
      P1 <- prod((1-group_prior)^group_size) * (1 + sum(group_size*(group_prior/(1-group_prior))))
      
      if (P1 < prob_single){
        regionlist2[[b]][[rn]] <- NULL
      }
    }
  }
  regionlist2
}

filter_regions_stable <- function(regionlist, group_prior, prob_single = 0.8){
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

