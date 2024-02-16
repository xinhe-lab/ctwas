#' Compute gene-gene and gene-SNP correlations
#'
#' @param z_snp A data frame with columns: "id", "A1", "A2", "z". giving the z scores for
#' snps. "A1" is effect allele. "A2" is the other allele. For harmonized data, A1 and A2 are not required.
#'
#' @param weight a string, pointing to a directory with the FUSION/TWAS format of weights, or a .db file in predictdb format.
#' A vector of multiple sets of weights in PredictDB format can also be specified; genes will have their filename appended
#' to their gene name to ensure IDs are unique.
#'
#' @param region_info a data frame of region definition and associated file names
#'
#' @param weight_format a string, the format of weight, PredictDB or FUSION
#'
#' @return region_cor_info: a data frame containing LD file names for each region
#'
#' @importFrom logging loginfo
#'
compute_cor <- function(z_snp,
                        weight,
                        region_info,
                        weight_format) {

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

  #load all weights info
  wgtall <- lapply(exprvarfs, function(x){load(paste0(strsplit(x, ".exprvar")[[1]], ".exprqc.Rd")); wgtlist})
  wgtlistall <- do.call(c, wgtall)
  names(wgtlistall) <- do.call(c, lapply(wgtall, names))
  rm(wgtall)

  regionlist <- list()
  for (b in 1:length(exprvarfs)){

    if (!is.null(pvarfs)){
      # get snp info (from pvarf file)
      pvarf <- pvarfs[b]
      snpinfo <- ctwas:::read_pvar(pvarf)
    } else {
      # get snp info (from LD R matrix file)
      ld_Rf <- ld_Rfs[b]
      snpinfo <- ctwas:::read_ld_Rvar(ld_Rf)
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


    ld_Rf <- ld_Rfs[b]
    ld_Rinfo <- as.data.frame(data.table::fread(ld_Rf, header = T))

    if (!isTRUE(merge) & nrow(regions) >=2){
    for (rn in 1:(nrow(regions)-1)){
      current <- regionlist[[b]][[as.character(rn)]]
      nextone <- regionlist[[b]][[as.character(rn+1)]]
      gnames <- regionlist[[b]][[as.character(rn)]][["gid"]]
      tmp_region <- regionlist[[b]][[as.character(rn)]]
      if(length(gnames>0)){
        ifreg <- ifelse(regionlist[[b]][[as.character(rn)]][["start"]] < ld_Rinfo[, "stop"] & regionlist[[b]][[as.character(rn)]][["stop"]] >= ld_Rinfo[, "start"], T, F)
        regRDS <- ld_Rinfo[ifreg, "RDS_file"]
        R_snp_anno <- as.data.frame(do.call(rbind, lapply(regRDS, ctwas:::read_ld_Rvar_RDS)))
        for (i in 1:length(gnames)){
          gname <- gnames[i]
          wgt <- wgtlistall[[gname]]
          snpnames <- rownames(wgt)
          ld.idx <- match(snpnames, R_snp_anno$id)
          if(anyNA(ld.idx)){
            thisindex <- !is.na(ld.idx)
            nextindex <- is.na(ld.idx)
            thisr2 <- sum(wgt[thisindex]^2)
            nextr2 <- sum(wgt[nextindex]^2)
            if(thisr2<nextr2){
              #modify weights file - drop weights in other regions
              tmp_wgt <- wgtlistall[[gname]][nextindex]
              if(length(tmp_wgt)==1){
                wgtlistall[[gname]] <- matrix(tmp_wgt,nrow = 1,ncol = 1)
                rownames(wgtlistall[[gname]]) <- snpnames[nextindex]
                colnames(wgtlistall[[gname]]) <- "weight"
              }
              else{
                wgtlistall[[gname]] <- matrix(tmp_wgt,nrow = length(tmp_wgt),ncol = 1)
                rownames(wgtlistall[[gname]]) <- snpnames[nextindex]
                colnames(wgtlistall[[gname]]) <- "weight"
              }
              #add gene to next region
              regionlist[[b]][[as.character(rn+1)]][["gidx"]] <- c(regionlist[[b]][[as.character(rn+1)]][["gidx"]],tmp_region[["gidx"]][which(gnames==gname)])
              regionlist[[b]][[as.character(rn+1)]][["gid"]] <- c(regionlist[[b]][[as.character(rn+1)]][["gid"]],gname)
              #remove gene from this region
              regionlist[[b]][[as.character(rn)]][["gidx"]] <- regionlist[[b]][[as.character(rn)]][["gidx"]][which(gnames!=gname)]
              regionlist[[b]][[as.character(rn)]][["gid"]] <- regionlist[[b]][[as.character(rn)]][["gid"]][!regionlist[[b]][[as.character(rn)]][["gid"]]==gname]
            }
            else{
              #modify weights file - drop weights in other regions
              tmp_wgt <- wgtlistall[[gname]][thisindex]
              if(length(tmp_wgt)==1){
                wgtlistall[[gname]] <- matrix(tmp_wgt,nrow = 1,ncol = 1)
                rownames(wgtlistall[[gname]]) <- snpnames[thisindex]
                colnames(wgtlistall[[gname]]) <- "weight"
              }
              else{
                wgtlistall[[gname]] <- matrix(tmp_wgt,nrow = length(tmp_wgt),ncol = 1)
                rownames(wgtlistall[[gname]]) <- snpnames[thisindex]
                colnames(wgtlistall[[gname]]) <- "weight"
              }
            }
          }
        }
      }
    }
  }
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
          #print("passed calculate correlations")
          #remove genes with NA in R_snp_gene (e.g. boundary spanning with merge=F) from analysis
          gene_not_NA <- which(apply(R_snp_gene, 2, function(x){!any(is.na(x))}))
          gene_NA <- which(apply(R_snp_gene, 2, function(x){any(is.na(x))}))
          outlist_core_region[["gidx"]] <- regionlist[[b]][[rn]][["gidx"]][gene_not_NA]
          outlist_core_region[["gid"]] <- regionlist[[b]][[rn]][["gid"]][gene_not_NA]
          R_snp_gene <- R_snp_gene[,gene_not_NA,drop=F]
          R_gene <- R_gene[gene_not_NA,gene_not_NA,drop=F]
          #print("passed remove NA genes")
          #save R_snp_gene and R_gene and R_SNP
          outlist_core_region[["regRDS"]] <- regRDS
          R_sg_file <- file.path(outputdir, paste0(outname, "_LDR"), paste0("chr", b, "_reg", rn, ".R_snp_gene.RDS"))
          R_g_file <- file.path(outputdir, paste0(outname, "_LDR"), paste0("chr", b, "_reg", rn, ".R_gene.RDS"))
          R_s_file <- file.path(outputdir, paste0(outname, "_LDR"), paste0("chr", b, "_reg", rn, ".R_snp.RDS"))
          R_snp <- R_snp[sidx, sidx, drop = F]
          saveRDS(R_snp, file=R_s_file)
          saveRDS(R_snp_gene, file=R_sg_file)
          saveRDS(R_gene, file=R_g_file)

          outlist_core_region[["R_sg_file"]] <- R_sg_file
          outlist_core_region[["R_g_file"]] <- R_g_file
          outlist_core_region[["R_s_file"]] <- R_s_file
          #print("saved R_snp")
          outlist_core[[length(outlist_core)+1]] <- outlist_core_region
          #print("laste step in dopar")
        }
      }
      #print("finish rn loop")
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

  # d.f. containing LD file names for each region
  region_cor_info <- data.frame(R_sg_file, R_g_file, R_s_file)

  return(region_cor_info)
}
