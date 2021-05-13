require(foreach)

impute_expr_z_LDR <- 
  function (z_snp, ld_pgenf = NULL, weight, method = "lasso", outputdir = getwd(), 
            outname = NULL, logfile = NULL, compress = T, harmonize = T,
            ldmat_file = NULL,
            ld_regions = c("EUR", "ASN", "AFR"),
            ld_regions_custom = NULL,
            merge = T,
            minvar = 1)
  {
    dir.create(outputdir, showWarnings=F)
    if (!is.null(logfile)) {
      addHandler(writeToFile, file = logfile, level = "DEBUG")
    }
    outname <- file.path(outputdir, outname)
    if (is.null(ldmat_file)){
      ld_pvarf <- prep_pvar(ld_pgenf, outputdir = outputdir)
      ld_snpinfo <- read_pvar(ld_pvarf)
    } else {
      dir.create(paste0(outputdir,"/LDR"), showWarnings=F)
      
      ld_snpinfo <- data.table::fread(paste0(ldmat_file,".hbim"))
      ld_snpinfo <- ld_snpinfo[,c("#CHROM", "ID", "POS", "ALT", "REF")]
      colnames(ld_snpinfo) <- c("chrom", "id", "pos", "alt", "ref")
    }
    b <- unique(ld_snpinfo$chrom)
    if (length(b) != 1) {
      stop("Input LD reference not split by chromosome")
    }
    if (isTRUE(harmonize)) {
      logging::loginfo("flipping z scores to match LD reference")
      z_snp <- harmonize_z_ld(z_snp, ld_snpinfo)
      logging::loginfo("will flip weights to match LD reference for each gene")
    }
    exprf <- paste0(outname, "_chr", b, ".expr")
    exprvarf <- paste0(outname, "_chr", b, ".exprvar")
    exprqcf <- paste0(outname, "_chr", b, ".exprqc.Rd")
    logging::loginfo("expression z score inmputation started for chr %s.", b)
    exprlist <- list()
    qclist <- list()
    wgtdir <- dirname(weight)
    wgtposfile <- file.path(wgtdir, paste0(basename(weight), ".pos"))
    wgtpos <- read.table(wgtposfile, header = T, stringsAsFactors = F)
    wgtpos <- transform(wgtpos, ID = ifelse(duplicated(ID) | duplicated(ID, fromLast = TRUE), 
                                            paste(ID, ave(ID, ID, FUN = seq_along), sep = "_ID"), ID))
    logging::loginfo("number of genes with weights provided: %s", nrow(wgtpos))
    wgtpos <- wgtpos[wgtpos$CHR==b,]
    logging::loginfo("number of genes on chromosome: %s", nrow(wgtpos))
    for (i in 1:nrow(wgtpos)) {
      wf <- file.path(wgtdir, wgtpos[i, "WGT"])
      load(wf)
      gname <- wgtpos[i, "ID"]
      if (isTRUE(harmonize)) {
        w <- harmonize_wgt_ld(wgt.matrix, snps, ld_snpinfo)
        wgt.matrix <- w[["wgt"]]
        snps <- w[["snps"]]
      }
      g.method = method
      if (g.method == "best") {
        g.method = names(which.max(cv.performance["rsq",]))
      }
      if (!(g.method %in% names(cv.performance[1, ])))
        next
      wgt.matrix <- wgt.matrix[abs(wgt.matrix[, g.method]) > 0, , drop = F]
      wgt.matrix <- wgt.matrix[complete.cases(wgt.matrix), , drop = F]
      if (nrow(wgt.matrix) == 0)
        next
      snpnames <- Reduce(intersect, list(rownames(wgt.matrix), ld_snpinfo$id, z_snp$id))
      if (length(snpnames) == 0)
        next
      wgt.idx <- match(snpnames, rownames(wgt.matrix))
      wgt <- wgt.matrix[wgt.idx, g.method, drop = F]
      
      exprlist[[gname]][["chrom"]] <- wgtpos[wgtpos$ID == gname, "CHR"]
      exprlist[[gname]][["p0"]] <- min(snps[snps[, "id"] %in% snpnames, "pos"])
      exprlist[[gname]][["p1"]] <- max(snps[snps[, "id"] %in% snpnames, "pos"])
      exprlist[[gname]][["wgt"]] <- wgt
      qclist[[gname]][["n"]] <- nrow(wgt.matrix)
      qclist[[gname]][["nmiss"]] <- nrow(wgt.matrix) - length(snpnames)
      qclist[[gname]][["missrate"]] <- qclist[[gname]][["nmiss"]]/qclist[[gname]][["n"]]
    }
    gnames <- names(exprlist)
    wgtlist <- lapply(exprlist, "[[", "wgt")
    chrom <- unlist(lapply(exprlist, "[[", "chrom"))
    p0 <- unlist(lapply(exprlist, "[[", "p0"))
    p1 <- unlist(lapply(exprlist, "[[", "p1"))
    geneinfo <- data.frame(chrom = chrom, id = gnames, p0 = p0, p1 = p1)
    if (is.null(ldmat_file)){
      logging::loginfo("ldmat_file not detected, using genotypes")
      ld_pgen <- prep_pgen(ld_pgenf, ld_pvarf)
      for (i in 1:length(gnames)){
        gname <- gnames[i]
        wgt <- wgtlist[[gname]]
        snpnames <- rownames(wgt)
        ld.idx <- match(snpnames, ld_snpinfo$id)
        z.idx <- match(snpnames, z_snp$id)
        X.g <- read_pgen(ld_pgen, variantidx = ld.idx)
        X.g <- scale(X.g)
        gexpr <- X.g %*% wgt
        if (abs(max(gexpr) - min(gexpr)) < 1e-08) #this check could cause errors later because the list has already been populated; consider setting "gexpr" and "z.g" to NA to avoid missing fields?
          next
        z.s <- as.matrix(z_snp[z.idx, "z"])
        var.s <- sqrt(apply(X.g, 2, var))
        Gamma.g <- cov(X.g)
        z.g <- (t(wgt) * var.s) %*% z.s/sqrt(t(wgt) %*% Gamma.g %*% wgt)
        exprlist[[gname]][["expr"]] <- gexpr
        exprlist[[gname]][["z.g"]] <- z.g
      }
    } else {
      if (is.null(ld_regions_custom)) {
        ld_regions <- match.arg(ld_regions)
        regionfile <- system.file("extdata", "ldetect", paste0(ld_regions, ".bed"), package = "ctwas")
      } else {
        regionfile <- ld_regions_custom
      }
      logging::loginfo("LD region file: %s", regionfile)
      regions <- read.table(regionfile, header = T, stringsAsFactors = F)
      if (is.character(regions$chr)) {
        regions$chr <- readr::parse_number(regions$chr)
      }
      regions <- regions[regions$chr==b,]
      regions <- regions[order(regions$start), ]
      ##########
      #now following ctwas:::index_regions with modification, removed info about SNPs
      #possible edge case in ctwas_rss - (length(gidx) < minvar) vs (length(gidx) + length(sidx) < minvar)
      regionlist <- list()
      for (rn in 1:nrow(regions)) {
        p0 <- regions[rn, "start"]
        p1 <- regions[rn, "stop"]
        if (isTRUE(merge)) {
          gidx <- which(geneinfo$chrom == b & geneinfo$p1 > p0 & geneinfo$p0 < p1)
        } else {
          gidx <- which(geneinfo$chrom == b & geneinfo$p0 > p0 & geneinfo$p0 < p1)
        }
        gid <- geneinfo[gidx, "id"]
        if (length(gidx) < minvar) {
          next
        }
        regionlist[[as.character(rn)]] <- list(gidx = gidx, gid = gid, start = p0, stop = p1)
      }
      logging::loginfo("No. regions with at least one gene for chr%s: %s", b, length(regionlist))
      if (isTRUE(merge) & nrow(regions) >= 2) {
        logging::loginfo("Merge regions for chr%s", b)
        for (rn in 2:nrow(regions)) {
          current <- regionlist[[as.character(rn)]]
          previous <- regionlist[[as.character(rn - 1)]]
          if (length(intersect(current[["gid"]], previous[["gid"]])) > 0) {
            gidx <- unique(c(previous[["gidx"]], current[["gidx"]]))
            gid <- geneinfo[gidx, "id"]
            regionlist[[as.character(rn)]] <- list(gidx = gidx,
                                                   gid = gid,
                                                   start = regionlist[[as.character(rn - 1)]]$start,
                                                   stop = regions[rn, "stop"])
            regionlist[[as.character(rn - 1)]] <- NULL
          }
        }
      }
      logging::loginfo("No. regions with at least one gene for chr%s after merging: %s", b, length(regionlist))
      ##########
      for (rn in 1:length(regionlist)){
        reg_current <- regions[regions$stop > regionlist[[rn]][["start"]] & regions$start < regionlist[[rn]][["stop"]],]
        R_snp <- lapply(1:nrow(reg_current), function(x){
          readRDS(paste0(ldmat_file, ".R_snp.", paste(reg_current[x,c("start","stop")], collapse="_"), ".RDS"))})
        R_snp_names <- unlist(sapply(R_snp, rownames))
        R_snp <- as.matrix(Matrix::bdiag(R_snp))
        rownames(R_snp) <- R_snp_names
        colnames(R_snp) <- R_snp_names
        gid <- regionlist[[rn]][["gid"]]
        R_snp_gene <- matrix(NA,nrow(R_snp),length(gid))
        rownames(R_snp_gene) <- rownames(R_snp)
        colnames(R_snp_gene) <- gid
        for (gname in gid){
          wgt <- exprlist[[gname]][["wgt"]]
          snpnames <- rownames(wgt)
          ld.idx <- match(snpnames, rownames(R_snp))
          zdf.idx <- match(snpnames, z_snp$id)
          R.s <- R_snp[ld.idx,ld.idx]
          z.s <-  as.matrix(z_snp[zdf.idx, "z"])
          z.g <- crossprod(wgt,z.s)/sqrt(crossprod(wgt,R.s)%*%wgt)
          R_snp_gene[,gname] <- sapply(1:nrow(R_snp), 
                                       function(x){crossprod(wgt,R_snp[ld.idx,x])/sqrt(crossprod(wgt,R.s)%*%wgt*R_snp[x,x])})
          exprlist[[gname]][["ld.idx"]] <- ld.idx
          exprlist[[gname]][["z.g"]] <- z.g
        }
        wgtlist <- lapply(exprlist[gid], "[[", "wgt")
        ld.idx.list <- lapply(exprlist[gid], "[[", "ld.idx")
        R_gene <- diag(length(gid))
        rownames(R_gene) <- gid
        colnames(R_gene) <- gid
        if (length(wgtlist)>1){
          gene_pairs <- combn(length(wgtlist), 2)
          gene_corrs <- apply(gene_pairs, 2, function(x){t(wgtlist[[x[1]]])%*%R_snp[ld.idx.list[[x[1]]], ld.idx.list[[x[2]]]]%*%wgtlist[[x[2]]]/(
            sqrt(t(wgtlist[[x[1]]])%*%R_snp[ld.idx.list[[x[1]]], ld.idx.list[[x[1]]]]%*%wgtlist[[x[1]]]) *
              sqrt(t(wgtlist[[x[2]]])%*%R_snp[ld.idx.list[[x[2]]], ld.idx.list[[x[2]]]]%*%wgtlist[[x[2]]]))})
          R_gene[t(gene_pairs)] <- gene_corrs
          R_gene[t(gene_pairs[c(2,1),])] <- gene_corrs
        }
        file_stem <- rev(unlist(strsplit(ldmat_file,split="/")))[1]
        saveRDS(R_snp_gene, file=paste0(out_dir, "/LDR/", file_stem, ".R_snp_gene.", regionlist[[rn]][["start"]], "_", regionlist[[rn]][["stop"]], ".RDS"))
        saveRDS(R_gene, file=paste0(out_dir, "/LDR/", file_stem, ".R_gene.", regionlist[[rn]][["start"]], "_", regionlist[[rn]][["stop"]], ".RDS"))
      }
    }
    z.g <- unlist(lapply(exprlist, "[[", "z.g"))
    expr <- do.call(cbind, lapply(exprlist, "[[", "expr"))
    if (length(exprlist) == 0) {
      expr <- data.table::data.table(NULL)
      geneinfo <- data.table::data.table(NULL)
    }
    logging::loginfo("Number of genes with imputed expression: %s for chr %s", length(z.g), b)
    if (is.null(ldmat_file)){
      data.table::fwrite(expr, file = exprf, row.names = F, col.names = F, sep = "\t", quote = F)
      if (isTRUE(compress)) {
        system(paste0("gzip -f ", exprf))
        exprf <- paste0(exprf, ".gz")
      }
    }
    data.table::fwrite(geneinfo, file = exprvarf, sep = "\t", quote = F)
    save(wgtlist, qclist, file = exprqcf)
    logging::loginfo("expression inmputation done for chr %s.", b)
    z_gene <- data.frame(id = gnames, z = z.g)
    return(list(z_gene = z_gene, ld_exprf = exprf))
  }

index_regions_LDR <- function (ldmat_files,
                               exprvarfs, regionfile, select = NULL, thin = 1, 
                               maxSNP = Inf, minvar = 1, merge = T) 
{
  reg <- read.table(regionfile, header = T, stringsAsFactors = F)
  if (is.character(reg$chr)) {
    reg$chr <- readr::parse_number(reg$chr)
  }
  reg <- reg[order(reg$chr, reg$start), ]
  logging::loginfo("No. LD regions: %s", nrow(reg))
  if (thin > 1 | thin <= 0) {
    stop("thin needs to be in (0,1]")
  }
  if (is.null(dim(select))) {
    selectid <- select
  }
  else {
    selectid <- select$id
  }
  regionlist <- list()
  for (b in 1:length(ldmat_files)) {
    ldmat_file <- ldmat_files[b]
    snpinfo <- data.table::fread(paste0(ldmat_file,".hbim"))
    snpinfo <- snpinfo[,c("#CHROM", "ID", "POS", "ALT", "REF")]
    colnames(snpinfo) <- c("chrom", "id", "pos", "alt", "ref")
    if (unique(snpinfo$chrom) != b) {
      stop("Input genotype file not split by chromosome or not in correct order")
    }
    snpinfo$keep <- 1
    if (!is.null(selectid)) {
      snpinfo[!(snpinfo$id %in% selectid), "keep"] <- 0
    }
    snpinfo$thin_tag <- 0
    nkept <- round(nrow(snpinfo) * thin)
    set.seed(99)
    snpinfo$thin_tag[sample(1:nrow(snpinfo), nkept)] <- 1
    exprvarf <- exprvarfs[b]
    geneinfo <- read_exprvar(exprvarf)
    if (nrow(geneinfo) != 0) {
      if (unique(geneinfo$chrom) != b) {
        stop("Imputed expression not by chromosome or not in correct order")
      }
      geneinfo$keep <- 1
      if (!is.null(selectid)) {
        geneinfo[!(geneinfo$id %in% selectid), "keep"] <- 0
      }
    }
    regions <- reg[reg$chr == b, ]
    regionlist[[b]] <- list()
    for (rn in 1:nrow(regions)) {
      p0 <- regions[rn, "start"]
      p1 <- regions[rn, "stop"]
      if (isTRUE(merge)) {
        gidx <- which(geneinfo$chrom == b & geneinfo$p1 > 
                        p0 & geneinfo$p0 < p1 & geneinfo$keep == 1)
      } else {
        gidx <- which(geneinfo$chrom == b & geneinfo$p0 > 
                        p0 & geneinfo$p0 < p1 & geneinfo$keep == 1)
      }
      sidx <- which(snpinfo$chrom == b & snpinfo$pos > 
                      p0 & snpinfo$pos < p1 & snpinfo$keep == 1 & snpinfo$thin_tag == 1)
      gid <- geneinfo[gidx, "id"]
      sid <- snpinfo[sidx, "id"]
      if (length(gidx) + length(sidx) < minvar) {
        next
      }
      regionlist[[b]][[as.character(rn)]] <- list(gidx = gidx, 
                                                  gid = gid, sidx = sidx, sid = sid, start = p0, 
                                                  stop = p1)
    }
    logging::loginfo("No. regions with at least one SNP/gene for chr%s: %s", b, length(regionlist[[b]]))
    if (isTRUE(merge) & nrow(regions) >= 2) {
      logging::loginfo("Merge regions for chr%s", b)
      for (rn in 2:nrow(regions)) {
        current <- regionlist[[b]][[as.character(rn)]]
        previous <- regionlist[[b]][[as.character(rn - 1)]]
        #if (length(intersect(current[["gid"]], previous[["gid"]])) > 0) {
        if (length(intersect(current[["gid"]]$id, previous[["gid"]]$id)) > 0) { #the previous line was not behaving correctly, something about data types when run out of the package
          gidx <- unique(c(previous[["gidx"]], current[["gidx"]]))
          sidx <- unique(c(previous[["sidx"]], current[["sidx"]]))
          gid <- geneinfo[gidx, "id"]
          sid <- snpinfo[sidx, "id"]
          regionlist[[b]][[as.character(rn)]] <- list(gidx = gidx, 
                                                      gid = gid, 
                                                      sidx = sidx, 
                                                      sid = sid, 
                                                      start = regionlist[[b]][[as.character(rn - 1)]]$start, 
                                                      stop = regions[rn, "stop"])
          regionlist[[b]][[as.character(rn - 1)]] <- NULL
        }
      }
    }
    logging::loginfo("No. regions with at least one SNP/gene for chr%s after merging: %s", 
                     b, length(regionlist[[b]]))
  }
  logging::loginfo("Trim regions with SNPs more than %s", maxSNP)
  if ("z" %in% colnames(select)) {
    for (b in 1:length(regionlist)) {
      for (rn in names(regionlist[[b]])) {
        if (length(regionlist[[b]][[rn]][["sid"]]) > 
            maxSNP) {
          idx <- match(regionlist[[b]][[rn]][["sid"]], select[, "id"])
          z.abs <- abs(select[idx, "z"])
          ifkeep <- rank(-z.abs) <= maxSNP
          regionlist[[b]][[rn]][["sidx"]] <- regionlist[[b]][[rn]][["sidx"]][ifkeep]
          regionlist[[b]][[rn]][["sid"]] <- regionlist[[b]][[rn]][["sid"]][ifkeep]
        }
      }
    }
  }
  else {
    for (b in 1:length(regionlist)) {
      for (rn in names(regionlist[[b]])) {
        if (length(regionlist[[b]][[rn]][["sid"]]) > 
            maxSNP) {
          n.ori <- length(regionlist[[b]][[rn]][["sid"]])
          ifkeep <- rep(F, n.ori)
          set.seed <- 99
          ifkeep[sample.int(n.ori, size = maxSNP)] <- T
          regionlist[[b]][[rn]][["sidx"]] <- regionlist[[b]][[rn]][["sidx"]][ifkeep]
          regionlist[[b]][[rn]][["sid"]] <- regionlist[[b]][[rn]][["sid"]][ifkeep]
        }
      }
    }
  }
  regionlist
}

anno_susie_LDR <- 
  function (susieres, exprvarf, 
            ldmat_file,
            gidx, sidx, region_tag1, 
            region_tag2) 
  {
    geneinfo <- read_exprvar(exprvarf)
    anno.gene <- NULL
    if (length(geneinfo) != 0) {
      anno.gene <- cbind(geneinfo[gidx, c("chrom", "id", "p0")], rep("gene", length(gidx)))
      colnames(anno.gene) <- c("chrom", "id", "pos", "type")
    }
    ##########
    snpinfo <- data.table::fread(paste0(ldmat_file,".hbim"))
    snpinfo <- snpinfo[,c("#CHROM", "ID", "POS", "ALT", "REF")]
    colnames(snpinfo) <- c("chrom", "id", "pos", "alt", "ref")
    ##########
    anno.SNP <- cbind(snpinfo[sidx, c("chrom", "id", "pos")], 
                      rep("SNP", length(sidx)))
    colnames(anno.SNP) <- c("chrom", "id", "pos", "type")
    anno <- rbind(anno.gene, anno.SNP)
    anno <- as.data.frame(anno)
    anno$region_tag1 <- region_tag1
    anno$region_tag2 <- region_tag2
    anno$cs_index <- 0
    if (!is.null(susieres$sets$cs)) {
      for (cs_i in susieres$sets$cs_index) {
        X.idx <- susieres$sets$cs[[paste0("L", cs_i)]]
        X.idx <- X.idx[X.idx != susieres$null_index]
        anno$cs_index[X.idx] <- cs_i
      }
    }
    outdf.rn <- cbind(anno, susieres$pip)
    colnames(outdf.rn)[8] <- "susie_pip"
    p <- length(gidx) + length(sidx)
    outdf.rn$mu2 <- colSums(susieres$mu2[, seq(1, p)[1:p != susieres$null_index], drop = F])
    outdf.rn
  }

susieI_rss_LDR <- function (zdf, ld_pgenfs=NULL, ld_exprfs, regionlist, niter = 20,
                            L = 1, z_ld_weight = 0, group_prior = NULL, group_prior_var = NULL,
                            estimate_group_prior = T, estimate_group_prior_var = T, use_null_weight = T,
                            coverage = 0.95, ncore = 1, outputdir = getwd(), outname = NULL,
                            ldmat_files=NULL,
                            ld_regions = c("EUR", "ASN", "AFR"),
                            ld_regions_custom=NULL) {
  outname <- file.path(outputdir, outname)
  ld_exprvarfs <- sapply(ld_exprfs, prep_exprvar)
  K <- 2
  group_prior_rec <- matrix(, nrow = K, ncol = niter)
  group_prior_var_rec <- matrix(, nrow = K, ncol = niter)
  prior.gene <- group_prior[1]
  prior.SNP <- group_prior[2]
  V.gene <- group_prior_var[1]
  V.SNP <- group_prior_var[2]
  for (iter in 1:niter) {
    logging::loginfo("run iteration %s", iter)
    snp.rpiplist <- list()
    gene.rpiplist <- list()
    cl <- parallel::makeCluster(ncore, outfile = "")
    doParallel::registerDoParallel(cl)
    corelist <- ctwas:::region2core(regionlist, ncore)
    
    outdf <- foreach(core = 1:length(corelist), .combine = "rbind", .packages = "ctwas") %dopar% {
      source("ctwas_LDR.R") #this line was necessary to load anno_susie_LDR on each core
      outdf.core.list <- list()
      regs <- corelist[[core]]
      for (reg in 1:nrow(regs)) {
        logging::loginfo(reg)
        b <- regs[reg, "b"]
        rn <- regs[reg, "rn"]
        if (is.null(ldmat_files)){
          ld_pgen <- prep_pgen(pgenf = ld_pgenfs[b], ld_pvarfs[b])
        }
        gidx <- regionlist[[b]][[rn]][["gidx"]]
        sidx <- regionlist[[b]][[rn]][["sidx"]]
        #gid <- regionlist[[b]][[rn]][["gid"]]
        #sid <- regionlist[[b]][[rn]][["sid"]]
        gid <- regionlist[[b]][[rn]][["gid"]][["id"]] #this was related to the same issue with data types when not run from the package
        sid <- regionlist[[b]][[rn]][["sid"]][["id"]] #this was related to the same issue with data types when not run from the package
        p <- length(gidx) + length(sidx)
        if (is.null(prior.gene) | is.null(prior.SNP)) {
          prior <- c(rep(1/p, length(gidx)), rep(1/p, length(sidx)))
        } else {
          prior <- c(rep(prior.gene, length(gidx)), rep(prior.SNP, length(sidx)))
        }
        if (is.null(V.gene) | is.null(V.SNP)) {
          V <- matrix(rep(50, L * p), nrow = L)
        } else {
          V <- c(rep(V.gene, length(gidx)), rep(V.SNP, 
                                                length(sidx)))
          V <- matrix(rep(V, each = L), nrow = L)
        }
        if (isTRUE(use_null_weight)) {
          nw <- max(0, 1 - sum(prior))
          prior <- prior/(1 - nw)
        } else {
          nw <- NULL
        }
        z.g <- zdf[match(gid, zdf$id), ][["z"]]
        z.s <- zdf[match(sid, zdf$id), ][["z"]]
        z <- c(z.g, z.s)
        if (is.null(ldmat_files)){
          X.g <- read_expr(ld_exprfs[b], variantidx = gidx)
          X.s <- read_pgen(ld_pgen, variantidx = sidx)
          X <- cbind(X.g, X.s)
          R <- Rfast::cora(X)
        } else {
          if (is.null(ld_regions_custom)) {
            #ld_regions <- match.arg(ld_regions) #try uncommenting this line, it gave me an error
            regionfile <- system.file("extdata", "ldetect", paste0(ld_regions, ".bed"), package = "ctwas")
          } else {
            regionfile <- ld_regions_custom
          }
          #logging::loginfo("LD region file: %s", regionfile)
          regions <- read.table(regionfile, header = T, stringsAsFactors = F)
          if (is.character(regions$chr)) {
            regions$chr <- readr::parse_number(regions$chr)
          }
          regions <- regions[regions$chr==b,]
          regions <- regions[order(regions$start), ]
          p0 <- regionlist[[b]][[rn]][["start"]]
          p1 <- regionlist[[b]][[rn]][["stop"]]
          reg_current <- regions[regions$stop > p0 & regions$start < p1,]
          R_snp <- lapply(1:nrow(reg_current), function(x){
            readRDS(paste0(ldmat_files[b], ".R_snp.", paste(reg_current[x,c("start","stop")], collapse="_"), ".RDS"))})
          R_snp_names <- unlist(sapply(R_snp, rownames))
          R_snp <- as.matrix(Matrix::bdiag(R_snp))
          rownames(R_snp) <- R_snp_names
          colnames(R_snp) <- R_snp_names
          R_snp <- R_snp[sid,sid,drop=F]
          file_stem <- rev(unlist(strsplit(ldmat_files[b],split="/")))[1]
          R_snp_gene_f <- paste0(outputdir,"/LDR/", file_stem, ".R_snp_gene.", p0, "_", p1, ".RDS")
          if (file.exists(R_snp_gene_f)){
            R_snp_gene <- readRDS(paste0(outputdir,"/LDR/", file_stem, ".R_snp_gene.", p0, "_", p1, ".RDS"))
            R_snp_gene <- R_snp_gene[sid,gid,drop=F]
            R_gene <- readRDS(paste0(outputdir,"/LDR/", file_stem, ".R_gene.", p0, "_", p1, ".RDS"))
            R_gene <- R_gene[gid,gid,drop=F]
            R <- rbind(cbind(R_gene, t(R_snp_gene)), 
                       cbind(R_snp_gene, R_snp))
          } else {
            R <- R_snp
          }
        }
        susieres <- susie_rss(z, R, z_ld_weight = z_ld_weight, 
                              L = L, prior_weights = prior, null_weight = nw, 
                              prior_variance = V, estimate_prior_variance = F, 
                              coverage = coverage)
        if (is.null(ldmat_files)){
          outdf.reg <- ctwas:::anno_susie(susieres, ld_exprvarfs[b], ld_pvarfs[b], gidx, sidx, b, rn)
        } else {
          outdf.reg <- anno_susie_LDR(susieres, ld_exprvarfs[b], ldmat_files[b], gidx, sidx, b, rn)
        }
        outdf.core.list[[reg]] <- outdf.reg
      }
      outdf.core <- do.call(rbind, outdf.core.list)
      outdf.core
    }
    if (isTRUE(estimate_group_prior)) {
      prior.SNP <- mean(outdf[outdf[, "type"] == "SNP", "susie_pip"])
      prior.gene <- mean(outdf[outdf[, "type"] == "gene", "susie_pip"])
      group_prior_rec[, iter] <- c(prior.gene, prior.SNP)
    }
    logging::loginfo("After iteration %s, gene prior %s:, SNP prior:%s", iter, prior.gene, prior.SNP)
    if (isTRUE(estimate_group_prior_var)) {
      outdf.g <- outdf[outdf[, "type"] == "gene", ]
      outdf.s <- outdf[outdf[, "type"] == "SNP", ]
      V.gene <- sum(outdf.g$susie_pip * outdf.g$mu2)/sum(outdf.g$susie_pip)
      V.SNP <- sum(outdf.s$susie_pip * outdf.s$mu2)/sum(outdf.s$susie_pip)
      group_prior_var_rec[, iter] <- c(V.gene, V.SNP)
    }
    save(group_prior_rec, group_prior_var_rec, file = paste0(outname, ".susieIrssres.Rd"))
    data.table::fwrite(outdf, file = paste0(outname, ".susieIrss.txt"), sep = "\t", quote = F)
    parallel::stopCluster(cl)
  }
  list(group_prior = c(prior.gene, prior.SNP), group_prior_var = c(V.gene, V.SNP))
}

ctwas_rss_LDR <- 
  function (z_snp, z_gene, ld_pgenfs=NULL, ld_exprfs, ld_regions = c("EUR", "ASN", "AFR"), ld_regions_custom = NULL, thin = 1, prob_single = 0.8, 
            rerun_gene_PIP = 0.8, niter1 = 3, niter2 = 30, L = 5, group_prior = NULL, 
            group_prior_var = NULL, estimate_group_prior = T, estimate_group_prior_var = T, 
            use_null_weight = T, coverage = 0.95, stardardize = T, harmonize = T, 
            max_snp_region = Inf, ncore = 1, ncore.rerun = 1, outputdir = getwd(), 
            outname = NULL, logfile = NULL,
            ldmat_files = NULL) 
  {
    if (!is.null(logfile)) {
      addHandler(writeToFile, file = logfile, level = "DEBUG")
    }
    logging::loginfo("ctwas started ... ")
    if (length(ld_pgenfs) != 22 & is.null(ldmat_files)) {
      stop("Not all pgen files for 22 chromosomes are provided.")
    }
    if (length(ld_exprfs) != 22) {
      stop("Not all imputed expression files for 22 chromosomes are provided.")
    }
    if (is.null(ld_regions_custom)) {
      ld_regions <- match.arg(ld_regions)
      regionfile <- system.file("extdata", "ldetect", paste0(ld_regions, ".bed"), package = "ctwas")
    } else {
      regionfile <- ld_regions_custom
    }
    logging::loginfo("LD region file: %s", regionfile)
    if (is.null(ldmat_files)){
      ld_pvarfs <- sapply(ld_pgenfs, prep_pvar, outputdir = outputdir)
    }
    ld_exprvarfs <- sapply(ld_exprfs, prep_exprvar)
    if (isTRUE(harmonize)) {
      logging::loginfo("flipping z scores to match LD reference")
      if (is.null(ldmat_files)){
        for (ld_pvarf in ld_pvarfs) {
          ld_snpinfo <- read_pvar(ld_pvarf)
          z_snp <- harmonize_z_ld(z_snp, ld_snpinfo)
        }
      } else {
        for (ldmat_file in ldmat_files){
          ld_snpinfo <- data.table::fread(paste0(ldmat_file,".hbim"))
          ld_snpinfo <- ld_snpinfo[,c("#CHROM", "ID", "POS", "ALT", "REF")]
          colnames(ld_snpinfo) <- c("chrom", "id", "pos", "alt", "ref")
          z_snp <- harmonize_z_ld(z_snp, ld_snpinfo)
        }
      }
    }
    zdf <- rbind(z_snp[, c("id", "z")], z_gene[, c("id", "z")])
    rm(z_snp)
    if (thin <= 0 | thin > 1) {
      stop("thin value needs to be in (0,1]")
    }
    if (is.null(ldmat_files)){
      regionlist <- ctwas:::index_regions(ld_pvarfs, ld_exprvarfs, regionfile, select = zdf$id, thin = thin, minvar = 2)
    } else {
      regionlist <- index_regions_LDR(ldmat_files, ld_exprvarfs, regionfile, select = zdf$id, thin = thin, minvar = 2)
    }
    
    regs <- do.call(rbind, lapply(1:22, function(x) cbind(x, unlist(lapply(regionlist[[x]], "[[", "start")), unlist(lapply(regionlist[[x]], "[[", "stop")))))
    write.table(regs, file = paste0(outputdir, "/", outname, ".regions.txt"), row.names = F, col.names = T, sep = "\t", quote = F)
    if (isTRUE(estimate_group_prior) | isTRUE(estimate_group_prior_var)) {
      logging::loginfo("Run susie iteratively, getting rough estimate ...")
      if (!is.null(group_prior)) {
        group_prior[2] <- group_prior[2]/thin
      }
      pars <- susieI_rss_LDR(zdf = zdf, ld_pgenfs = ld_pgenfs, 
                             ld_exprfs = ld_exprfs, regionlist = regionlist, niter = niter1, 
                             L = 1, z_ld_weight = 0, group_prior = group_prior, 
                             group_prior_var = group_prior_var, estimate_group_prior = estimate_group_prior, 
                             estimate_group_prior_var = estimate_group_prior_var, 
                             use_null_weight = use_null_weight, coverage = coverage, 
                             ncore = ncore, outputdir = outputdir, outname = paste0(outname,".s1"),
                             ldmat_files = ldmat_files,
                             ld_regions = ld_regions)
      group_prior <- pars[["group_prior"]]
      group_prior_var <- pars[["group_prior_var"]]
      regionlist2 <- ctwas:::filter_regions(regionlist, group_prior, prob_single = prob_single)
      logging::loginfo("Blocks are filtered: %s blocks left", sum(unlist(lapply(regionlist2, length))))
      logging::loginfo("Run susie iteratively, getting accurate estimate ...")
      pars <- susieI_rss_LDR(zdf = zdf, ld_pgenfs = ld_pgenfs, 
                             ld_exprfs = ld_exprfs, regionlist = regionlist2, 
                             niter = niter2, L = 1, z_ld_weight = 0, group_prior = group_prior, 
                             group_prior_var = group_prior_var, estimate_group_prior = estimate_group_prior, 
                             estimate_group_prior_var = estimate_group_prior_var, 
                             use_null_weight = use_null_weight, coverage = coverage, ncore = ncore, outputdir = outputdir, outname = paste0(outname, ".s2"),
                             ldmat_files = ldmat_files,
                             ld_regions = ld_regions)
      group_prior <- pars[["group_prior"]]
      group_prior_var <- pars[["group_prior_var"]]
    }
    logging::loginfo("Run susie for all regions.")
    pars <- susieI_rss_LDR(zdf = zdf, ld_pgenfs = ld_pgenfs, ld_exprfs = ld_exprfs, 
                           regionlist = regionlist, niter = 1, L = L, z_ld_weight = 0, 
                           group_prior = group_prior, group_prior_var = group_prior_var, 
                           estimate_group_prior = estimate_group_prior, estimate_group_prior_var = estimate_group_prior_var, 
                           use_null_weight = use_null_weight, coverage = coverage, 
                           ncore = ncore, outputdir = outputdir, outname = paste0(outname, ".temp"),
                           ldmat_files = ldmat_files,
                           ld_regions = ld_regions)
    group_prior[2] <- group_prior[2] * thin
    if (thin == 1) {
      file.rename(paste0(file.path(outputdir, outname), ".temp.susieIrss.txt"), 
                  paste0(file.path(outputdir, outname), ".susieIrss.txt"))
    } else {
      regionlist <- index_regions_LDR(ldmat_files, ld_exprvarfs, 
                                      regionfile, select = zdf, thin = 1, maxSNP = max_snp_region, 
                                      minvar = 2)
      res <- data.table::fread(paste0(file.path(outputdir, 
                                                outname), ".temp.susieIrss.txt"))
      res.keep <- NULL
      for (b in 1:length(regionlist)) {
        for (rn in names(regionlist[[b]])) {
          gene_PIP <- max(res[res$type == "gene" & res$region_tag1 == 
                                b & res$region_tag2 == rn, ]$susie_pip, 0)
          if (gene_PIP < rerun_gene_PIP) {
            regionlist[[b]][[rn]] <- NULL
            res.keep <- rbind(res.keep, res[res$region_tag1 == b & res$region_tag2 == rn, ])
          }
        }
      }
      nreg <- sum(unlist(lapply(regionlist, length)))
      loginfo("Number of regions that contains strong gene signals: %s", nreg)
      if (nreg == 0) {
        file.rename(paste0(file.path(outputdir, outname), 
                           ".temp.susieIrss.txt"), paste0(file.path(outputdir, 
                                                                    outname), ".susieIrss.txt"))
      } else {
        loginfo("Rerun susie for regions with strong gene signals using full SNPs.")
        pars <- susieI_rss_LDR(zdf = zdf, ld_pgenfs = ld_pgenfs, 
                               ld_exprfs = ld_exprfs, regionlist = regionlist, 
                               niter = 1, L = L, z_ld_weight = 0, group_prior = group_prior, 
                               group_prior_var = group_prior_var, estimate_group_prior = estimate_group_prior, 
                               estimate_group_prior_var = estimate_group_prior_var, 
                               use_null_weight = use_null_weight, coverage = coverage, 
                               ncore = ncore.rerun, outputdir = outputdir, outname = paste0(outname, ".s3"),
                               ldmat_files = ldmat_files,
                               ld_regions = ld_regions)
        res.rerun <- data.table::fread(paste0(file.path(outputdir, outname), ".s3.susieIrss.txt"))
        res <- rbind(res.keep, res.rerun)
        data.table::fwrite(res, file = paste0(file.path(outputdir, outname), ".susieIrss.txt"), sep = "\t", quote = F)
        file.remove((paste0(file.path(outputdir, outname), ".temp.susieIrss.txt")))
      }
    }
    list(group_prior = group_prior, group_prior_var = group_prior_var)
  }
