#' annotate susie results. this will generate out put txt files
#'
anno_susie <- function(susieres,
                       exprvarf,
                       pvarf,
                       gidx,
                       sidx,
                       region_tag1,
                       region_tag2) {

  geneinfo <- read_exprvar(exprvarf)

  anno.gene <- cbind(geneinfo[gidx,  c("chrom", "id", "p0")],
                     rep("gene", length(gidx)))
  colnames(anno.gene) <-  c("chrom", "id", "pos", "type")

  snpinfo <- read_pvar(pvarf)

  anno.SNP <- cbind(snpinfo[sidx, c("chrom", "id", "pos")],
                    rep("SNP", length(sidx)))
  colnames(anno.SNP) <-  c("chrom", "id", "pos", "type")

  anno <- rbind(anno.gene, anno.SNP)

  anno <- as.data.frame(anno)

  anno$region_tag1 <- region_tag1
  anno$region_tag2 <- region_tag2

  anno$cs_index <- 0
  if (!is.null(susieres$sets$cs)){
    for (cs_i in susieres$sets$cs_index){
      X.idx <- susieres$sets$cs[[paste0("L", cs_i)]]
      X.idx <- X.idx[X.idx != susieres$null_index] # susie_rss' bug
      anno$cs_index[X.idx] <- cs_i
      #TODO: note this ignore the fact that some variant can belong to multiple CS
    }
  }
  outdf.rn <- cbind(anno, susieres$pip)
  colnames(outdf.rn)[8] <- "susie_pip"
  outdf.rn$mu2 <- colSums(susieres$mu2[ ,
                      seq(1, ncol(X))[1:ncol(X)!=susieres$null_index], drop = F])
  #WARN: not sure for L>1
  outdf.rn
}

