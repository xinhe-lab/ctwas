#' annotate susie results. this will generate out put txt files
#'
anno_susie <- function(susieres,
                       geneinfo,
                       snpinfo,
                       gidx,
                       sidx,
                       region_tag1,
                       region_tag2,
                       type) {

  anno.gene <- NULL
  if (length(geneinfo) !=0){
    anno.gene <- cbind(geneinfo[gidx,  c("chrom", "id", "p0")],
                       type)
    colnames(anno.gene) <-  c("chrom", "id", "pos", "type")
  }

  anno.SNP <- cbind(snpinfo[sidx, c("chrom", "id", "pos")],
                    "SNP")
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
      #TODO: note this ignores the fact that some variants can belong to multiple CS
    }
  }
  outdf.rn <- cbind(anno, susieres$pip)
  colnames(outdf.rn)[8] <- "susie_pip"
  p <- length(gidx) + length(sidx)
  outdf.rn$mu2 <- colSums(susieres$mu2[ ,
                                        seq(1, p)[1:p!=susieres$null_index], drop = F]) #WARN: not sure for L>1

  outdf.rn
}

