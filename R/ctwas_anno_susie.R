# annotate susie results
anno_susie <- function(susie_res,
                       gene_info,
                       snp_info,
                       gidx,
                       sidx,
                       region_tag,
                       type,
                       QTLtype) {

  anno.gene <- NULL
  if (length(gene_info) !=0){
    anno.gene <- cbind(gene_info[gidx,  c("chrom", "id", "p0")], type, QTLtype)
    colnames(anno.gene) <-  c("chrom", "id", "pos", "type", "QTLtype")
  }

  anno.SNP <- cbind(snp_info[sidx, c("chrom", "id", "pos")], type = "SNP", QTLtype = "SNP")
  colnames(anno.SNP) <-  c("chrom", "id", "pos", "type", "QTLtype")

  res <- as.data.frame(rbind(anno.gene, anno.SNP))

  res$region_tag <- region_tag

  res$cs_index <- 0
  if (!is.null(susie_res$sets$cs)){
    for (cs_i in susie_res$sets$cs_index){
      X.idx <- susie_res$sets$cs[[paste0("L", cs_i)]]
      X.idx <- X.idx[X.idx != susie_res$null_index] # susie_rss' bug
      res$cs_index[X.idx] <- cs_i
      #TODO: note this ignores the fact that some variants can belong to multiple CS
    }
  }

  res$susie_pip <- susie_res$pip

  p <- length(gidx) + length(sidx)
  res$mu2 <- colSums(susie_res$mu2[ ,
                                    seq(1, p)[1:p!=susie_res$null_index], drop = F]) #WARN: not sure for L>1

  return(res)
}

