# identify cross-boundary genes, adjust regionlist and update weigh_list
adjust_boundary_genes <- function(regioninfo, weights, regionlist){
  boundary_genes <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(boundary_genes) <- c("gene","chrom","region1","region2")

  for (i in 1:(nrow(regioninfo)-1)){
    region_tag_current <- regioninfo$region_tag[i]
    region_tag_next <- regioninfo$region_tag[i+1]
    # current <- regionlist[[region_tag_current]]
    # nextone <- regionlist[[region_tag_next]]
    gnames <- regionlist[[region_tag_current]][["gid"]]
    # tmp_region <- regionlist[[region_tag_current]]
    if(length(gnames>0)){
      snpinfo_file <- regioninfo[regioninfo$region_tag == region_tag_current, "SNP_info"]
      ld_snpinfo <- read_LD_SNP_files(snpinfo_file)
      for (i in 1:length(gnames)){
        gname <- gnames[i]
        wgt <- weights[[gname]][["wgt"]]
        snpnames <- rownames(wgt)
        ld.idx <- match(snpnames, ld_snpinfo$id)
        if(anyNA(ld.idx)){ # QTLs are across boundary
          boundary_genes <- rbind(boundary_genes,
                                  data.frame("gene"=gname,
                                             "region1"=region_tag_current,
                                             "region2"=region_tag_next))
          thisindex <- !is.na(ld.idx)
          nextindex <- is.na(ld.idx)
          thisr2 <- sum(wgt[thisindex]^2)
          nextr2 <- sum(wgt[nextindex]^2)
          if(thisr2<nextr2){
            #modify weights file - drop weights in other regions
            tmp_wgt <- weights[[gname]][["wgt"]][nextindex]
            if(length(tmp_wgt)==1){
              weights[[gname]][["wgt"]] <- matrix(tmp_wgt,nrow = 1,ncol = 1)
              rownames(weights[[gname]][["wgt"]]) <- snpnames[nextindex]
              colnames(weights[[gname]][["wgt"]]) <- "weight"
            }
            else{
              weights[[gname]][["wgt"]] <- matrix(tmp_wgt,nrow = length(tmp_wgt),ncol = 1)
              rownames(weights[[gname]][["wgt"]]) <- snpnames[nextindex]
              colnames(weights[[gname]][["wgt"]]) <- "weight"
            }
            #add gene to next region
            #regionlist[[region_tag_next]][["gidx"]] <- c(regionlist[[region_tag_next]][["gidx"]],tmp_region[["gidx"]][which(gnames==gname)])
            regionlist[[region_tag_next]][["gid"]] <- c(regionlist[[region_tag_next]][["gid"]],gname)
            #remove gene from this region
            #regionlist[[region_tag_current]][["gidx"]] <- regionlist[[region_tag_current]][["gidx"]][which(gnames!=gname)]
            regionlist[[region_tag_current]][["gid"]] <- regionlist[[region_tag_current]][["gid"]][regionlist[[region_tag_current]][["gid"]]!=gname]
          }
          else{
            #modify weights file - drop weights in other regions
            tmp_wgt <- weights[[gname]][["wgt"]][thisindex]
            if(length(tmp_wgt)==1){
              weights[[gname]][["wgt"]] <- matrix(tmp_wgt,nrow = 1,ncol = 1)
              rownames(weights[[gname]][["wgt"]]) <- snpnames[thisindex]
              colnames(weights[[gname]][["wgt"]]) <- "weight"
            }
            else{
              weights[[gname]][["wgt"]] <- matrix(tmp_wgt,nrow = length(tmp_wgt),ncol = 1)
              rownames(weights[[gname]][["wgt"]]) <- snpnames[thisindex]
              colnames(weights[[gname]][["wgt"]]) <- "weight"
            }
          }
        }
      }
    }
  }
  return(list(regionlist=regionlist,weights=weights,boundary_genes=boundary_genes))
}
