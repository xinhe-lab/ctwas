# identify cross-boundary genes, adjust regionlist and update weigh_list
adjust_boundary <- function(regioninfo, weight_list, regionlist){
  boundary_genes <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(boundary_genes) <- c("gene","chrom","region1","region2")

  for (i in 1:(nrow(regioninfo)-1)){
    region_tag_current <- regioninfo$region_tag[i]
    region_tag_next <- regioninfo$region_tag[i+1]
    current <- regionlist[[region_tag_current]]
    nextone <- regionlist[[region_tag_next]]
    gnames <- regionlist[[region_tag_current]][["gid"]]
    tmp_region <- regionlist[[region_tag_current]]
    if(length(gnames>0)){
      ld_snpinfo <- read_LD_SNP_file(regionlist[[region_tag_current]][["SNP_info"]]) #ctwas:::
      for (i in 1:length(gnames)){
        gname <- gnames[i]
        wgt <- weight_list[[gname]]
        snpnames <- rownames(wgt)
        ld.idx <- match(snpnames, ld_snpinfo$id)
        if(anyNA(ld.idx)){ # QTLs are across boundary
          boundary_genes <- rbind(boundary_genes,
                                  data.frame("gene"=gname,
                                             "region1"=paste0(regioninfo$chrom[i],"_",region_tag_current),
                                             "region2"=paste0(regioninfo$chrom[i],"_",region_tag_next)))
          thisindex <- !is.na(ld.idx)
          nextindex <- is.na(ld.idx)
          thisr2 <- sum(wgt[thisindex]^2)
          nextr2 <- sum(wgt[nextindex]^2)
          if(thisr2<nextr2){
          #modify weights file - drop weights in other regions
            tmp_wgt <- weight_list[[gname]][nextindex]
            if(length(tmp_wgt)==1){
              weight_list[[gname]] <- matrix(tmp_wgt,nrow = 1,ncol = 1)
              rownames(weight_list[[gname]]) <- snpnames[nextindex]
              colnames(weight_list[[gname]]) <- "weight"
            }
            else{
              weight_list[[gname]] <- matrix(tmp_wgt,nrow = length(tmp_wgt),ncol = 1)
              rownames(weight_list[[gname]]) <- snpnames[nextindex]
              colnames(weight_list[[gname]]) <- "weight"
            }
            #add gene to next region
            #regionlist[[region_tag_next]][["gidx"]] <- c(regionlist[[region_tag_next]][["gidx"]],tmp_region[["gidx"]][which(gnames==gname)])
            regionlist[[region_tag_next]][["gid"]] <- c(regionlist[[region_tag_next]][["gid"]],gname)
            #remove gene from this region
            #regionlist[[region_tag_current]][["gidx"]] <- regionlist[[region_tag_current]][["gidx"]][which(gnames!=gname)]
            regionlist[[region_tag_current]][["gid"]] <- regionlist[[region_tag_current]][["gid"]][!regionlist[[region_tag_current]][["gid"]]==gname]
          }
          else{
            #modify weights file - drop weights in other regions
            tmp_wgt <- weight_list[[gname]][thisindex]
            if(length(tmp_wgt)==1){
              weight_list[[gname]] <- matrix(tmp_wgt,nrow = 1,ncol = 1)
              rownames(weight_list[[gname]]) <- snpnames[thisindex]
              colnames(weight_list[[gname]]) <- "weight"
            }
            else{
              weight_list[[gname]] <- matrix(tmp_wgt,nrow = length(tmp_wgt),ncol = 1)
              rownames(weight_list[[gname]]) <- snpnames[thisindex]
              colnames(weight_list[[gname]]) <- "weight"
            }
          }
        }
      }
    }
  }
  return(list(regionlist=regionlist,weight_list=weight_list,boundary_genes=boundary_genes))
}
