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
#' @param thin  A scalar in (0,1]. The proportion of SNPs
#' left after down sampling. Only applied on SNPs after selecting variants.
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
#' @param minvar minimum number of variatns in a region
#'
#' @param merge TRUE/FALSE. If TRUE, merge regions when a gene spans a region boundary (i.e. belongs to multiple regions.)
#' 
#' @param outputdir a string, the directory to store output
#' 
#' @param outname a string, the output name
#' 
#' @param ncore the number of cores used to parallelize region indexing
#' 
#' @param reuse_R_gene an option to reuse the R_gene matrix when indexing for the final rerun step
#'
#' @return A list. Items correspond to each pvarf/exprvarf. Each Item is
#'  also a list, the items in this list are for each region.
#'
#' @importFrom logging loginfo
#'
adjust_boundary <- function(regions, weight_list, regionlist){
  boundary_genes <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(boundary_genes) <- c("gene","chr","region1","region2")
  for (rn in 1:(nrow(regions)-1)){
    current <- regionlist[[as.character(rn)]]
    nextone <- regionlist[[as.character(rn+1)]]
    gnames <- regionlist[[as.character(rn)]][["gid"]]
    tmp_region <- regionlist[[as.character(rn)]]
    if(length(gnames>0)){
      ld_snpinfo <- ctwas:::read_LD_SNP_file(regionlist[[as.character(rn)]][["SNP_info"]])
      for (i in 1:length(gnames)){
        gname <- gnames[i]
        wgt <- weight_list[[gname]] 
        snpnames <- rownames(wgt)
        ld.idx <- match(snpnames, ld_snpinfo$id)
        if(anyNA(ld.idx)){ # QTLs are across boundary
          boundary_genes <- rbind(boundary_genes,data.frame("gene"=gname,"chr"=regions$chr[rn],"region1"=rn,"region2"=rn+1)) 
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
            regionlist[[as.character(rn+1)]][["gidx"]] <- c(regionlist[[as.character(rn+1)]][["gidx"]],tmp_region[["gidx"]][which(gnames==gname)])
            regionlist[[as.character(rn+1)]][["gid"]] <- c(regionlist[[as.character(rn+1)]][["gid"]],gname)
            #remove gene from this region
            regionlist[[as.character(rn)]][["gidx"]] <- regionlist[[as.character(rn)]][["gidx"]][which(gnames!=gname)]
            regionlist[[as.character(rn)]][["gid"]] <- regionlist[[as.character(rn)]][["gid"]][!regionlist[[as.character(rn)]][["gid"]]==gname]
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