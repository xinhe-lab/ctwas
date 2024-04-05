
#' adjust regionlist for boundary genes
adjust_boundary_genes <- function(boundary_genes, region_info, weights, regionlist){

  for (i in 1:nrow(boundary_genes)){
    gname <- boundary_genes[i, "id"]
    region_tags <- unlist(strsplit(boundary_genes[i, "region_tag"], split = ";"))
    wgt <- weights[[gname]][["wgt"]]

    region_r2 <- sapply(region_tags, function(region_tag){
      ld_snpinfo <- ctwas:::read_LD_SNP_files(region_info[region_info$region_tag == region_tag, "SNP_info"])
      sum(wgt[which(rownames(wgt) %in% ld_snpinfo$id)]^2)})

    # assign boundary gene to the region with max r2, and remove it from other regions
    selected_region_tag <- region_tags[which.max(region_r2)]
    unselected_region_tags <- setdiff(region_tags, selected_region_tag)
    regionlist[[selected_region_tag]][["gid"]] <- unique(c(regionlist[[selected_region_tag]][["gid"]],gname))
    for(unselected_region_tag in unselected_region_tags){
      regionlist[[unselected_region_tag]][["gid"]] <- regionlist[[unselected_region_tag]][["gid"]][regionlist[[unselected_region_tag]][["gid"]]!=gname]
    }
  }

  return(regionlist)
}

