#' Select single effect regions
compute_region_p_single_effect <- function(regionlist, group_prior){
  region_ids <- names(regionlist)
  p1 <- sapply(region_ids, function(x){
    group_size <- table(regionlist[[x]][["gs_group"]])
    group_size <- group_size[names(group_prior)]
    group_size[is.na(group_size)] <- 0
    p1 <- prod((1-group_prior)^group_size) * (1 + sum(group_size*(group_prior/(1-group_prior))))
    p1
  })
  return(data.frame(region_id = region_ids, p1 = p1))

}

#' Select regions with high non-SNP PIPs
compute_region_nonSNP_PIPs <- function(finemap_res){
  region_ids <- unique(finemap_res$region_id)
  nonSNP_PIPs <- sapply(region_ids, function(x){
    finemap_region_res <- finemap_res[finemap_res$region_id == x,]
    nonSNP_PIP <- sum(finemap_region_res$susie_pip[finemap_region_res$type != "SNP"])
    nonSNP_PIP[is.na(nonSNP_PIP)] <- 0 # 0 if nonSNP_PIP is NA
    nonSNP_PIP
  })
  return(data.frame(region_id = region_ids, nonSNP_PIP = nonSNP_PIPs))
}

#' assign regions to cores
region2core <- function(regionlist, ncore = 1){
  region_ids <- names(regionlist)
  if (ncore > 1) {
    d <- cut(1:length(region_ids), ncore, labels = FALSE)
    corelist <- split(region_ids,d)
  } else {
    corelist <- list("1" = region_ids)
  }
  return(corelist)
}
