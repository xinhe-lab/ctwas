#' Select single effect regions
select_single_effect_regions <- function(regionlist, group_prior, p_single_effect = 0.8){
  region_ids <- names(regionlist)
  p1 <- sapply(region_ids, function(x){
    group_size <- table(regionlist[[x]][["gs_group"]])
    group_size <- group_size[names(group_prior)]
    group_size[is.na(group_size)] <- 0
    p1 <- prod((1-group_prior)^group_size) * (1 + sum(group_size*(group_prior/(1-group_prior))))
    return(p1)
  })
  selected_region_ids <- region_ids[p1 >= p_single_effect]
  return(selected_region_ids)
}

#' Select regions with high non-SNP PIPs
select_highPIP_regions <- function(finemap_res, min_nonSNP_PIP = 0.5){
  region_ids <- unique(finemap_res$region_id)
  nonSNP_PIPs <- sapply(region_ids, function(x){
    finemap_region_res <- finemap_res[finemap_res$region_id == x,]
    nonSNP_PIP <- sum(finemap_region_res$susie_pip[finemap_region_res$type != "SNP"])
    nonSNP_PIP[is.na(nonSNP_PIP)] <- 0 # 0 if nonSNP_PIP is NA
    return(nonSNP_PIP)
  })
  selected_region_ids <- region_ids[nonSNP_PIPs >= min_nonSNP_PIP]
  return(selected_region_ids)
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
