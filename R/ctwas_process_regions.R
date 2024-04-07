#' Select single effect regions
select_single_effect_regions <- function(regionlist, group_prior, p_single_effect = 0.8){
  region_tags <- names(regionlist)
  p1 <- sapply(region_tags, function(x){
    group_size <- table(regionlist[[x]][["gs_type"]])
    group_size <- group_size[names(group_prior)]
    group_size[is.na(group_size)] <- 0
    p1 <- prod((1-group_prior)^group_size) * (1 + sum(group_size*(group_prior/(1-group_prior))))
    return(p1)
  })
  selected_region_tags <- region_tags[p1 >= p_single_effect]
  return(selected_region_tags)
}

#' Select regions with high non-SNP PIPs
select_highPIP_regions <- function(finemap_res, region_tags, min_nonSNP_PIP = 0.5){
  nonSNP_PIPs <- sapply(region_tags, function(x){
    finemap_region_res <- finemap_res[finemap_res$region_tag == x,]
    nonSNP_PIP <- sum(finemap_region_res$susie_pip[finemap_region_res$type != "SNP"])
    nonSNP_PIP[is.na(nonSNP_PIP)] <- 0 # 0 if nonSNP_PIP is NA
    return(nonSNP_PIP)
  })
  selected_region_tags <- region_tags[nonSNP_PIPs >= min_nonSNP_PIP]
  return(selected_region_tags)
}

#' assign regions to cores
region2core <- function(regionlist, ncore = 1){
  region_tags <- names(regionlist)
  if (ncore > 1) {
    d <- cut(1:length(region_tags), ncore, labels = FALSE)
    corelist <- split(region_tags,d)
  } else {
    corelist <- list("1" = region_tags)
  }
  return(corelist)
}
