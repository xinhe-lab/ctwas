#' select single effect regions
select_single_effect_regions <- function(regionlist, group_prior, p_single_effect = 0.8, ncore=ncore){
  loginfo("Select regions with P(single effect) >= %s", p_single_effect)

  region_tags <- names(regionlist)
  cl <- parallel::makeCluster(ncore, outfile = "", type = "FORK")
  doParallel::registerDoParallel(cl)

  selected_region_tags <- foreach(region_tag = region_tags, .combine = "c") %dopar% {
    gid <- regionlist[[region_tag]][["gid"]]
    sid <- regionlist[[region_tag]][["sid"]]
    gs_type <- regionlist[[region_tag]][["gs_type"]]

    group_size <- table(gs_type)[names(group_prior)]
    group_size[is.na(group_size)] <- 0

    P1 <- prod((1-group_prior)^group_size) * (1 + sum(group_size*(group_prior/(1-group_prior))))
    if (P1 >= p_single_effect){
      region_tag
    }
  }
  parallel::stopCluster(cl)

  loginfo("%d regions selected", length(selected_region_tags))
  selected_regionlist <- regionlist[selected_region_tags]
  return(selected_regionlist)
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


