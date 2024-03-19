
#' filter regions based on probability of at most 1 causal effect
filter_regions <- function(regionlist, zdf, group_prior, prob_single = 0.8){
  regionlist2 <- regionlist
  for (region_tag in names(regionlist)){
      gid <- regionlist[[region_tag]][["gid"]]
      sid <- regionlist[[region_tag]][["sid"]]
      gs_type <- zdf$type[match(c(gid,sid), zdf$id)]

      group_size <- table(gs_type)[names(group_prior)]
      group_size[is.na(group_size)] <- 0

      P1 <- prod((1-group_prior)^group_size) * (1 + sum(group_size*(group_prior/(1-group_prior))))

      if (P1 < prob_single){
        regionlist2[[region_tag]] <- NULL
      }
  }
  regionlist2
}

#' assign regions to cores
#'
#' @param ncore integer, numeber of cores, at least 1
#' regions allocated to given number of cores
#' regionlist need to contain at least 1 non-empty
#'
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


