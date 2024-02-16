
#' Merge regions.
#'
#' @param region_info a data frame of region definition and associated file names
#'
#' @param regionlist the list of regions to be merged.
#'
#' @importFrom logging loginfo
#'
#' @return a list of merged regions.
#'
merge_regions <- function(region_info,
                          regionlist) {

}

#' filter regions based on probability of at most 1 causal effect
filter_regions <- function(regionlist, group_prior, prob_single = 0.8, zdf){
  regionlist2 <- regionlist
  for (b in 1: length(regionlist)){
    for (rn in names(regionlist[[b]])) {
      gid <- regionlist[[b]][[rn]][["gid"]]
      sid <- regionlist[[b]][[rn]][["sid"]]
      gs_type <- zdf$type[match(c(gid,sid), zdf$id)]

      group_size <- table(gs_type)[names(group_prior)]
      group_size[is.na(group_size)] <- 0

      P1 <- prod((1-group_prior)^group_size) * (1 + sum(group_size*(group_prior/(1-group_prior))))

      if (P1 < prob_single){
        regionlist2[[b]][[rn]] <- NULL
      }
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
  dflist <- list()
  for (b in 1:length(regionlist)){
    if (length(regionlist[[b]]) > 0){
      dflist[[b]] <- data.frame("b" = b, "rn"= names(regionlist[[b]]), stringsAsFactors = FALSE)
    }
  }
  df <- do.call(rbind, dflist)
  if (ncore > 1) {
    d <- cut(1:nrow(df), ncore, labels = FALSE)
    corelist <- split(df,d)
  } else {
    corelist <- list("1" = df)
  }
  return(corelist)
}
