
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

get_region_tags <- function(regionlist){
  region_tags <- NULL
  for (b in 1:length(regionlist)){
    if (length(regionlist[[b]]) > 0){
      rn <- names(regionlist[[b]])
      region_tags <- c(region_tags, paste0(b, "_", rn))
    }
  }
  return(region_tags)
}

# select and assemble a subset of regionlist by region_tags
subset_regionlist <- function(regionlist, region_tags = NULL){

  if (length(region_tags) == 0){
    stop("No regions selected!")
  }

  # subset regionlist by region_tags
  region_tags_df <- data.frame(
    b = sapply(strsplit(region_tags, "_"), "[[", 1),
    rn = sapply(strsplit(region_tags, "_"), "[[", 2))

  regionlist_subset <- vector("list", length = 22)
  for (b in region_tags_df$b) {
    rn <- region_tags_df[region_tags_df$b == b, "rn"]
    regionlist_subset[[b]] <- regionlist[[b]][as.character(rn)]
  }

  # coordinates of the subset regions
  temp_regs <- lapply(1:22, function(x) cbind(x,
                                              unlist(lapply(regionlist_subset[[x]], "[[", "start")),
                                              unlist(lapply(regionlist_subset[[x]], "[[", "stop"))))

  regs_subset <- do.call(rbind, lapply(temp_regs, function(x) if (ncol(x) == 3){x}))
  loginfo("Subset %d regions from the regionlist", nrow(regs_subset))

  return(list("regionlist" = regionlist_subset,
              "regs" = regs))
}


