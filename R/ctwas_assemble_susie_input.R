#' add z-scores from z_snp and z_gene to regionlist
add_z_to_regionlist <- function(regionlist,
                                z_snp,
                                z_gene,
                                ncore = 1){

  # Combine z-scores from z_snp and z_gene
  zdf <- combine_z(z_snp, z_gene)

  cl <- parallel::makeCluster(ncore, outfile = "", type = "FORK")
  doParallel::registerDoParallel(cl)
  corelist <- region2core(regionlist, ncore)

  region_tags <- names(regionlist)
  loginfo("Add z-scores to regionlist for %d regions", length(region_tags))

  regionlist2 <- foreach (core = 1:length(corelist), .combine = "c") %dopar% {
    regionlist2.core <- list()
    region_tags.core <- corelist[[core]]
    for (region_tag in region_tags.core) {
      # add z-scores and types of the region to the regionlist
      regionlist2.core[[region_tag]] <- regionlist[[region_tag]]
      sid <- regionlist[[region_tag]][["sid"]]
      gid <- regionlist[[region_tag]][["gid"]]
      g_idx <- match(gid, zdf$id)
      s_idx <- match(sid, zdf$id)
      gs_idx <- c(g_idx, s_idx)
      regionlist2.core[[region_tag]][["z"]] <- zdf[gs_idx, "z"]
      regionlist2.core[[region_tag]][["g_type"]] <- zdf$type[g_idx]
      regionlist2.core[[region_tag]][["g_QTLtype"]] <- zdf$QTLtype[g_idx]
      regionlist2.core[[region_tag]][["gs_type"]] <- zdf$type[gs_idx]
    }
    regionlist2.core
  }
  parallel::stopCluster(cl)

  return(regionlist2)
}


# set prior and prior variance values for the region
set_region_susie_priors <- function(pi_prior, V_prior, gs_type, L, use_null_weight = TRUE){

  p <- length(gs_type)

  if (any(is.na(pi_prior))){
    prior <- rep(1/p, p)
  } else {
    prior <- unname(pi_prior[gs_type])
  }

  if (any(is.na(V_prior))){
    V <- matrix(rep(50, L * p), nrow = L)
    # following the default in susieR::susie_rss
  } else{
    V <- unname(V_prior[gs_type])
    V <- matrix(rep(V, each = L), nrow=L)
  }

  if (isTRUE(use_null_weight)){
    null_weight <- max(0, 1 - sum(prior))
    prior <- prior/(1-null_weight)
  } else {
    null_weight <- NULL
  }

  return(list(prior = prior,
              V = V,
              null_weight = null_weight))

}
