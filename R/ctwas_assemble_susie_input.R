#' assemble susie input data for all regions
assemble_susie_input_list <- function(zdf,
                                      regionlist,
                                      L,
                                      group_prior = NULL,
                                      group_prior_var = NULL,
                                      use_null_weight = TRUE,
                                      ncore = 1){
  # set pi_prior and V_prior based on group_prior and group_prior_var
  types <- unique(zdf$type)
  if (is.null(group_prior)){
    group_prior <- structure(as.numeric(rep(NA,length(types))), names=types)
  }
  if (is.null(group_prior_var)){
    group_prior_var <- structure(as.numeric(rep(NA,length(types))), names=types)
  }
  pi_prior <- list()
  V_prior <- list()
  for (type in types){
    pi_prior[[type]] <- unname(group_prior[type])
    V_prior[[type]] <- unname(group_prior_var[type])
  }
  pi_prior <- unlist(pi_prior)
  V_prior <- unlist(V_prior)

  # start running EM iterations
  cl <- parallel::makeCluster(ncore, outfile = "")
  doParallel::registerDoParallel(cl)

  corelist <- region2core(regionlist, ncore)

  susie_input_list <- foreach (core = 1:length(corelist), .combine = "c", .packages = "ctwas") %dopar% {
    susie_input.core <- list()
    region_tags.core <- corelist[[core]]
    for (region_tag in region_tags.core) {
      susie_input <- assemble_region_susie_input(sid = regionlist[[region_tag]][["sid"]],
                                                 gid = regionlist[[region_tag]][["gid"]],
                                                 zdf = zdf,
                                                 pi_prior,
                                                 V_prior,
                                                 L = L,
                                                 use_null_weight = use_null_weight)
      susie_input.core[[region_tag]] <- susie_input
    }
    susie_input.core
  }

  return(susie_input_list)
}

#' assemble susie input data for a region
assemble_region_susie_input <- function(sid, gid, zdf, pi_prior, V_prior, L, use_null_weight=TRUE) {

  types <- unique(zdf$type)

  # assemble zscores
  # keep only GWAS SNPs and imputed genes
  sid <- intersect(sid, zdf$id)
  gid <- intersect(gid, zdf$id)
  g_idx <- match(gid, zdf$id)
  s_idx <- match(sid, zdf$id)
  gs_idx <- c(g_idx, s_idx)
  z <- zdf[gs_idx, "z"]
  g_type <- zdf$type[g_idx]
  g_QTLtype <- zdf$QTLtype[g_idx]
  gs_type <- zdf$type[gs_idx]

  if (anyNA(z))
    loginfo("Warning: z-scores contains missing values!")

  # set prior and prior variance values for the region
  res <- set_region_susie_priors(pi_prior, V_prior, gs_type, L, use_null_weight = use_null_weight)
  prior <- res$prior
  V <- res$V
  null_weight <- res$null_weight
  rm(res)

  susie_input <- list(sid = sid,
                      gid = gid,
                      z = z,
                      g_type = g_type,
                      g_QTLtype = g_QTLtype,
                      gs_type = gs_type,
                      prior = prior,
                      V = V,
                      null_weight = null_weight)

  return(susie_input)
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
