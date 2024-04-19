
#' set pi_prior and V_prior based on init_group_prior and init_group_prior_var
initiate_group_priors <- function(group_prior = NULL, group_prior_var = NULL, groups) {

  if (is.null(group_prior)){
    group_prior <- structure(as.numeric(rep(NA,length(groups))), names=groups)
  }

  if (is.null(group_prior_var)){
    group_prior_var <- structure(as.numeric(rep(NA,length(groups))), names=groups)
  }

  if ( !setequal(names(group_prior), groups) || !setequal(names(group_prior_var), groups) ) {
    stop("names of group_prior or group_prior_var do not match with groups")
  }

  pi_prior <- list()
  V_prior <- list()
  for (group in groups){
    pi_prior[[group]] <- unname(group_prior[group])
    V_prior[[group]] <- unname(group_prior_var[group])
  }
  pi_prior <- unlist(pi_prior)
  V_prior <- unlist(V_prior)

  return(list(pi_prior = pi_prior,
              V_prior = V_prior))
}


# set prior and prior variance values for the region
set_region_susie_priors <- function(pi_prior, V_prior, gs_group, L, use_null_weight = TRUE){

  p <- length(gs_group)

  if (any(is.na(pi_prior))){
    prior <- rep(1/p, p)
  } else {
    prior <- unname(pi_prior[gs_group])
  }

  if (any(is.na(V_prior))){
    V <- matrix(rep(50, L * p), nrow = L)
    # following the default in susieR::susie_rss
  } else{
    V <- unname(V_prior[gs_group])
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
