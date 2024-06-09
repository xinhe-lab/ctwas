# @title Runs EM to estimate parameters.
#
# @param region_data a list object with the susie input data for each region
#
# @param niter the number of iterations of the E-M algorithm to perform
#
# @param init_group_prior a vector of initial prior inclusion probabilities for SNPs and genes.
#
# @param init_group_prior_var a vector of initial prior variances for SNPs and gene effects.
#
# @param group_prior_var_structure a string indicating the structure to put on the prior variance parameters.
# "shared_type" allows all groups in one molecular QTL type to share the same variance parameter.
# "shared_context" allows all groups in one context (tissue, cell type, condition) to share the same variance parameter.
# "shared_nonSNP" allows all non-SNP groups to share the same variance parameter.
# "shared_all" allows all groups to share the same variance parameter.
# "independent" allows all groups to have their own separate variance parameters.
#
# @param use_null_weight TRUE/FALSE. If TRUE, allow for a probability of no effect in susie
#
# @param ncore The number of cores used to parallelize susie over regions
#
# @param verbose TRUE/FALSE. If TRUE, print detail messages
#
# @return a list of parameters
#
#' @importFrom logging loginfo
#' @importFrom parallel mclapply
EM_est_param <- function(
    region_data,
    niter = 20,
    init_group_prior = NULL,
    init_group_prior_var = NULL,
    group_prior_var_structure = c("shared_type", "shared_context", "shared_nonSNP", "shared_all", "independent"),
    use_null_weight = TRUE,
    ncore = 1,
    verbose = FALSE,
    ...){

  # get groups, types and contexts from region_data
  groups <- unique(unlist(lapply(region_data, "[[", "gs_group")))
  types <- unique(unlist(lapply(region_data, "[[", "gs_type")))
  contexts <- unique(unlist(lapply(region_data, "[[", "gs_context")))

  # set pi_prior and V_prior based on init_group_prior and init_group_prior_var
  res <- initiate_group_priors(init_group_prior, init_group_prior_var, groups)
  pi_prior <- res$pi_prior
  V_prior <- res$V_prior
  rm(res)

  # store estimated group priors from all iterations
  group_prior_iters <- matrix(NA, nrow = length(groups), ncol = niter)
  rownames(group_prior_iters) <- groups
  colnames(group_prior_iters) <- paste0("iter", 1:ncol(group_prior_iters))

  group_prior_var_iters <- matrix(NA, nrow = length(groups), ncol = niter)
  rownames(group_prior_var_iters) <- groups
  colnames(group_prior_var_iters) <- paste0("iter", 1:ncol(group_prior_var_iters))

  region_ids <- names(region_data)
  for (iter in 1:niter){
    if (verbose){
      loginfo("Start EM iteration %d ...", iter)
    }

    EM_susie_res_list <- mclapply(region_ids, function(region_id){
      ctwas_susie_region_L1(region_data, region_id,
                            pi_prior, V_prior,
                            use_null_weight = use_null_weight,
                            ...)
    }, mc.cores = ncore)
    if (length(EM_susie_res_list) != length(region_ids)) {
      stop("Not all cores returned results. Try rerun with bigger memory or fewer cores")
    }
    EM_susie_res <- do.call(rbind, EM_susie_res_list)

    # update estimated group_prior from the current iteration
    pi_prior <- sapply(names(pi_prior), function(x){mean(EM_susie_res$susie_pip[EM_susie_res$group==x])})
    group_prior_iters[names(pi_prior),iter] <- pi_prior

    if (verbose){
      loginfo("Iteration %d, group_prior {%s}: {%s}", iter, names(pi_prior), format(pi_prior, digits = 4))
    }

    # update estimated group_prior_var from the current iteration
    # in susie, mu2 is under standardized scale (if performed)
    # e.g. X = matrix(rnorm(n*p),nrow=n,ncol=p)
    # y = X %*% beta + rnorm(n)
    # res = susie(X,y,L=10)
    # X2 = 5*X
    # res2 = susie(X2,y,L=10)
    # res$mu2 is identical to res2$mu2 but coefficients are on diff scale.
    if (group_prior_var_structure=="independent") {
      V_prior <- sapply(names(V_prior),
                        function(x){
                          tmp_EM_susie_res <- EM_susie_res[EM_susie_res$group==x,];
                          sum(tmp_EM_susie_res$susie_pip*tmp_EM_susie_res$mu2)/sum(tmp_EM_susie_res$susie_pip)})
    } else if (group_prior_var_structure=="shared_nonSNP") {
      tmp_EM_susie_res <- EM_susie_res[EM_susie_res$group=="SNP",]
      V_prior["SNP"] <- sum(tmp_EM_susie_res$susie_pip*tmp_EM_susie_res$mu2)/sum(tmp_EM_susie_res$susie_pip)
      tmp_EM_susie_res <- EM_susie_res[EM_susie_res$group!="SNP",]
      V_prior[names(V_prior)!="SNP"] <- sum(tmp_EM_susie_res$susie_pip*tmp_EM_susie_res$mu2)/sum(tmp_EM_susie_res$susie_pip)
    } else if (group_prior_var_structure=="shared_all") {
      V_prior[names(V_prior)] <- sum(tmp_EM_susie_res$susie_pip*tmp_EM_susie_res$mu2)/sum(tmp_EM_susie_res$susie_pip)
    } else if (group_prior_var_structure=="shared_type") {
      tmp_EM_susie_res <- EM_susie_res[EM_susie_res$group=="SNP",]
      V_prior["SNP"] <- sum(tmp_EM_susie_res$susie_pip*tmp_EM_susie_res$mu2)/sum(tmp_EM_susie_res$susie_pip)
      nonSNP_types <- setdiff(types, "SNP")
      for(type in nonSNP_types){
        tmp_EM_susie_res <- EM_susie_res[EM_susie_res$type==type,]
        V_prior[sapply(names(V_prior), function(x){
          unlist(strsplit(x, "[|]"))[1]})==type] <-
          sum(tmp_EM_susie_res$susie_pip*tmp_EM_susie_res$mu2)/sum(tmp_EM_susie_res$susie_pip)
      }
    } else if (group_prior_var_structure=="shared_context") {
      tmp_EM_susie_res <- EM_susie_res[EM_susie_res$group=="SNP",]
      V_prior["SNP"] <- sum(tmp_EM_susie_res$susie_pip*tmp_EM_susie_res$mu2)/sum(tmp_EM_susie_res$susie_pip)
      nonSNP_contexts <- setdiff(contexts, "SNP")
      for(context in nonSNP_contexts){
        tmp_EM_susie_res <- EM_susie_res[EM_susie_res$context==context,]
        V_prior[sapply(names(V_prior), function(x){
          unlist(strsplit(x, "[|]"))[2]})==context] <-
          sum(tmp_EM_susie_res$susie_pip*tmp_EM_susie_res$mu2)/sum(tmp_EM_susie_res$susie_pip)
      }
    }
    group_prior_var_iters[names(V_prior), iter] <- V_prior

    if (verbose){
      loginfo("Iteration %d, group_prior_var {%s}: {%s}", iter, names(V_prior), format(V_prior, digits = 4))
    }
  }

  group_size <- table(EM_susie_res$group)
  group_size <- group_size[rownames(group_prior_iters)]

  return(list("group_prior"= pi_prior,
              "group_prior_var" = V_prior,
              "group_prior_iters" = group_prior_iters,
              "group_prior_var_iters" = group_prior_var_iters,
              "group_prior_var_structure" = group_prior_var_structure,
              "group_size" = group_size))
}

# run susie for one region with L = 1 without LD matrix
ctwas_susie_region_L1 <- function(region_data, region_id,
                                  pi_prior, V_prior,
                                  use_null_weight = TRUE, ...){
  # load susie input data
  regiondata <- region_data[[region_id]]
  sid <- regiondata[["sid"]]
  gid <- regiondata[["gid"]]
  z <- regiondata[["z"]]
  gs_group <- regiondata[["gs_group"]]
  g_type <- regiondata[["g_type"]]
  g_context <- regiondata[["g_context"]]
  g_group <- regiondata[["g_group"]]

  if (length(z) < 2) {
    stop(paste(length(z), "variables in the region. At least two variables in a region are needed to run susie"))
  }

  # update priors, prior variances and null_weight based on the estimated group_prior and group_prior_var from the previous iteration
  res <- set_region_susie_priors(pi_prior, V_prior, gs_group, L = 1, use_null_weight = use_null_weight)
  prior <- res$prior
  V <- res$V
  null_weight <- res$null_weight
  rm(res)

  # Use an identity matrix as LD, R does not matter for susie when L = 1
  R <- diag(length(z))

  # in susie, prior_variance is under standardized scale (if performed)
  susie_res <- ctwas_susie_rss(z = z,
                               R = R,
                               prior_weights = prior,
                               prior_variance = V,
                               L = 1,
                               null_weight = null_weight,
                               max_iter = 1,
                               ...)
  # annotate susie result
  susie_res_df <- anno_susie(susie_res,
                             gids = gid,
                             sids = sid,
                             g_type = g_type,
                             g_context = g_context,
                             g_group = g_group,
                             region_id = region_id,
                             include_cs_index = FALSE)

  return(susie_res_df)
}
