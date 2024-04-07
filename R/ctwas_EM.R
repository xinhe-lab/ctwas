#' Run EM to estimate parameters.
#' Iteratively run susie and estimate parameters - RSS version
#'
#' @param regionlist a list object with the susie input data for each region
#'
#' @param niter the number of iterations of the E-M algorithm to perform
#'
#' @param init_group_prior a vector of initial prior inclusion probabilities for SNPs and genes.
#'
#' @param init_group_prior_var a vector of initial prior variances for SNPs and gene effects.
#'
#' @param group_prior_var_structure a string indicating the structure to put on the prior variance parameters.
#' "independent" is the default and allows all groups to have their own separate variance parameters.
#' "shared_all" allows all groups to share the same variance parameter.
#' "shared_type" allows all groups in one molecular QTL type to share the same variance parameter.
#'
#' @param use_null_weight TRUE/FALSE. If TRUE, allow for a probability of no effect in susie
#'
#' @param coverage A number between 0 and 1 specifying the \dQuote{coverage} of the estimated confidence sets
#'
#' @param min_abs_corr Minimum absolute correlation allowed in a
#'   credible set. The default, 0.5, corresponds to a squared
#'   correlation of 0.25, which is a commonly used threshold for
#'   genotype data in genetic studies.
#'
#' @param max_iter Maximum number of IBSS iterations to perform.
#'
#' @param ncore The number of cores used to parallelize susie over regions
#'
#' @param ncore The number of cores used to parallelize susie over regions
#'
#' @param verbose TRUE/FALSE. If TRUE, print detail messages
#'
#' @importFrom logging loginfo
#' @importFrom foreach %dopar% foreach
#'
#' @return a list of parameters
#'
EM <- function(regionlist,
               niter = 20,
               init_group_prior = NULL,
               init_group_prior_var = NULL,
               group_prior_var_structure = c("independent","shared_all","shared_QTLtype"),
               use_null_weight = TRUE,
               max_iter = 1,
               ncore = 1,
               verbose = FALSE,
               ...){

  # get types and QTLtypes from regionlist
  types <- unique(unlist(lapply(regionlist, "[[", "gs_type")))
  QTLtypes <- unique(unlist(lapply(regionlist, "[[", "gs_QTLtype")))

  # set pi_prior and V_prior based on init_group_prior and init_group_prior_var
  res <- initiate_group_priors(init_group_prior, init_group_prior_var, types)
  pi_prior <- res$pi_prior
  V_prior <- res$V_prior
  rm(res)

  # store estimated group priors from each iteration
  group_prior_rec <- matrix(NA, nrow = length(types), ncol = niter)
  rownames(group_prior_rec) <- types
  group_prior_var_rec <- matrix(NA, nrow = length(types), ncol = niter)
  rownames(group_prior_var_rec) <- types

  # start running EM iterations
  cl <- parallel::makeCluster(ncore, outfile = "", type = "FORK")
  doParallel::registerDoParallel(cl)

  corelist <- region2core(regionlist, ncore)

  for (iter in 1:niter){
    loginfo("Start EM iteration %d", iter)
    EM_susie_res <- foreach (core = 1:length(corelist), .combine = "rbind", .packages = "ctwas") %dopar% {
      susie_res.core.list <- list()
      # run susie for each region
      region_tags.core <- corelist[[core]]
      for (region_tag in region_tags.core) {
        # load susie input data
        if (verbose)
          loginfo("load susie input data for region %s", region_tag)
        sid <- regionlist[[region_tag]][["sid"]]
        gid <- regionlist[[region_tag]][["gid"]]
        z <- regionlist[[region_tag]][["z"]]
        gs_type <- regionlist[[region_tag]][["gs_type"]]
        g_type <- regionlist[[region_tag]][["g_type"]]
        g_QTLtype <- regionlist[[region_tag]][["g_QTLtype"]]

        # update priors, prior variances and null_weight based on the estimated group_prior and group_prior_var from the previous iteration
        if (verbose)
          loginfo("update priors, prior variances for region %s", region_tag)
        res <- set_region_susie_priors(pi_prior, V_prior, gs_type, L = 1, use_null_weight = use_null_weight)
        prior <- res$prior
        V <- res$V
        null_weight <- res$null_weight
        rm(res)

        # Use an identity matrix as LD, R does not matter for susie when L = 1
        R <- diag(length(z))

        # in susie, prior_variance is under standardized scale (if performed)
        if (verbose)
          loginfo("run susie for region %s", region_tag)
        susie_res <- ctwas_susie_rss(z = z,
                                     R = R,
                                     prior_weights = prior,
                                     prior_variance = V,
                                     L = 1,
                                     null_weight = null_weight,
                                     max_iter = max_iter,
                                     ...)
        if (verbose)
          loginfo("annotate susie result for region %s", region_tag)
        # annotate susie result
        susie_res_df <- anno_susie(susie_res,
                                   gid = gid,
                                   sid = sid,
                                   g_type = g_type,
                                   g_QTLtype = g_QTLtype,
                                   region_tag = region_tag,
                                   include_cs_index = FALSE)

        susie_res.core.list[[region_tag]] <- susie_res_df
      }

      susie_res.core <- do.call(rbind, susie_res.core.list)
      susie_res.core
    }

    # update estimated group_prior from the current iteration
    pi_prior <- sapply(names(pi_prior), function(x){mean(EM_susie_res$susie_pip[EM_susie_res$type==x])})
    group_prior_rec[names(pi_prior),iter] <- pi_prior

    loginfo("After iteration %d, priors {%s}: {%s}", iter, names(pi_prior), pi_prior)

    # update estimated group_prior_var from the current iteration
    # in susie, mu2 is under standardized scale (if performed)
    # e.g. X = matrix(rnorm(n*p),nrow=n,ncol=p)
    # y = X %*% beta + rnorm(n)
    # res = susie(X,y,L=10)
    # X2 = 5*X
    # res2 = susie(X2,y,L=10)
    # res$mu2 is identical to res2$mu2 but coefficients are on diff scale.

    if (group_prior_var_structure=="independent"){
      V_prior <- sapply(names(V_prior),
                        function(x){
                          tmp_EM_susie_res <- EM_susie_res[EM_susie_res$type==x,];
                          sum(tmp_EM_susie_res$susie_pip*tmp_EM_susie_res$mu2)/sum(tmp_EM_susie_res$susie_pip)})
    } else if (group_prior_var_structure=="shared_all"){
      tmp_EM_susie_res <- EM_susie_res[EM_susie_res$type=="SNP",]
      V_prior["SNP"] <- sum(tmp_EM_susie_res$susie_pip*tmp_EM_susie_res$mu2)/sum(tmp_EM_susie_res$susie_pip)
      tmp_EM_susie_res <- EM_susie_res[EM_susie_res$type!="SNP",]
      V_prior[names(V_prior)!="SNP"] <- sum(tmp_EM_susie_res$susie_pip*tmp_EM_susie_res$mu2)/sum(tmp_EM_susie_res$susie_pip)
    } else if (group_prior_var_structure=="shared_QTLtype"){
      tmp_EM_susie_res <- EM_susie_res[EM_susie_res$QTLtype=="SNP",]
      V_prior["SNP"] <- sum(tmp_EM_susie_res$susie_pip*tmp_EM_susie_res$mu2)/sum(tmp_EM_susie_res$susie_pip)
      for(i in QTLtypes){
        if (i != "SNP"){
          tmp_EM_susie_res <- EM_susie_res[EM_susie_res$QTLtype==i,]
          V_prior[sapply(names(V_prior), function(x){
            unlist(strsplit(x, "[_]"))[2]})==i] <-
            sum(tmp_EM_susie_res$susie_pip*tmp_EM_susie_res$mu2)/sum(tmp_EM_susie_res$susie_pip)
        }
      }
    }
    group_prior_var_rec[names(V_prior), iter] <- V_prior
  }
  parallel::stopCluster(cl)

  group_size <- table(EM_susie_res$type)

  return(list("group_prior"= pi_prior,
              "group_prior_var" = V_prior,
              "group_prior_rec" = group_prior_rec,
              "group_prior_var_rec" = group_prior_var_rec,
              "group_prior_var_structure" = group_prior_var_structure,
              "group_size" = group_size))
}




