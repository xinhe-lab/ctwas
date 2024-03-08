#' Run EM to estimate parameters.
#' Iteratively run susie and estimate parameters - RSS version
#'
#' @param zdf A data frame with three columns: "id", "z", "type".
#' This data frame gives the the z scores for SNPs and genes, denoted in "type".
#' The "type" column can also be use to specify multiple sets of weights
#'
#' @param regionlist a list object indexing regions, variants and genes.
#'
#' @param region_info a data frame of region definition and associated file names
#'
#' @param gene_info a data frame of gene information
#'
#' @param niter the number of iterations of the E-M algorithm to perform
#'
#' @param group_prior a vector of two prior inclusion probabilities for SNPs and genes.
#'
#' @param group_prior_var a vector of two prior variances for SNPs and gene effects.
#'
#' @param group_prior_var_structure a string indicating the structure to put on the prior variance parameters.
#' "independent" is the default and allows all groups to have their own separate variance parameters.
#' "shared_all" allows all groups to share the same variance parameter.
#' "shared+snps" allows all groups to share the same variance parameter, and this variance parameter is also shared with SNPs.
#' "inv_gamma" places an inverse-gamma prior on the variance parameters for each group, with shape and rate hypeparameters.
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
#' @param ncore The number of cores used to parallelize susie over regions
#'
#' @importFrom logging loginfo
#' @importFrom foreach %dopar% foreach
#'
#' @return a list of parameters
#'
#' @export
#'
ctwas_EM <- function(zdf,
                     regionlist,
                     region_info,
                     gene_info,
                     niter = 20,
                     group_prior = NULL,
                     group_prior_var = NULL,
                     group_prior_var_structure = c("independent","shared_all","shared+snps","shared_QTLtype"),
                     use_null_weight = TRUE,
                     coverage = 0.95,
                     min_abs_corr = 0.5,
                     max_iter = 1,
                     ncore = 1,
                     verbose = TRUE){

  types <- unique(zdf$type)
  QTLtypes <- unique(zdf$QTLtype)
  K <- length(types)

  # store group priors from each iteration
  group_prior_rec <- matrix(NA, nrow = K , ncol =  niter)
  group_prior_var_rec <- matrix(NA, nrow = K , ncol =  niter)
  rownames(group_prior_rec) <- types
  rownames(group_prior_var_rec) <- types

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

  for (iter in 1:niter){

    if (verbose)
      loginfo("run iteration %s", iter)

    cl <- parallel::makeCluster(ncore, outfile = "")
    doParallel::registerDoParallel(cl)

    corelist <- region2core(regionlist, ncore)

    susie_res_regions <- foreach (core = 1:length(corelist), .combine = "rbind", .packages = "ctwas") %dopar% {
      susie_res.core.list <- list()

      # run susie for each region
      regs <- corelist[[core]]
      for (reg in 1: nrow(regs)) {
        b <- regs[reg, "b"]
        rn <- regs[reg, "rn"]

        region_idx <- regionlist[[b]][[rn]]
        gidx <- region_idx[["gidx"]]
        sidx <- region_idx[["sidx"]]
        gid <- region_idx[["gid"]]
        sid <- region_idx[["sid"]]
        g_type <- zdf$type[match(gid, zdf$id)]
        s_type <- zdf$type[match(sid, zdf$id)]
        gs_type <- c(g_type, s_type)
        # g_QTLtype <- zdf$QTLtype[match(gid, zdf$id)]

        p <- length(gidx) + length(sidx)

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

        z.g <- zdf[match(gid, zdf$id), ][["z"]]
        z.s <- zdf[match(sid, zdf$id), ][["z"]]
        z <- c(z.g, z.s)

        # R does not matter for susie when L = 1
        R <- diag(length(z))

        # in susie, prior_variance is under standardized scale (if performed)
        susie_res <- ctwas_susie_rss(z = z,
                                     R = R,
                                     prior_weights = prior,
                                     prior_variance = V,
                                     L = 1,
                                     null_weight = null_weight,
                                     coverage = coverage,
                                     min_abs_corr = min_abs_corr,
                                     max_iter = max_iter,
                                     ...)

        # annotate susie result with SNP and gene information
        gene_info_chr <- gene_info[gene_info$chrom == b, ]
        ld_snpinfo <- read_LD_SNP_file(region_idx[["SNP_info"]])

        susie_res <- anno_susie(susie_res,
                                gene_info = gene_info_chr,
                                snp_info = ld_snpinfo,
                                region_idx = region_idx,
                                zdf = zdf)

        susie_res.core.list[[reg]] <- susie_res
      }

      susie_res.core <- do.call(rbind, susie_res.core.list)
      susie_res.core
    }

    # update estimated group_prior from the current iteration
    pi_prior <- sapply(names(pi_prior), function(x){mean(susie_res_regions$susie_pip[susie_res_regions$type==x])})
    group_prior_rec[names(pi_prior),iter] <- pi_prior

    if (verbose){
      loginfo("After iteration %s, priors {%s}: {%s}",
              iter, names(pi_prior), pi_prior)
    }

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
                          temp_susie_res_regions <- susie_res_regions[susie_res_regions$type==x,];
                          sum(temp_susie_res_regions$susie_pip*temp_susie_res_regions$mu2)/sum(temp_susie_res_regions$susie_pip)})
    } else if (group_prior_var_structure=="shared_all"){
      temp_susie_res_regions <- susie_res_regions[susie_res_regions$type=="SNP",]
      V_prior["SNP"] <- sum(temp_susie_res_regions$susie_pip*temp_susie_res_regions$mu2)/sum(temp_susie_res_regions$susie_pip)

      temp_susie_res_regions <- susie_res_regions[susie_res_regions$type!="SNP",]
      V_prior[names(V_prior)!="SNP"] <- sum(temp_susie_res_regions$susie_pip*temp_susie_res_regions$mu2)/sum(temp_susie_res_regions$susie_pip)

    } else if (group_prior_var_structure=="shared+snps"){
      V_prior[names(V_prior)] <- sum(susie_res_regions$susie_pip*susie_res_regions$mu2)/sum(susie_res_regions$susie_pip)
    } else if (group_prior_var_structure=="shared_QTLtype"){
      temp_susie_res_regions <- susie_res_regions[susie_res_regions$QTLtype=="SNP",]
      V_prior["SNP"] <- sum(temp_susie_res_regions$susie_pip*temp_susie_res_regions$mu2)/sum(temp_susie_res_regions$susie_pip)
      for(i in QTLtypes){
        if (i != "SNP"){
          temp_susie_res_regions <- susie_res_regions[susie_res_regions$QTLtype==i,]
          V_prior[sapply(names(V_prior), function(x){
            unlist(strsplit(x, "[_]"))[2]})==i] <-
            sum(temp_susie_res_regions$susie_pip*temp_susie_res_regions$mu2)/sum(temp_susie_res_regions$susie_pip)
        }
      }

    }

    group_prior_var_rec[names(V_prior), iter] <- V_prior

    parallel::stopCluster(cl)
  }

  return(list("group_prior"= pi_prior,
              "group_prior_var" = V_prior,
              "group_prior_rec" = group_prior_rec,
              "group_prior_var_rec" = group_prior_var_rec,
              "group_prior_var_structure" = group_prior_var_structure))
}




