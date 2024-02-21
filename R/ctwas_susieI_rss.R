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
#' @param niter the number of iterations of the E-M algorithm to perform
#'
#' @param L the number of effects for susie
#'
#' @param group_prior a vector of two prior inclusion probabilities for SNPs and genes. This is ignored
#' if \code{estimate_group_prior = T}
#'
#' @param group_prior_var a vector of two prior variances for SNPs and gene effects. This is ignored
#' if \code{estimate_group_prior_var = T}
#'
#' @param group_prior_var_structure a string indicating the structure to put on the prior variance parameters.
#' "independent" is the default and allows all groups to have their own separate variance parameters.
#' "shared_all" allows all groups to share the same variance parameter.
#' "shared+snps" allows all groups to share the same variance parameter, and this variance parameter is also shared with SNPs.
#' "inv_gamma" places an inverse-gamma prior on the variance parameters for each group, with shape and rate hypeparameters.
#' "shared_type" allows all groups in one molecular QTL type to share the same variance parameter.
#'
#' @param estimate_group_prior TRUE/FALSE. If TRUE, the prior inclusion probabilities for SNPs and genes are estimated
#' using the data. If FALSE, \code{group_prior} must be specified
#'
#' @param estimate_group_prior_var TRUE/FALSE. If TRUE, the prior variances for SNPs and genes are estimated
#' using the data. If FALSE, \code{group_prior_var} must be specified
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
#' @param outputdir a string, the directory to store output
#'
#' @param outname a string, the output name
#'
#' @param report_parameters TRUE/FALSE. If TRUE, estimated parameters are reported at the end of iteration
#'
#' @importFrom logging loginfo
#' @importFrom foreach %dopar% foreach
#'
#' @return a list of susie result,including PIPs and credible sets, and parameters
#'
#' @export
#'
ctwas_susieI_rss <- function(zdf,
                             regionlist,
                             region_info,
                             niter = 20,
                             L = 1,
                             group_prior = NULL,
                             group_prior_var = NULL,
                             group_prior_var_structure = c("independent","shared_all","shared+snps","inv_gamma","shared_QTLtype"),
                             estimate_group_prior = T,
                             estimate_group_prior_var = T,
                             use_null_weight = T,
                             coverage = 0.95,
                             min_abs_corr = 0.5,
                             ncore = 1,
                             verbose = T){

  types <- unique(zdf$type)
  QTLtypes <- unique(zdf$QTLtype)
  K <- length(types)

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

    if (verbose){
      loginfo("run iteration %s", iter)
    }

    cl <- parallel::makeCluster(ncore, outfile = "")
    doParallel::registerDoParallel(cl)

    corelist <- region2core(regionlist, ncore)

    susieI_res <- foreach (core = 1:length(corelist), .combine = "rbind", .packages = "ctwas") %dopar% {
      susieI_res.core.list <- list()

      # run susie for each region
      regs <- corelist[[core]]
      for (reg in 1: nrow(regs)) {
        b <- regs[reg, "b"]
        rn <- regs[reg, "rn"]

        region_data <- extract_region_data(zdf,
                                           regionlist[[b]][[rn]],
                                           b = b,
                                           rn = rn,
                                           pi_prior = pi_prior,
                                           V_prior = V_prior,
                                           use_null_weight = use_null_weight)

        z <- region_data$z
        R <- region_data$R
        prior <- region_data$prior
        V <- region_data$V
        nw <- region_data$nw

        # in susie, prior_variance is under standardized scale (if performed)
        susie_region_res <- ctwas_susie_rss(z = z,
                                            R = R,
                                            region_data = region_data,
                                            prior_weights = prior,
                                            prior_variance = V,
                                            L = L,
                                            null_weight = nw,
                                            coverage = coverage,
                                            min_abs_corr = min_abs_corr)

        susieI_res.core.list[[reg]] <- susie_region_res
      }

      susieI_res.core <- do.call(rbind, susieI_res.core.list)
      susieI_res.core
    }

    # update estimated group_prior from the current iteration
    if (isTRUE(estimate_group_prior)){
      pi_prior <- sapply(names(pi_prior), function(x){mean(susieI_res$susie_pip[susieI_res$type==x])})
      group_prior_rec[names(pi_prior),iter] <- pi_prior
    }

    if (verbose){
      loginfo("After iteration %s, priors {%s}: {%s}",
              iter, names(pi_prior), pi_prior)
    }

    # update estimated group_prior_var from the current iteration
    if (isTRUE(estimate_group_prior_var)){
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
                            temp_susieI_res <- susieI_res[susieI_res$type==x,];
                            sum(temp_susieI_res$susie_pip*temp_susieI_res$mu2)/sum(temp_susieI_res$susie_pip)})
      } else if (group_prior_var_structure=="shared_all"){
        temp_susieI_res <- susieI_res[susieI_res$type=="SNP",]
        V_prior["SNP"] <- sum(temp_susieI_res$susie_pip*temp_susieI_res$mu2)/sum(temp_susieI_res$susie_pip)

        temp_susieI_res <- susieI_res[susieI_res$type!="SNP",]
        V_prior[names(V_prior)!="SNP"] <- sum(temp_susieI_res$susie_pip*temp_susieI_res$mu2)/sum(temp_susieI_res$susie_pip)

        rm(temp_susieI_res)
      } else if (group_prior_var_structure=="shared+snps"){
        V_prior[names(V_prior)] <- sum(susieI_res$susie_pip*susieI_res$mu2)/sum(susieI_res$susie_pip)
      } else if (group_prior_var_structure=="shared_QTLtype"){
        temp_susieI_res <- susieI_res[susieI_res$QTLtype=="SNP",]
        V_prior["SNP"] <- sum(temp_susieI_res$susie_pip*temp_susieI_res$mu2)/sum(temp_susieI_res$susie_pip)
        for(i in QTLtypes){
          if (i != "SNP"){
            temp_susieI_res <- susieI_res[susieI_res$QTLtype==i,]
            V_prior[sapply(names(V_prior), function(x){
              unlist(strsplit(x, "[_]"))[2]})==i] <-
              sum(temp_susieI_res$susie_pip*temp_susieI_res$mu2)/sum(temp_susieI_res$susie_pip)
          }
        }

      }

      group_prior_var_rec[names(V_prior), iter] <- V_prior

    }

    parallel::stopCluster(cl)
  }

  param <- list("group_prior"= pi_prior,
                "group_prior_var" = V_prior,
                "group_prior_rec" = group_prior_rec,
                "group_prior_var_rec" = group_prior_var_rec,
                "group_prior_var_structure" = group_prior_var_structure,
                "L" = L,
                "coverage" = coverage,
                "min_abs_corr" = min_abs_corr)

  return(list("susieI_res" = susieI_res,
              "param" = param))
}


#' Run susie_rss for a single region
#'
#' @param z A p-vector of z scores.
#'
#' @param R A p by p symmetric, positive semidefinite correlation
#' matrix.
#'
#' @param region_data
#'
#' @param prior_weights A vector of length p, in which each entry
#'   gives the prior probability that SNP j has non-zero effect.
#'
#' @param prior_variance The prior variance. It is either a scalar or
#'   a vector of length L.
#'
#' @param L Number of components (nonzero coefficients) in the susie
#'   regression model.
#'
#' @param z_ld_weight the z_ld_weight parameter for susie_rss
#'
#' @param null_weight Prior probability of no effect (a number between
#'   0 and 1, and cannot be exactly 1).
#'
#' @param coverage A number between 0 and 1 specifying the
#'   \dQuote{coverage} of the estimated confidence sets.
#'
#' @param min_abs_corr Minimum absolute correlation allowed in a
#'   credible set. The default, 0.5, corresponds to a squared
#'   correlation of 0.25, which is a commonly used threshold for
#'   genotype data in genetic studies.
#'
#' @return a data frame of annotated susie result
#'
#' @export
#'
ctwas_susie_rss <- function(z,
                            R,
                            region_data,
                            prior_weights = NULL,
                            prior_variance = NULL,
                            L = 5,
                            z_ld_weight = 0,
                            null_weight = NULL,
                            coverage = 0.95,
                            min_abs_corr = 0.5){

  # in susie, prior_variance is under standardized scale (if performed)
  susie_res <- susie_rss(z,
                         R,
                         prior_weights = prior_weights,
                         prior_variance = prior_variance,
                         estimate_prior_variance = F,
                         L = L,
                         z_ld_weight = z_ld_weight,
                         null_weight = null_weight,
                         coverage = coverage,
                         min_abs_corr = min_abs_corr)

  # annotate susie results with SNP, gene information
  susie_region_res <- anno_susie(susie_res, region_data)

  return(susie_region_res)
}

extract_region_data <- function(zdf,
                                region_info,
                                b,
                                rn,
                                ld_exprvarfs,
                                ld_pgenfs = NULL,
                                ld_Rfs = NULL,
                                pi_prior = NULL,
                                V_prior = NULL,
                                use_null_weight = TRUE){

  gidx <- region_info[["gidx"]]
  sidx <- region_info[["sidx"]]
  gid <- region_info[["gid"]]
  sid <- region_info[["sid"]]
  g_type <- zdf$type[match(gid, zdf$id)]
  s_type <- zdf$type[match(sid, zdf$id)]
  gs_type <- c(g_type, s_type)
  g_QTLtype <- zdf$QTLtype[match(gid, zdf$id)]

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
    nw <- max(0, 1 - sum(prior))
    prior <- prior/(1-nw)
  } else {
    nw <- NULL
  }

  z.g <- zdf[match(gid, zdf$id), ][["z"]]
  z.s <- zdf[match(sid, zdf$id), ][["z"]]
  z <- c(z.g, z.s)

  if (!(is.null(ld_pgenfs))){
    # prepare LD genotype data
    # ***** DO WE STILL KEEP THIS ?
    ld_pgen <- prep_pgen(pgenf = ld_pgenfs[b], ld_pvarfs[b])

    X.g <- read_expr(ld_exprfs[b], variantidx = gidx)
    X.s <- read_pgen(ld_pgen, variantidx = sidx)
    X <- cbind(X.g, X.s)
    R <- Rfast::cora(X)
  } else {
    # prepare R matrix
    if (!("regRDS" %in% names())){
      stop("R matrix info not available for region", b, ",", rn)
    }

    regRDS <- region_info[["regRDS"]]
    R_snp <- readRDS(region_info[["R_s_file"]])
    R_snp_gene <- readRDS(region_info[["R_sg_file"]])
    R_snp_gene <- R_snp_gene[sidx, , drop = F]
    R_gene <- readRDS(region_info[["R_g_file"]])

    # gene first then SNPs
    R <- rbind(cbind(R_gene, t(R_snp_gene)),
               cbind(R_snp_gene, R_snp))
  }

  geneinfo <- read_exprvar(ld_exprvarfs[b])

  if (!is.null(ld_pgenfs)){
    snpinfo <-  read_pvar(ld_pvarfs[b])
  } else {
    snpinfo <- do.call(rbind, lapply(regRDS, read_ld_Rvar_RDS))
  }

  region_data <- list(z = z,
                      R = R,
                      geneinfo = geneinfo,
                      snpinfo = snpinfo,
                      gidx = gidx,
                      gid = gid,
                      sidx = sidx,
                      sid = sid,
                      g_type,
                      s_type,
                      gs_type,
                      g_QTLtype,
                      p = p,
                      region_tag1 = region_tag1,
                      region_tag2 = region_tag2,
                      type = type,
                      QTLtype = QTLtype,
                      prior = prior,
                      V = V,
                      nw = nw)

  return(region_data)

}


