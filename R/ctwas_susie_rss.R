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
#' @param param a list of estimated parameters, including \code{group_prior}, \code{group_prior_var}, and \code{group_prior_var_structure}.
#'
#' @param niter the number of iterations of the E-M algorithm to perform
#'
#' @param L the number of effects for susie
#'
#' @param z_ld_weight the z_ld_weight parameter for susie_rss
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
ctwas_susieI_rss <- function(zdf,
                             regionlist,
                             region_info,
                             param = NULL,
                             niter = 20,
                             L= 1,
                             z_ld_weight = 0,
                             estimate_group_prior = T,
                             estimate_group_prior_var = T,
                             use_null_weight = T,
                             coverage = 0.95,
                             min_abs_corr = 0.5,
                             ncore = 1,
                             report_parameters=T){

  if (is.null(ld_pgenfs) & is.null(ld_Rfs)){
    stop("Error: need to provide either .pgen file or ld_R file")
  }

  if (!is.null(ld_pgenfs)){
    ld_pvarfs <- sapply(ld_pgenfs, prep_pvar, outputdir = outputdir)
  }

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

    loginfo("run iteration %s", iter)

    cl <- parallel::makeCluster(ncore, outfile = "")
    doParallel::registerDoParallel(cl)

    corelist <- region2core(regionlist, ncore)

    outdf <- foreach (core = 1:length(corelist), .combine = "rbind",
                      .packages = "ctwas") %dopar% {
                        outdf.core.list <- list()

                        # run susie for each region
                        regs <- corelist[[core]]
                        for (reg in 1: nrow(regs)) {
                          b <- regs[reg, "b"]
                          rn <- regs[reg, "rn"]

                          gidx <- regionlist[[b]][[rn]][["gidx"]]
                          sidx <- regionlist[[b]][[rn]][["sidx"]]
                          gid <- regionlist[[b]][[rn]][["gid"]]
                          sid <- regionlist[[b]][[rn]][["sid"]]
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
                            ld_pgen <- prep_pgen(pgenf = ld_pgenfs[b], ld_pvarfs[b])

                            X.g <- read_expr(ld_exprfs[b], variantidx = gidx)
                            X.s <- read_pgen(ld_pgen, variantidx = sidx)
                            X <- cbind(X.g, X.s)
                            R <- Rfast::cora(X)
                          } else {
                            # prepare R matrix
                            if (!("regRDS" %in% names(regionlist[[b]][[rn]]))){
                              stop("R matrix info not available for region", b, ",", rn)
                            }

                            regRDS <- regionlist[[b]][[rn]][["regRDS"]]
                            R_snp <- readRDS(regionlist[[b]][[rn]][["R_s_file"]])
                            R_snp_gene <- readRDS(regionlist[[b]][[rn]][["R_sg_file"]])
                            R_snp_gene <- R_snp_gene[sidx, , drop = F]
                            R_gene <- readRDS(regionlist[[b]][[rn]][["R_g_file"]])

                            # gene first then SNPs
                            R <- rbind(cbind(R_gene, t(R_snp_gene)),
                                       cbind(R_snp_gene, R_snp))
                          }

                          # in susie, prior_variance is under standardized scale (if performed)
                          susieres <- susie_rss(z, R,
                                                z_ld_weight = z_ld_weight,
                                                L = L, prior_weights = prior,
                                                null_weight = nw,
                                                prior_variance = V,
                                                estimate_prior_variance = F,
                                                coverage = coverage,
                                                min_abs_corr = min_abs_corr)

                          geneinfo <- read_exprvar(ld_exprvarfs[b])

                          if (!is.null(ld_pgenfs)){
                            snpinfo <-  read_pvar(ld_pvarfs[b])
                          } else {
                            snpinfo <- do.call(rbind, lapply(regRDS, read_ld_Rvar_RDS))
                          }

                          outdf.reg <- anno_susie(susieres,
                                                  geneinfo,
                                                  snpinfo,
                                                  gidx,
                                                  sidx,
                                                  b,
                                                  rn,
                                                  g_type,
                                                  g_QTLtype)

                          outdf.core.list[[reg]] <- outdf.reg
                        }

                        outdf.core <- do.call(rbind, outdf.core.list)
                        outdf.core
                      }

    if (isTRUE(estimate_group_prior)){
      pi_prior <- sapply(names(pi_prior), function(x){mean(outdf$susie_pip[outdf$type==x])})
      group_prior_rec[names(pi_prior),iter] <- pi_prior
    }

    if (report_parameters){
      loginfo("After iteration %s, priors {%s}: {%s}",
              iter, names(pi_prior), pi_prior)
    }

    if (isTRUE(estimate_group_prior_var)){
      # in susie, mu2 is under standardized scale (if performed)
      # e.g. X = matrix(rnorm(n*p),nrow=n,ncol=p)
      # y = X %*% beta + rnorm(n)
      # res = susie(X,y,L=10)
      # X2 = 5*X
      # res2 = susie(X2,y,L=10)
      # res$mu2 is identical to res2$mu2 but coefficients are on diff scale.

      if (group_prior_var_structure=="independent"){
        V_prior <- sapply(names(V_prior), function(x){outdf_temp <- outdf[outdf$type==x,]; sum(outdf_temp$susie_pip*outdf_temp$mu2)/sum(outdf_temp$susie_pip)})
      } else if (group_prior_var_structure=="shared_all"){
        outdf_temp <- outdf[outdf$type=="SNP",]
        V_prior["SNP"] <- sum(outdf_temp$susie_pip*outdf_temp$mu2)/sum(outdf_temp$susie_pip)

        outdf_temp <- outdf[outdf$type!="SNP",]
        V_prior[names(V_prior)!="SNP"] <- sum(outdf_temp$susie_pip*outdf_temp$mu2)/sum(outdf_temp$susie_pip)

        rm(outdf_temp)
      } else if (group_prior_var_structure=="shared+snps"){
        V_prior[names(V_prior)] <- sum(outdf$susie_pip*outdf$mu2)/sum(outdf$susie_pip)
      } else if (group_prior_var_structure=="shared_QTLtype"){
        outdf_temp <- outdf[outdf$QTLtype=="SNP",]
        V_prior["SNP"] <- sum(outdf_temp$susie_pip*outdf_temp$mu2)/sum(outdf_temp$susie_pip)
        for(i in QTLtypes){
          if (i != "SNP"){
            outdf_temp <- outdf[outdf$QTLtype==i,]
            V_prior[sapply(names(V_prior), function(x){unlist(strsplit(x, "[_]"))[2]})==i] <- sum(outdf_temp$susie_pip*outdf_temp$mu2)/sum(outdf_temp$susie_pip)
          }
        }
      }

      group_prior_var_rec[names(V_prior), iter] <- V_prior
    }

    parallel::stopCluster(cl)
  }

  return(list("susieIrss_res" = outdf,
              "group_prior"= pi_prior,
              "group_prior_var" = V_prior,
              "group_prior_rec" = group_prior_rec,
              "group_prior_var_rec" = group_prior_var_rec))
}
