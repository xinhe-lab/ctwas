#' Iteratively run susie and estimate parameters
#' 
#' @param pgenfs A character vector of .pgen or .bed files. One file for one
#' chromosome, in the order of 1 to 22. Therefore, the length of this vector
#' needs to be 22. If .pgen files are given, then .pvar and .psam are assumed
#' to present in the same directory. If .bed files are given, then .bim and
#' .fam files are assumed to present in the same directory.
#' 
#' @param A character vector of .`expr` or `.expr.gz` files. One file for
#' one chromosome, in the order of 1 to 22. Therefore, the length of this vector
#' needs to be 22.  `.expr.gz` file is gzip compressed `.expr` files. `.expr` is
#' a matrix of imputed expression values, row is for each sample, column is for
#' each gene. Its sample order is same as in files provided by `.pgenfs`. We also
#' assume corresponding `.exprvar` files are present in the same directory.
#' `.exprvar` files are just tab delimited text files, with columns:
#' \describe{
#'   \item{chrom}{chromosome number, numeric}
#'   \item{p0}{gene boundary position, the smaller value}
#'   \item{p1}{gene boundary position, the larger value}
#'   \item{id}{gene id}
#' }
#' Its rows should be in the same order as the columns for corresponding `.expr`
#' files.
#' 
#' @param Y a vector of length n, phenotype, the same order as provided
#' by `.pgenfs` (defined in .psam or .fam files).
#' 
#' @param regionlist a list object indexing regions, variants and genes. The output of 
#' \code{index_regions}
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
#' @param standardize TRUE/FALSE. If TRUE, all variables are standardized to unit variance
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
susieI <- function(pgenfs,
                   exprfs,
                   Y,
                   regionlist,
                   niter = 20,
                   L= 1,
                   group_prior = NULL,
                   group_prior_var = NULL,
                   estimate_group_prior = T,
                   estimate_group_prior_var = T,
                   use_null_weight = T,
                   coverage = 0.95,
                   standardize = T,
                   ncore = 1,
                   outputdir = getwd(),
                   outname = NULL,
                   report_parameters=T) {

  outname <- file.path(outputdir, outname)

  pvarfs <- sapply(pgenfs, prep_pvar, outputdir = outputdir)
  exprvarfs <- sapply(exprfs, prep_exprvar)

  varY <- var(Y)

  K <- 2

  group_prior_rec <- matrix(, nrow = K , ncol =  niter)
  group_prior_var_rec <- matrix(, nrow = K , ncol =  niter)

  prior.gene <- group_prior[1]
  prior.SNP <-  group_prior[2]

  V.gene <- group_prior_var[1]
  V.SNP <- group_prior_var[2]

  cl <- parallel::makeCluster(ncore, outfile = "")
  doParallel::registerDoParallel(cl)

  for (iter in 1:niter){

    loginfo("run iteration %s", iter)

    snp.rpiplist <- list()
    gene.rpiplist <- list()

    corelist <- region2core(regionlist, ncore)

    outdf <- foreach (core = 1:length(corelist), .combine = "rbind",
                      .packages = "ctwas") %dopar% {
    # outdf <- NULL
    # for (core in 1) {

        outdf.core.list <- list()

                        # run susie for each region
        regs <- corelist[[core]]
        for (reg in 1: nrow(regs)) {
            b <- regs[reg, "b"]
            rn <- regs[reg, "rn"]

            # prepare genotype data
            pgen <- prep_pgen(pgenf = pgenfs[b], pvarfs[b])

            gidx <- regionlist[[b]][[rn]][["gidx"]]
            sidx <- regionlist[[b]][[rn]][["sidx"]]
            p <- length(gidx) + length(sidx)

            if (is.null(prior.gene) | is.null(prior.SNP)){
              prior <- c(rep(1/p, length(gidx)),
                         rep(1/p, length(sidx)))
            } else {
              prior <- c(rep(prior.gene, length(gidx)),
                         rep(prior.SNP, length(sidx)))
            }

            if (is.null(V.gene) | is.null(V.SNP)){
              V.scaled <- matrix(rep(0.2, L * p), nrow = L)
            } else{
              V.scaled <- c(rep(V.gene/varY, length(gidx)),
                            rep(V.SNP/varY, length(sidx)))
              V.scaled <- matrix(rep(V.scaled, each = L), nrow=L)
            }

            if (isTRUE(use_null_weight)){
              nw <- max(0, 1 - sum(prior))
              prior <- prior/(1-nw)
            } else {
              nw <- NULL
            }

            X.g <- read_expr(exprfs[b], variantidx = gidx)
            X.s <- read_pgen(pgen, variantidx = sidx)
            X <- cbind(X.g, X.s)

            # in susie, prior_variance is under standardized scale (if performed)
            susieres <- susie(X, Y, L = L, prior_weights = prior,
                              null_weight = nw, scaled_prior_variance = V.scaled,
                              standardize = standardize,
                              estimate_prior_variance = F, coverage = coverage)

            geneinfo <- read_exprvar(exprvarfs[b])
            snpinfo <- read_pvar(pvarfs[b])

            outdf.reg <- anno_susie(susieres,
                                   geneinfo,
                                   snpinfo,
                                   gidx,
                                   sidx,
                                   b, rn)

            outdf.core.list[[reg]] <- outdf.reg
        }

        outdf.core <- do.call(rbind, outdf.core.list)
        outdf.core
        # outdf <- rbind(outdf, outdf.core)
    }

    if (isTRUE(estimate_group_prior)){
      prior.SNP <- mean(outdf[outdf[ , "type"] == "SNP", "susie_pip"])
      prior.gene <- mean(outdf[outdf[ , "type"] == "gene", "susie_pip"])
      group_prior_rec[, iter] <- c(prior.gene, prior.SNP)
    }

    if (report_parameters){
      loginfo("After iteration %s, gene prior %s:, SNP prior:%s",
              iter, prior.gene, prior.SNP)
    }

    if (isTRUE(estimate_group_prior_var)){
      outdf.g <- outdf[outdf[ , "type"] == "gene", ]
      outdf.s <- outdf[outdf[ , "type"] == "SNP", ]

      # in susie, mu2 is under standardized scale (if performed)
      # e.g. X = matrix(rnorm(n*p),nrow=n,ncol=p)
      # y = X %*% beta + rnorm(n)
      # res = susie(X,y,L=10)
      # X2 = 5*X
      # res2 = susie(X2,y,L=10)
      # res$mu2 is identifical to res2$mu2 but coefficients are on diff scale.
      V.gene <- sum(outdf.g$susie_pip * outdf.g$mu2)/sum(outdf.g$susie_pip)
      V.SNP <- sum(outdf.s$susie_pip * outdf.s$mu2)/sum(outdf.s$susie_pip)
      group_prior_var_rec[, iter] <- c(V.gene, V.SNP)
    }

    save(group_prior_rec, group_prior_var_rec,
         file = paste0(outname, ".susieIres.Rd"))
    data.table::fwrite(outdf, file= paste0(outname, ".susieI.txt"),
                      sep="\t", quote = F)
  }

  parallel::stopCluster(cl)

  list("group_prior"= c(prior.gene, prior.SNP),
       "group_prior_var" = c(V.gene, V.SNP))

}
