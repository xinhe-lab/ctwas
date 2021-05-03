#' SuSiE iteraitve
#' @param standardize True/False. Follows susie convention: If
#' standardize = TRUE, standardize the columns of X to unit variance
#'  prior to fitting. Note that scaled_prior_variance' specifies the prior
#'  on the coefficients of X after standardization (if it is performed).
#'  If you do not standardize, you may need to think more carefully about
#'  specifying scaled_prior_variance. Whatever your choice, the coefficients
#'  returned by coef are given for X on the original input scale. Any
#'  column of X that has zero variance is not standardized, but left
#'  as is.
#'
#' @importFrom logging loginfo
#' @importFrom foreach %dopar% foreach
#' @export
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
                   outname = NULL) {

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
        # for (core in 1:1) {

        outdf.core.list <- list()

        # run susie for each region
        regs <- corelist[[core]]
        for (reg in 1: nrow(regs)) {
            b <- regs[reg, "b"]
            rn <- regs[reg, "rn"]

            print(c(b,rn))

            # prepare genotype data
            pgen <- prep_pgen(pgenf = pgenfs[b], pvarfs[b])

            gidx <- regionlist[[b]][[rn]][["gidx"]]
            sidx <- regionlist[[b]][[rn]][["sidx"]]
            prop <-  regionlist[[b]][[rn]][["prop"]]
            if (is.null(prop)) prop <- 1

            p <- length(gidx) + length(sidx)

            if (is.null(prior.gene) | is.null(prior.SNP)){
              prior <- c(rep(1/p, length(gidx)),
                         rep(1/p, length(sidx)))
            } else {
              prior <- c(rep(prior.gene, length(gidx)),
                         rep(prior.SNP, length(sidx)))
            }

            if (length(sidx) >=1){
                prior[(length(gidx) + 1) : p] <- prior[(length(gidx) + 1) : p]/prop
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

            outdf.reg <- anno_susie(susieres,
                                   exprvarfs[b],
                                   pvarfs[b],
                                   gidx,
                                   sidx,
                                   b, rn)

            outdf.core.list[[reg]] <- outdf.reg
        }

        outdf.core <- do.call(rbind, outdf.core.list)
        outdf.core
    }

    if (isTRUE(estimate_group_prior)){
      prior.SNP <- mean(outdf[outdf[ , "type"] == "SNP", "susie_pip"])
      prior.gene <- mean(outdf[outdf[ , "type"] == "gene", "susie_pip"])
      group_prior_rec[, iter] <- c(prior.gene, prior.SNP)
    }

    loginfo("After iteration %s, gene prior %s:, SNP prior:%s",
                 iter, prior.gene, prior.SNP)

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


