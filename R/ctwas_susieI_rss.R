#' SuSiE iteraitve rss version
#'
#' @importFrom logging loginfo
#' @importFrom foreach %dopar% foreach
#' @export
susieI_rss <- function(zdf,
                       regionlist,
                       ld_exprfs,
                       ld_pgenfs = NULL,
                       ld_Rfs = NULL,
                       niter = 20,
                       L= 1,
                       z_ld_weight = 0,
                       group_prior = NULL,
                       group_prior_var = NULL,
                       estimate_group_prior = T,
                       estimate_group_prior_var = T,
                       use_null_weight = T,
                       coverage = 0.95,
                       ncore = 1,
                       outputdir = getwd(),
                       outname = NULL
                      ) {

  outname <- file.path(outputdir, outname)

  ld_exprvarfs <- sapply(ld_exprfs, prep_exprvar)

  if (is.null(ld_pgenfs) & is.null(ld_Rfs)){
    stop("Error: need to provide either .pgen file or ld_R file")
  }

  if (!is.null(ld_pgenfs)){
    ld_pvarfs <- sapply(ld_pgenfs, prep_pvar, outputdir = outputdir)
  }

  K <- 2

  group_prior_rec <- matrix(, nrow = K , ncol =  niter)
  group_prior_var_rec <- matrix(, nrow = K , ncol =  niter)

  prior.gene <- group_prior[1]
  prior.SNP <-  group_prior[2]

  V.gene <- group_prior_var[1]
  V.SNP <- group_prior_var[2]

  for (iter in 1:niter){

    loginfo("run iteration %s", iter)

    snp.rpiplist <- list()
    gene.rpiplist <- list()

    cl <- parallel::makeCluster(ncore, outfile = "")
    doParallel::registerDoParallel(cl)

    corelist <- region2core(regionlist, ncore)

    outdf <- foreach (core = 1:length(corelist), .combine = "rbind",
                      .packages = "ctwas") %dopar% {
      # outdf <- NULL
      # for (core in 2:10) {

          outdf.core.list <- list()

          # run susie for each region
          regs <- corelist[[core]]
          for (reg in 1: nrow(regs)) {
          # for (reg in 1:2){
            b <- regs[reg, "b"]
            rn <- regs[reg, "rn"]

            gidx <- regionlist[[b]][[rn]][["gidx"]]
            sidx <- regionlist[[b]][[rn]][["sidx"]]
            gid <- regionlist[[b]][[rn]][["gid"]]
            sid <- regionlist[[b]][[rn]][["sid"]]

            p <- length(gidx) + length(sidx)

            if (is.null(prior.gene) | is.null(prior.SNP)){
              prior <- c(rep(1/p, length(gidx)),
                         rep(1/p, length(sidx)))
            } else {
              prior <- c(rep(prior.gene, length(gidx)),
                         rep(prior.SNP, length(sidx)))
            }

            if (is.null(V.gene) | is.null(V.SNP)){
              V <- matrix(rep(50, L * p), nrow = L)
              # following the default in susieR::susie_rss
            } else{
              V <- c(rep(V.gene, length(gidx)), rep(V.SNP, length(sidx)))
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

              if (is.null(regionlist[[b]][[rn]][["R_s_file"]])){
                R_snp <- lapply(regRDS, readRDS)
                R_snp <- as.matrix(Matrix::bdiag(R_snp))
                R_snp <- R_snp[sidx, sidx, drop = F]
              } else {
                R_snp <- readRDS(regionlist[[b]][[rn]][["R_s_file"]])
              }
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
                                  coverage = coverage)

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
         file = paste0(outname, ".susieIrssres.Rd"))
    data.table::fwrite(outdf, file= paste0(outname, ".susieIrss.txt"),
                       sep="\t", quote = F)

    parallel::stopCluster(cl)
  }


  list("group_prior"= c(prior.gene, prior.SNP),
       "group_prior_var" = c(V.gene, V.SNP))

}


