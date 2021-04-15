## code to prepare `Y` dataset goes here
library(ctwas)

pgenfs <- system.file("extdata/example_genotype_files", paste0("example_chr", 1:22, ".pgen"), package = "ctwas")
weight <- system.file("extdata/example_fusion_weights", "Tissue", package = "ctwas")

outputdir <- "~"
exprfs <- vector()
for (i in 1:22){
  pgenf <- pgenfs[i]
  exprfs[i] <- impute_expr(pgenf = pgenf, weight = weight,
                           method = "lasso", outputdir = outputdir,
                           outname = "test")
}

# one causal gene
X.g <- read_expr(exprfs[1])
X.g <- scale(X.g)

# 10 causal snps
X.s.list <- list()
for (i in 1:10){
  pgen <- prep_pgen(pgenfs[i], prep_pvar(pgenfs[i]))
  X.s.list[[i]] <- read_pgen(pgen, variantidx = 200)
}
X.s <- do.call(cbind, X.s.list)

X.s <- scale(X.s)
nsample <- nrow(X.g)

# effect size
beta.gene <- 0.5
beta.snp <- rep(0.2, 10)

# get phenotype
set.seed(100)
Y <- X.g %*% beta.gene + X.s%*% beta.snp + rnorm(nsample, 0, 1)

# get summary statistics
z_snp <- NULL
for (i in 1:22){
  pvar <- prep_pvar(pgenfs[i])
  pgen <- prep_pgen(pgenfs[i], pvar)
  snp <- read_pvar(pvar)
  X <- read_pgen(pgen)
  ss <- univariate_regression(X, Y)
  zhat <- with(ss, betahat/sebetahat)
  z_snp <- rbind(z_snp, cbind(snp[, c("id", "alt", "ref")], zhat))
}
colnames(z_snp) <- c("id", "A1", "A2", "z")

usethis::use_data(Y, z_snp, overwrite = T)

