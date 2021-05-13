setwd("/project2/xinhe/wcrouse/ctwas_UKBB")

####################

library(VariantAnnotation)
library(gwasvcf)
library(ctwas)

source("ctwas_LDR.R")

####################
# import GWAS VCF file and compute z-scores
# https://mrcieu.github.io/gwasvcf/articles/guide.html
# https://github.com/MRCIEU/gwas-vcf-specification

# vcffile <- "/project2/xinhe/wcrouse/GWAS_sumstats/ukb-d-30780_irnt.vcf"
# vcf <- VariantAnnotation::readVcf(vcffile)
# vcf
# 
# gwas_df <- gwasvcf::vcf_to_tibble(vcf)
# gwas_df$Z <- gwas_df$ES/gwas_df$SE
# gwas_df <- gwas_df[,c("rsid", "REF", "ALT", "Z", "SS")]
# 
# save(gwas_df, file="/project2/xinhe/wcrouse/GWAS_sumstats/ukb-d-30780_irnt.RData")
# 
# load("/project2/xinhe/wcrouse/GWAS_sumstats/ukb-d-30780_irnt.RData")
# 
# gwas_df <- gwas_df[,c("rsid", "ALT", "REF",  "Z")]
# 
# colnames(gwas_df) <- c("id", "A1", "A2", "z")
# 
# write.table(gwas_df, quote=F, sep="\t", row.names=F, file="/project2/xinhe/wcrouse/GWAS_sumstats/ukb-d-30780_irnt.sumstats")

####################

# specify GWAS summary stats
gwas_sumstat_file <- "/project2/xinhe/wcrouse/GWAS_sumstats/ukb-d-30780_irnt.sumstats"

# specify LD regions
ld_regions <- "EUR"
regionfile <- system.file("extdata", "ldetect", paste0(ld_regions, ".bed"), package = "ctwas")

# specify LD matrices
ldmat_chr <- "/project2/xinhe/wcrouse/1000G_EUR_LDMAT/LDR/1000G.EUR."
ldmat_files <- paste0(ldmat_chr, 1:22)

# specify LD reference
# ldref_chr="/project2/xinhe/shared_data/TWAS/FUSION/LDREF/1000G.EUR."
# ldref_files <- paste0(ldref_chr, 1:22, ".bed")

# specify weights; FUSION Liver
fusion_weights_dir <- "/project2/xinhe/shared_data/TWAS/FUSION/weights/GTEx_v7/"
tissue <- "Liver"
weight_dir <- paste0(fusion_weights_dir, "GTEx.", tissue, ".P01/")
weight <- paste0(weight_dir, tissue, ".P01")

# additional specifications
prediction_model <- "lasso"
out_dir <- "/project2/xinhe/wcrouse/ctwas_UKBB/output"

outname <- "ctwas_LDR"
ld_regions <- "EUR"
thin <- 1
ncore <- 1
niter2 <- 30
niter1 <- 3

####################

# STEP1: IMPUTE GENE EXPRESSION Z-SCORES

z_snp <- data.table::fread(gwas_sumstat_file)

ptm <- proc.time()
ld_exprfs <- vector()
z_gene <- NULL

for (i in 1:22){
  res <- impute_expr_z_LDR(z_snp = z_snp,
                           weight = weight,
                           method = prediction_model, outputdir = out_dir,
                           outname = paste("impute_exr_z", outname, sep = "."),
                           ldmat_file = ldmat_files[i])
  ld_exprfs[i] <- res$ld_exprf
  z_gene <- rbind(z_gene, res$z_gene)
}

timing <- proc.time() - ptm
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# save(z_gene , file=paste0(out_dir, "/", outname, "_z_gene.RData"))
# save(ld_exprfs, file=paste0(out_dir, "/", outname, "_ld_exprfs.RData"))
# 
# load(paste0(out_dir, "/", outname, "_z_gene.RData"))
# load(paste0(out_dir, "/", outname, "_ld_exprfs.RData"))

####################

# # STEP3: RUN ctwas_rss
# # ----------------------
# # Perform the causal TWAS algorithm. The algorithm will run susie iteratively for parameter estimation and
# # lastly provide PIPs for all genes and SNPs included in the analysis.

ptm <- proc.time()

pars <- ctwas_rss_LDR(z_snp = z_snp,
                      z_gene = z_gene,
                      ld_regions = ld_regions,
                      niter2 = niter2,
                      niter1 = niter1,
                      thin = thin,
                      ncore = ncore,
                      outputdir = out_dir,
                      outname = outname,
                      ldmat_files = ldmat_files,
                      ld_exprfs = ld_exprfs)

timing <- proc.time() - ptm

#save(pars, file=paste0(out_dir, "/pars.RData"))
#load(paste0(out_dir, "/pars.RData"))

cat("ctwas_rss completed. \n")
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

cat("pars: \n")
print(pars)

# R sessionInfo
sessionInfo()



