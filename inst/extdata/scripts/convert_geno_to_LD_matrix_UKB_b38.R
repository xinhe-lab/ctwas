
# load packages
library(ctwas)
library(Rfast)
# note: Rfast package is required for running the convert_geno_to_LD_matrix() function

# PLINK binary genotype data in .pgen or .bed format
ldref_dir <- "/gpfs/data/xhe-lab/ukb_LDR/genotype_data_0.1"
genotype_files <- file.path(ldref_dir, paste0("ukb_chr", 1:22, ".pgen"))

# PLINK reference variant information files in .pvar or .bim format.
varinfo_files <- "/gpfs/data/xhe-lab/ukb_LDR/neale_lab/neale_variants_hg38.bim"

# get region_info that match the population and genome build
region_file <- system.file("extdata/ldetect", "EUR.b38.ldetect.regions.RDS", package = "ctwas")
region_info <- readRDS(region_file)

# specify output
outputdir <- "/gpfs/data/xhe-lab/shared_data/ctwas/LDR/UKB_b38/"
outname <- "ukb_b38_0.1"

# convert genotype to LD matrix
region_metatable <- convert_geno_to_LD_matrix(region_info,
                                              genotype_files,
                                              varinfo_files,
                                              chrom = 1:22,
                                              outputdir = outputdir,
                                              outname = outname)
