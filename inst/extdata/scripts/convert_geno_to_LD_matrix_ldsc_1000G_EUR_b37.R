
# load packages
library(ctwas)
library(Rfast)
# note: Rfast package is required for running the convert_geno_to_LD_matrix() function

# PLINK binary genotype data in .pgen or .bed format
ldref_dir <- "./1000G_EUR_Phase3_plink"
genotype_files <- file.path(ldref_dir, paste0("1000G.EUR.QC.", 1:22, ".bed"))

# PLINK reference variant information files in .pvar or .bim format.
varinfo_files <- file.path(ldref_dir, paste0("1000G.EUR.QC.", 1:22, ".bim"))

# get region_info matching the population and genome build
region_file <- system.file("extdata/ldetect", "EUR.b37.ldetect.regions.RDS", package = "ctwas")
region_info <- readRDS(region_file)

# specify output
outputdir <- "./ctwas_LD/1000G_EUR_Phase3_b37/"
outname <- "1000G_EUR_Phase3_b37"

# convert genotype to LD matrix
region_metatable <- convert_geno_to_LD_matrix(region_info,
                                              genotype_files,
                                              varinfo_files,
                                              chrom = 1:22,
                                              outputdir = outputdir,
                                              outname = outname)
