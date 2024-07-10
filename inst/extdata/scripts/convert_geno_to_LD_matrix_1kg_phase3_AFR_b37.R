
# load packages
library(ctwas)
# note: Rfast package is required for running the convert_geno_to_LD_matrix() function

# specify LD reference
ldref_dir <- "/project2/xinhe/1kg/1000G_Phase3_plink/b37/1000G_AFR_phase3_b37"
genotype_files <- file.path(ldref_dir, paste0("1000G_AFR_phase3_b37_chr", 1:22, ".bed"))
# the output Rvar files will use the positions and allele information from these files
varinfo_files <- file.path(ldref_dir, paste0("1000G_AFR_phase3_b37_chr", 1:22, ".bim"))

# get region_info that match the population and genome build
region_file <- system.file("extdata/ldetect", "AFR.b37.ldetect.regions.RDS", package = "ctwas")
region_info <- readRDS(region_file)

# specify output
outputdir <- "/project2/xinhe/shared_data/ctwas/LDR/1000G_AFR_phase3_b37/"
outname <- "1000G_AFR_phase3_b37"

# convert genotype to LD matrix
updated_region_info <- convert_geno_to_LD_matrix(region_info,
                                                 genotype_files,
                                                 varinfo_files,
                                                 chrom = 1:22,
                                                 outputdir = outputdir,
                                                 outname = outname,
                                                 verbose = TRUE,
                                                 show_progress_bar = FALSE)
