
##### Settings #####
# require pgenlibr and Rfast packages installed
library(ctwas)

# specify LD reference
ldref_dir <- "/gpfs/data/xhe-lab/ukb_LDR/genotype_data_0.1"
genotype_files <- file.path(ldref_dir, paste0("ukb_chr", 1:22, ".pgen"))

# the output Rvar files use the positions and allele information in varinfo_files
varinfo_files <- "/gpfs/data/xhe-lab/ukb_LDR/neale_lab/neale_variants.bim"

# prepare a data frame region_info for LD regions with columns "chr", "start", and "stop"
# the positions should match those in varinfo_files
ldref_population = "EUR"
genome_version = "b37"
region_file <- system.file("extdata/ldetect", paste0(ldref_population, ".", genome_version, ".bed"), package = "ctwas")
region_info <- read.table(region_file, header = T, stringsAsFactors = F)

# specify output
outputdir <- "/gpfs/data/xhe-lab/shared_data/ctwas/LDR/UKB_b37/"
outname <- "ukb_b37_0.1"

##### start #####
updated_region_info <- convert_geno_to_LD_matrix(region_info, genotype_files, varinfo_files,
                                                 chrom = 1:22,
                                                 outputdir = outputdir, outname = outname)
