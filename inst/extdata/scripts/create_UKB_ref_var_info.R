library(ctwas)

# Create variant info in b38 version
cat("Creating UKB reference variant info in b38 version ...\n")
genome_version <- "b38"
LD_R_dir <- "/project2/mstephens/wcrouse/UKB_LDR_0.1/"
outputdir <- "/project2/xinhe/shared_data/cTWAS/ref_snp_info/UKB"
dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)
region_file <- system.file("extdata/ldetect", "EUR.b38.ldetect.regions.RDS", package = "ctwas")
region_metatable <- readRDS(region_file)
filestem <- paste0("ukb_", genome_version, "_0.1")
LD_filestem <- sprintf("%s_chr%s.R_snp.%s_%s", filestem, region_metatable$chrom, region_metatable$start, region_metatable$stop)
region_metatable$LD_file <- file.path(LD_R_dir, paste0(LD_filestem, ".RDS"))
region_metatable$SNP_file <- file.path(LD_R_dir, paste0(LD_filestem, ".Rvar"))
stopifnot(all(file.exists(region_metatable$LD_file)))
stopifnot(all(file.exists(region_metatable$SNP_file)))
write.table(region_metatable, file.path(outputdir, paste0(filestem, "_region_metatable.txt")),
            quote = F, col.names = T, row.names = F, sep = "\t")

ref_snp_info <- read_snp_info_files(region_metatable$SNP_file)
write.table(ref_snp_info, gzfile(file.path(outputdir, paste0(filestem, "_var_info.Rvar.gz"))),
            quote = F, col.names = T, row.names = F, sep = "\t")


# Create reference variant info in b37 version
cat("Creating UKB reference variant info in b37 version ...\n")
genome_version <- "b37"
LD_R_dir <- "/project2/mstephens/wcrouse/UKB_LDR_0.1_b37/"
outputdir <- "/project2/xinhe/shared_data/cTWAS/ref_snp_info/UKB"
dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)
region_file <- system.file("extdata/ldetect", "EUR.b37.ldetect.regions.RDS", package = "ctwas")
region_metatable <- readRDS(region_file)
filestem <- paste0("ukb_", genome_version, "_0.1")
LD_filestem <- sprintf("%s_chr%s.R_snp.%s_%s", filestem, region_metatable$chrom, region_metatable$start, region_metatable$stop)
region_metatable$LD_file <- file.path(LD_R_dir, paste0(LD_filestem, ".RDS"))
region_metatable$SNP_file <- file.path(LD_R_dir, paste0(LD_filestem, ".Rvar"))
stopifnot(all(file.exists(region_metatable$LD_file)))
stopifnot(all(file.exists(region_metatable$SNP_file)))
write.table(region_metatable, file.path(outputdir, paste0(filestem, "_region_metatable.txt")),
            quote = F, col.names = T, row.names = F, sep = "\t")

ref_snp_info <- read_snp_info_files(region_metatable$SNP_file)
write.table(ref_snp_info, gzfile(file.path(outputdir, paste0(filestem, "_var_info.Rvar.gz"))),
            quote = F, col.names = T, row.names = F, sep = "\t")
