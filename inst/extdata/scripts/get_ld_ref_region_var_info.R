
# get variant info in LD reference for b38 version
genome_version <- "b38"
ld_R_dir <- "/project2/mstephens/wcrouse/UKB_LDR_0.1/"
outputdir <- "/project2/xinhe/shared_data/multigroup_ctwas/LD_reference/"
# get region_info that match the population and genome build
region_file <- system.file("extdata/ldetect", "EUR.b38.ldetect.regions.RDS", package = "ctwas")
region_info <- readRDS(region_file)
filestem <- paste0("ukb_", genome_version, "_0.1")
ld_filestem <- sprintf("%s_chr%s.R_snp.%s_%s", filestem, region_info$chrom, region_info$start, region_info$stop)
region_info$LD_matrix <- file.path(ld_R_dir, paste0(ld_filestem, ".RDS"))
region_info$SNP_info <- file.path(ld_R_dir, paste0(ld_filestem, ".Rvar"))
stopifnot(all(file.exists(region_info$LD_matrix)))
stopifnot(all(file.exists(region_info$SNP_info)))
write.table(region_info, file.path(outputdir, paste0(filestem, "_region_info.txt")),
            quote = F, col.names = T, row.names = F, sep = "\t")

ld_snpinfo <- read_LD_SNP_files(region_info$SNP_info)
write.table(ld_snpinfo, gzfile(file.path(outputdir, paste0(filestem, "_var_info.Rvar.gz"))),
            quote = F, col.names = T, row.names = F, sep = "\t")


# get variant info in LD reference for b37 version
genome_version <- "b37"
ld_R_dir <- "/project2/mstephens/wcrouse/UKB_LDR_0.1_b37/"
outputdir <- "/project2/xinhe/shared_data/multigroup_ctwas/LD_reference/"
# get region_info that match the population and genome build
region_file <- system.file("extdata/ldetect", "EUR.b37.ldetect.regions.RDS", package = "ctwas")
region_info <- readRDS(region_file)
filestem <- paste0("ukb_", genome_version, "_0.1")
ld_filestem <- sprintf("%s_chr%s.R_snp.%s_%s", filestem, region_info$chrom, region_info$start, region_info$stop)
region_info$LD_matrix <- file.path(ld_R_dir, paste0(ld_filestem, ".RDS"))
region_info$SNP_info <- file.path(ld_R_dir, paste0(ld_filestem, ".Rvar"))
stopifnot(all(file.exists(region_info$LD_matrix)))
stopifnot(all(file.exists(region_info$SNP_info)))
write.table(region_info, file.path(outputdir, paste0(filestem, "_region_info.txt")),
            quote = F, col.names = T, row.names = F, sep = "\t")

ld_snpinfo <- read_LD_SNP_files(region_info$SNP_info)
write.table(ld_snpinfo, gzfile(file.path(outputdir, paste0(filestem, "_var_info.Rvar.gz"))),
            quote = F, col.names = T, row.names = F, sep = "\t")
