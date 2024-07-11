test_that("preprocess_region_LD_snp_info (no-LD version) works", {

  region_file <- system.file("extdata/ldetect", "EUR.b38.ldetect.regions.RDS", package = "ctwas")
  region_info <- readRDS(region_file)

  example_chrom <- 16
  ref_snp_info_file <- system.file("extdata/sample_data", "ukb_b38_0.1_chr16_var_info.Rvar.gz", package = "ctwas")
  ref_snp_info <- read.table(gzfile(ref_snp_info_file), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  capture.output({
    res <- preprocess_region_LD_snp_info(region_info,
                                         ref_snp_info = ref_snp_info,
                                         chrom = example_chrom,
                                         use_LD = FALSE)
  })
  region_info <- res$region_info
  snp_info <- res$snp_info

  expected_region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  expected_snp_info <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_info.RDS", package = "ctwas"))

  expect_equal(region_info, expected_region_info)
  expect_equal(snp_info, expected_snp_info)

})

test_that("preprocess_region_LD_snp_info (LD version) works", {

  region_file <- system.file("extdata/ldetect", "EUR.b38.ldetect.regions.RDS", package = "ctwas")
  region_info <- readRDS(region_file)

  example_chrom <- 16

  ld_R_dir <- "/project2/mstephens/wcrouse/UKB_LDR_0.1"
  filestem <- paste0("ukb_b38_0.1")
  ld_filestem <- sprintf("%s_chr%s.R_snp.%s_%s", filestem, region_info$chrom, region_info$start, region_info$stop)
  region_info$LD_matrix <- file.path(ld_R_dir, paste0(ld_filestem, ".RDS"))
  region_info$SNP_info <- file.path(ld_R_dir, paste0(ld_filestem, ".Rvar"))
  skip_if_no_LD_matrix(region_info$LD_matrix)

  capture.output({
    res <- preprocess_region_LD_snp_info(region_info,
                                         chrom = example_chrom,
                                         use_LD = TRUE)
  })
  region_info <- res$region_info
  snp_info <- res$snp_info
  LD_info <- res$LD_info

  expected_region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  expected_snp_info <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_info.RDS", package = "ctwas"))
  expected_LD_info <- readRDS(system.file("extdata/sample_data", "LDL_example.LD_info.RDS", package = "ctwas"))

  expect_equal(region_info, expected_region_info)
  expect_equal(snp_info, expected_snp_info)
  expect_equal(LD_info, expected_LD_info)

})
