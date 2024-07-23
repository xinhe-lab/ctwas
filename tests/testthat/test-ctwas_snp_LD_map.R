test_that("create_snp_map works", {

  region_file <- system.file("extdata/ldetect", "EUR.b38.ldetect.regions.RDS", package = "ctwas")
  region_info <- readRDS(region_file)
  example_chrom <- 16
  region_info <- subset(region_info, chrom == example_chrom)

  ref_snp_info_file <- system.file("extdata/sample_data", "ukb_b38_0.1_chr16_var_info.Rvar.gz", package = "ctwas")
  ref_snp_info <- read.table(gzfile(ref_snp_info_file), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  capture.output({
    res <- create_snp_map(region_info, ref_snp_info)
  })
  region_info <- res$region_info
  snp_map <- res$snp_map

  expected_region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  expected_snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))

  expect_equal(region_info, expected_region_info)
  expect_equal(snp_map, expected_snp_map)

})

test_that("create_snp_LD_map works", {

  region_file <- system.file("extdata/ldetect", "EUR.b38.ldetect.regions.RDS", package = "ctwas")
  region_metatable <- readRDS(region_file)
  example_chrom <- 16
  region_metatable <- subset(region_metatable, chrom == example_chrom)

  ld_R_dir <- "/project2/mstephens/wcrouse/UKB_LDR_0.1"
  filestem <- paste0("ukb_b38_0.1")
  ld_filestem <- sprintf("%s_chr%s.R_snp.%s_%s",
                         filestem, region_metatable$chrom, region_metatable$start, region_metatable$stop)
  region_metatable$LD_file <- file.path(ld_R_dir, paste0(ld_filestem, ".RDS"))
  region_metatable$SNP_file <- file.path(ld_R_dir, paste0(ld_filestem, ".Rvar"))
  skip_if_no_LD_file(region_metatable$LD_file)

  capture.output({
    res <- create_snp_LD_map(region_metatable)
  })
  region_info <- res$region_info
  snp_map <- res$snp_map
  LD_map <- res$LD_map

  expected_region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  expected_snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  expected_LD_map <- readRDS(system.file("extdata/sample_data", "LDL_example.LD_map.RDS", package = "ctwas"))

  expect_equal(region_info, expected_region_info)
  expect_equal(snp_map, expected_snp_map)
  expect_equal(LD_map, expected_LD_map)

})
