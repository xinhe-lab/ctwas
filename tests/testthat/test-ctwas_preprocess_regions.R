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

  # Load region info
  expected_region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  # Load SNP info
  expected_snp_info <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_info.RDS", package = "ctwas"))

  expect_equal(region_info, expected_region_info)
  expect_equal(snp_info, expected_snp_info)

})
