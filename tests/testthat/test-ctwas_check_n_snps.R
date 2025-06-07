test_that("check_n_snps works", {

  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))
  weights <- readRDS("LDL_example.preprocessed.liver.weights.RDS")

  expected_count_n_snp_res <- readRDS("LDL_example.check_n_snps_res.RDS")

  capture.output({
    check_n_snps_res <- check_n_snps(snp_map, z_snp, weights)
  })

  # saveRDS(check_n_snps_res, "LDL_example.check_n_snps_res.RDS")

  expect_equal(check_n_snps_res, expected_count_n_snp_res)
})
