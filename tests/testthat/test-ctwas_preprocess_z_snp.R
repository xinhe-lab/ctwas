test_that("preprocess_z_snp works", {

  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.z_snp.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  preprocessed_z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))

  capture.output({
    z_snp <- preprocess_z_snp(z_snp, snp_map)
  })

  expect_equal(z_snp, preprocessed_z_snp)

})
