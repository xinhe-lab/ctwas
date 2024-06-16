test_that("preprocess_z_snp works", {

  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.z_snp.RDS", package = "ctwas"))
  snp_info <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_info.RDS", package = "ctwas"))
  preprocessed_z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))

  capture.output({
    z_snp <- preprocess_z_snp(z_snp,
                              snp_info,
                              drop_multiallelic = TRUE,
                              drop_strand_ambig = TRUE)
  })

  expect_equal(z_snp, preprocessed_z_snp)

})
