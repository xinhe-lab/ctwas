test_that("preprocess_z_snp works", {

  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.z_snp.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  expected_preprocessed_z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))

  capture.output({
    preprocessed_z_snp <- preprocess_z_snp(z_snp, snp_map,
                                           drop_multiallelic = TRUE,
                                           drop_strand_ambig = TRUE,
                                           varID_converter_fun = convert_to_ukb_varIDs)
  })

  expect_equal(preprocessed_z_snp, expected_preprocessed_z_snp)

})
