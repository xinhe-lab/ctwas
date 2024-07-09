test_that("preprocess_weights works", {

  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_info <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_info.RDS", package = "ctwas"))
  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))
  precomputed_weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))

  weight_path <- system.file("extdata/sample_data", "expression_Liver.db", package = "ctwas")

  capture.output({
    weights <- preprocess_weights(weight_path,
                                  region_info,
                                  z_snp$id,
                                  snp_info,
                                  weight_format = "PredictDB")
  })

  expect_equal(weights, precomputed_weights)

})
