test_that("load_weights works for predictDB weights", {

  weight_file <- system.file("extdata/sample_data", "expression_Liver.db", package = "ctwas")
  expected_weights <- readRDS("expression_Liver.loaded.predictdb.weights.RDS")

  capture.output({
    loaded_weights <- load_weights(weight_file,
                                   weight_format = "PredictDB")
  })

  expect_equal(loaded_weights, expected_weights)

})

test_that("load_weights works for FUSION weights", {

  weight_dir <- "/project/mstephens/causalTWAS/apa_models/Heart_Atrial_Appendage/Heart_Atrial_Appendage/"
  skip_if_no_weight(weight_dir)

  expected_weights <- readRDS("Heart_Atrial_Appendage.loaded.fusion.weights.RDS")

  capture.output({
    loaded_weights <- load_weights(weight_dir,
                                   weight_format = "FUSION",
                                   fusion_method = "enet",
                                   fusion_genome_version = "b38")
  })

  expect_equal(loaded_weights, expected_weights)

})
