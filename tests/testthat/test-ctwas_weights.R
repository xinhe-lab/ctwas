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

test_that("trim_weights works", {

  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))

  weights <- readRDS("LDL_example.preprocessed.liver.weights.RDS")
  expected_trimmed_weights <- readRDS("LDL_example.preprocessed.trimmed.liver.weights.RDS")

  capture.output({
    trimmed_weights <- trim_weights(weights, snp_map, top_n_snps = 2)
  })

  # saveRDS(trimmed_weights, "LDL_example.preprocessed.trimmed.liver.weights.RDS")

  expect_equal(trimmed_weights, expected_trimmed_weights)

})

test_that("subset_weights works", {

  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))
  expected_weights_subset <- readRDS("LDL_example.preprocessed.subset.liver.expression.weights.RDS")

  capture.output({
    weights_subset <- subset_weights(weights, names = "liver|expression", select_by = "group")
  })

  saveRDS(weights_subset, "LDL_example.preprocessed.subset.liver.expression.weights.RDS")

  expect_equal(weights_subset, expected_weights_subset)

})

