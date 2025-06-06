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
