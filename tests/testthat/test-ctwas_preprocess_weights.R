test_that("preprocess_weights works", {

  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))

  expected_weights <- readRDS("LDL_example.preprocessed.liver.weights.RDS")

  weight_dir <- system.file("extdata/sample_data", "expression_Liver.db", package = "ctwas")

  capture.output({
    weights <- preprocess_weights(weight_dir,
                                  region_info,
                                  z_snp$id,
                                  snp_map,
                                  type = "expression",
                                  context = "liver",
                                  weight_format = "PredictDB",
                                  drop_strand_ambig = TRUE,
                                  scale_predictdb_weights = TRUE,
                                  load_predictdb_LD = TRUE,
                                  filter_protein_coding_genes = TRUE)
  })

  expect_equal(weights, expected_weights)

})
