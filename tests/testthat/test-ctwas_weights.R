
test_that("load_predictdb_weights works", {

  weight_path <- system.file("extdata/sample_data", "expression_Liver.db", package = "ctwas")
  loaded_weights <- readRDS("expression_Liver.loaded.predictdb.weights.RDS")

  capture.output({
    res <- load_predictdb_weights(weight_path,
                                  filter_protein_coding_genes = TRUE,
                                  load_predictdb_LD = TRUE)
  })

  expect_equal(res, loaded_weights)
})

test_that("load_weights works", {

  weight_path <- system.file("extdata/sample_data", "expression_Liver.db", package = "ctwas")
  loaded_weights <- readRDS("expression_Liver.loaded.predictdb.weights.RDS")

  capture.output({
    res <- load_weights(weight_path,
                        weight_format = "PredictDB",
                        filter_protein_coding_genes = TRUE,
                        load_predictdb_LD = TRUE)
  })

  expect_equal(res, loaded_weights)
})


# test_that("load_fusion_weights works", {
#
#   weight_path <- "/project/mstephens/causalTWAS/apa_models/Heart_Atrial_Appendage/Heart_Atrial_Appendage/"
#   loaded_weights <- readRDS("Heart_Atrial_Appendage.loaded.fusion.weights.RDS")
#
#   capture.output({
#     res <- load_fusion_weights(weight_path, fusion_method = "enet")
#   })
#
#   expect_equal(res$weight_table, loaded_weights$weight_table)
# })
