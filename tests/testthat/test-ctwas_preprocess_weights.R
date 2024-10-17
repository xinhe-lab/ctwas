test_that("preprocess_weights with PredictDB weights works", {

  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))

  expected_weights <- readRDS("LDL_example.preprocessed.liver.weights.RDS")

  weight_file <- system.file("extdata/sample_data", "expression_Liver.db", package = "ctwas")

  capture.output({
      weights <- preprocess_weights(weight_file,
                                    region_info,
                                    z_snp$id,
                                    snp_map,
                                    type = "expression",
                                    context = "liver",
                                    weight_format = "PredictDB",
                                    drop_strand_ambig = TRUE,
                                    scale_predictdb_weights = TRUE,
                                    load_predictdb_LD = TRUE,
                                    filter_protein_coding_genes = TRUE,
                                    varID_converter_fun = convert_to_ukb_varIDs)
  })

  expect_equal(weights, expected_weights)

})

test_that("preprocess_weights with multiple PredictDB weights works", {

  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))

  expected_weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))

  capture.output({
    weight_liver_file <- system.file("extdata/sample_data", "expression_Liver.db", package = "ctwas")
    weights_liver <- preprocess_weights(weight_liver_file,
                                        region_info,
                                        z_snp$id,
                                        snp_map,
                                        type = "expression",
                                        context = "liver",
                                        weight_name = "liver_expression",
                                        varID_converter_fun = convert_to_ukb_varIDs)

    weight_adipose_file <- system.file("extdata/sample_data", "expression_Adipose_Subcutaneous.db", package = "ctwas")
    weights_adipose <- preprocess_weights(weight_adipose_file,
                                          region_info,
                                          z_snp$id,
                                          snp_map,
                                          type = "expression",
                                          context = "adipose",
                                          weight_name = "adipose_expression",
                                          varID_converter_fun = convert_to_ukb_varIDs)
    weights <- c(weights_liver, weights_adipose)
  })

  expect_equal(weights, expected_weights)

})

test_that("preprocess_weights with FUSION weights works", {

  weight_dir <- "/project/mstephens/causalTWAS/apa_models/Heart_Atrial_Appendage/Heart_Atrial_Appendage/"
  skip_if_no_weight(weight_dir)

  LD_map <- readRDS(system.file("extdata/sample_data", "LDL_example.LD_map.RDS", package = "ctwas"))
  skip_if_no_LD_file(LD_map$LD_file)

  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))

  expected_weights <- readRDS("APA.Heart_Atrial_Appendage.preprocessed.fusion.weights.RDS")

  capture.output({
    weights <- preprocess_weights(weight_dir,
                                  region_info,
                                  z_snp$id,
                                  snp_map,
                                  LD_map,
                                  type = "APA",
                                  context = "Heart_Atrial_Appendage",
                                  weight_name = "Heart_Atrial_Appendage_APA",
                                  weight_format = "FUSION",
                                  fusion_method = "enet",
                                  fusion_genome_version = "b38",
                                  drop_strand_ambig = TRUE,
                                  scale_predictdb_weights = FALSE,
                                  load_predictdb_LD = FALSE,
                                  filter_protein_coding_genes = FALSE,
                                  varID_converter_fun = convert_to_ukb_varIDs,
                                  ncore = 2)
  })

  expect_equal(weights, expected_weights)

})
