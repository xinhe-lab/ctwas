test_that("assemble_region_data works", {

  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))
  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))
  z_gene <- readRDS(system.file("extdata/sample_data", "LDL_example.z_gene.RDS", package = "ctwas"))

  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_noLD_res.RDS", package = "ctwas"))
  preassembled_region_data <- ctwas_res$region_data
  preassembled_boundary_genes <- ctwas_res$boundary_genes

  capture.output({
    res <- assemble_region_data(region_info,
                                z_snp,
                                z_gene,
                                weights,
                                snp_map,
                                thin = 0.1,
                                maxSNP = 20000)
  })
  region_data <- res$region_data
  boundary_genes <- res$boundary_genes

  expect_equal(region_data, preassembled_region_data)
  expect_equal(boundary_genes, preassembled_boundary_genes)

})

test_that("expand_region_data works", {

  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))
  z_gene <- readRDS(system.file("extdata/sample_data", "LDL_example.z_gene.RDS", package = "ctwas"))

  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_res.RDS", package = "ctwas"))
  preassembled_region_data <- ctwas_res$region_data
  expected_screened_region_data <- ctwas_res$screened_region_data
  expected_screened_region_ids <- names(expected_screened_region_data)
  preassembled_region_data <- preassembled_region_data[expected_screened_region_ids]

  capture.output({
    screened_region_data <- expand_region_data(preassembled_region_data,
                                               snp_map,
                                               z_snp,
                                               trim_by = "z",
                                               maxSNP = 20000,
                                               ncore = 1)
  })

  expect_equal(screened_region_data, expected_screened_region_data)

})
