test_that("ctwas_sumstats works", {

  LD_map <- readRDS(system.file("extdata/sample_data", "LDL_example.LD_map.RDS", package = "ctwas"))
  skip_if_no_LD_file(LD_map$LD_file)

  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))
  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))

  expected_ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_res.RDS", package = "ctwas"))

  capture.output({
    ctwas_res <- ctwas_sumstats(z_snp,
                                weights,
                                region_info,
                                LD_map,
                                snp_map,
                                thin = 0.1,
                                maxSNP = 20000,
                                ncore = 2)
  })

  expect_equal(ctwas_res, expected_ctwas_res)
  # expect_equal(ctwas_res$z_gene, expected_ctwas_res$z_gene)
  # expect_equal(ctwas_res$param, expected_ctwas_res$param)
  # expect_equal(ctwas_res$finemap_res, expected_ctwas_res$finemap_res)
  # expect_equal(ctwas_res$region_data, expected_ctwas_res$region_data)
  # expect_equal(ctwas_res$boundary_genes, expected_ctwas_res$boundary_genes)
  # expect_equal(ctwas_res$screen_res, expected_ctwas_res$screen_res)

})
