test_that("ctwas_sumstats_noLD works", {

  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))
  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))

  precomputed_ctwas_res <- readRDS("LDL_example.ctwas_sumstats_noLD_res.RDS")

  capture.output({
    ctwas_res <- ctwas_sumstats_noLD(z_snp,
                                     weights,
                                     region_info,
                                     snp_map,
                                     thin = 0.1,
                                     maxSNP = 20000,
                                     min_nonSNP_PIP = 0.5,
                                     ncore = 6)
  })

  expect_equal(ctwas_res, precomputed_ctwas_res)

})
