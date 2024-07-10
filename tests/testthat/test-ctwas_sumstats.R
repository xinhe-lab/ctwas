test_that("ctwas_sumstats works", {

  LD_info <- readRDS(system.file("extdata/sample_data", "LDL_example.LD_info.RDS", package = "ctwas"))
  skip_if_no_LD_matrix(LD_info$LD_matrix)

  # Load z_snp
  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))
  # Load weights
  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))
  # Load region_info, snp_info, and LD_info
  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_info <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_info.RDS", package = "ctwas"))

  precomputed_ctwas_res <- readRDS("LDL_example.ctwas_sumstats_res.RDS")

  capture.output({
    ctwas_res <- ctwas_sumstats(z_snp,
                                weights,
                                region_info,
                                snp_info,
                                LD_info,
                                thin = 0.1,
                                maxSNP = 20000)
  })

  expect_equal(ctwas_res, precomputed_ctwas_res)

})
