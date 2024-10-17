test_that("est_param works", {

  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_noLD_res.RDS", package = "ctwas"))
  region_data <- ctwas_res$region_data
  expected_param <- ctwas_res$param

  capture.output({
    param <- est_param(region_data,
                       niter_prefit = 3,
                       niter = 30,
                       min_gene = 1,
                       ncore = 2)
  })

  expect_equal(param, expected_param)
})
