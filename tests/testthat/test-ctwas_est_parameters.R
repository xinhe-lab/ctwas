test_that("est_param works", {

  ctwas_res <- readRDS("LDL_example.ctwas_sumstats_noLD_res.RDS")
  region_data <- ctwas_res$region_data
  expected_param <- ctwas_res$param

  capture.output({
    param <- est_param(region_data,
                       niter_prefit = 3,
                       niter = 30)
  })

  expect_equal(param, expected_param)
})
