test_that("est_param works", {

  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_noLD_v0.5_res.RDS", package = "ctwas"))
  region_data <- ctwas_res$region_data
  expected_param <- ctwas_res$param

  capture.output({
    suppressWarnings({
      param <- est_param(region_data,
                         niter_prefit = 3,
                         niter = 30,
                         group_prior_var_structure = "shared_all",
                         null_method = "ctwas",
                         run_enrichment_test = TRUE,
                         include_loglik = FALSE,
                         ncore = 2)
    })
  })

  expect_equal(param, expected_param)

})
