test_that("est_param works", {

  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_noLD_res.RDS", package = "ctwas"))
  region_data <- ctwas_res$region_data
  expected_param <- ctwas_res$param
  rm(ctwas_res)

  capture.output({
    suppressWarnings({
      param <- est_param(region_data,
                         niter_prefit = 3,
                         niter = 50,
                         group_prior_var_structure = "shared_all",
                         null_method = "ctwas",
                         ncore = 2)
    })
  })

  expect_equal(param$group_prior, expected_param$group_prior)
  expect_equal(param$group_prior_var, expected_param$group_prior_var)
  expect_equal(param$group_size, expected_param$group_size)
  expect_equal(param$loglik_iters, expected_param$loglik_iters)

})
