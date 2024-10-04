test_that("summarize_param works", {

  gwas_n <- 343621
  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_noLD_res.RDS", package = "ctwas"))
  param <- ctwas_res$param

  expected_ctwas_parameters <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_parameters.RDS", package = "ctwas"))

  ctwas_parameters <- summarize_param(param, gwas_n)

  expect_equal(ctwas_parameters, expected_ctwas_parameters)

})
