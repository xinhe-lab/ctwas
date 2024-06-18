test_that("summarize_param works", {

  gwas_n <- 343621
  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_res.RDS", package = "ctwas"))
  param <- ctwas_res$param
  precomputed_ctwas_parameters <- readRDS("LDL_example.ctwas_parameters_no_plots.RDS")

  ctwas_parameters <- summarize_param(param, gwas_n)
  expect_equal(ctwas_parameters, precomputed_ctwas_parameters)

})
