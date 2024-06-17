# test_that("est_param works", {
#
#   ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_noLD_res.RDS", package = "ctwas"))
#   region_data <- ctwas_res$region_data
#   estimated_param <- ctwas_res$param
#
#   capture.output({
#     param <- suppressWarnings({
#       est_param(region_data,
#                 group_prior_var_structure = "shared_type",
#                 niter_prefit = 3,
#                 niter = 30)
#     })
#   })
#
#   expect_equal(param, estimated_param)
# })
