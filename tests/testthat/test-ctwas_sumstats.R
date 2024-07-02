# test_that("ctwas_sumstats works", {
#
#   # Load z_snp
#   z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))
#   # Load weights
#   weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))
#   # Load region_info, snp_info, and LD_info
#   region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
#   snp_info <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_info.RDS", package = "ctwas"))
#   LD_info <- readRDS(system.file("extdata/sample_data", "LDL_example.LD_info.RDS", package = "ctwas"))
#
#   precomputed_ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_res.RDS", package = "ctwas"))
#
#   capture.output({
#     suppressWarnings({
#       ctwas_res <- ctwas_sumstats(z_snp,
#                                   weights,
#                                   region_info,
#                                   snp_info,
#                                   LD_info,
#                                   thin = 0.1,
#                                   niter_prefit = 3,
#                                   niter = 30,
#                                   L = 5,
#                                   group_prior_var_structure = "shared_type",
#                                   screen_method = "nonSNP_PIP",
#                                   maxSNP = 20000,
#                                   min_nonSNP_PIP = 0.5)
#     })
#   })
#
#   expect_equal(ctwas_res, precomputed_ctwas_res)
#
# })
