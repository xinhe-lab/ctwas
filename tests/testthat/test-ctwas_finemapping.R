test_that("finemap_regions works with no-LD", {

  # finemap a single region with no-LD version
  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))
  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_noLD_res.RDS", package = "ctwas"))
  precomputed_finemap_res <- ctwas_res$finemap_res
  screened_region_data <- ctwas_res$screened_region_data
  param <- ctwas_res$param
  group_prior <- param$group_prior
  group_prior_var <- param$group_prior_var
  rm(ctwas_res)

  region_id <- sample(names(screened_region_data),1)
  precomputed_finemap_res <- precomputed_finemap_res[precomputed_finemap_res$region_id == region_id,]
  rownames(precomputed_finemap_res) <- NULL

  capture.output({
    finemap_res <- suppressWarnings({
      finemap_regions(region_data = screened_region_data[region_id],
                      use_LD = FALSE,
                      L = 1,
                      weights = weights,
                      group_prior = group_prior,
                      group_prior_var = group_prior_var)
    })
  })
  expect_equal(finemap_res, precomputed_finemap_res)
})

# test_that("finemap_regions works with LD", {
#
#   # finemap a single region with LD
#   snp_info <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_info.RDS", package = "ctwas"))
#   LD_info <- readRDS(system.file("extdata/sample_data", "LDL_example.LD_info.RDS", package = "ctwas"))
#   weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))
#
#   ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_res.RDS", package = "ctwas"))
#   precomputed_finemap_res <- ctwas_res$finemap_res
#   screened_region_data <- ctwas_res$screened_region_data
#   param <- ctwas_res$param
#   group_prior <- param$group_prior
#   group_prior_var <- param$group_prior_var
#   rm(ctwas_res)
#
#   region_id <- sample(names(screened_region_data),1)
#   precomputed_finemap_res <- precomputed_finemap_res[precomputed_finemap_res$region_id == region_id,]
#   rownames(precomputed_finemap_res) <- NULL
#
#   finemap_res <- suppressWarnings({
#     finemap_regions(region_data = screened_region_data[region_id],
#                     use_LD = TRUE,
#                     LD_info = LD_info,
#                     snp_info = snp_info,
#                     weights = weights,
#                     group_prior = group_prior,
#                     group_prior_var = group_prior_var,
#                     L = 5)
#   })
#   expect_equal(finemap_res, precomputed_finemap_res)
# })
