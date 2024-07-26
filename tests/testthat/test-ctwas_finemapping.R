test_that("finemap_regions works with no-LD", {

  # finemap a single region with no-LD version
  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))
  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_noLD_res.RDS", package = "ctwas"))

  precomputed_finemap_res <- ctwas_res$finemap_res
  screened_region_data <- ctwas_res$screened_region_data
  param <- ctwas_res$param
  group_prior <- param$group_prior
  group_prior_var <- param$group_prior_var

  region_id <- sample(names(screened_region_data),1)
  precomputed_finemap_res <- precomputed_finemap_res[precomputed_finemap_res$region_id == region_id,]
  rownames(precomputed_finemap_res) <- NULL

  capture.output({
    finemap_res <- finemap_regions(region_data = screened_region_data[region_id],
                                   use_LD = FALSE,
                                   weights = weights,
                                   group_prior = group_prior,
                                   group_prior_var = group_prior_var)
  })
  expect_equal(finemap_res, precomputed_finemap_res)
})

test_that("finemap_regions works with LD", {

  # finemap a single region with LD
  LD_map <- readRDS(system.file("extdata/sample_data", "LDL_example.LD_map.RDS", package = "ctwas"))
  skip_if_no_LD_file(LD_map$LD_file)

  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))

  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_res.RDS", package = "ctwas"))

  precomputed_finemap_res <- ctwas_res$finemap_res
  screened_region_data <- ctwas_res$screened_region_data
  L <- ctwas_res$L
  param <- ctwas_res$param
  group_prior <- param$group_prior
  group_prior_var <- param$group_prior_var

  region_id <- sample(names(screened_region_data),1)
  precomputed_finemap_res <- precomputed_finemap_res[precomputed_finemap_res$region_id == region_id,]
  rownames(precomputed_finemap_res) <- NULL

  capture.output({
    finemap_res <- finemap_regions(region_data = screened_region_data[region_id],
                                   use_LD = TRUE,
                                   LD_map = LD_map,
                                   snp_map = snp_map,
                                   weights = weights,
                                   group_prior = group_prior,
                                   group_prior_var = group_prior_var,
                                   L = L[region_id])
  })

  expect_equal(finemap_res, precomputed_finemap_res)
})
