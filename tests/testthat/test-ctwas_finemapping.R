test_that("finemap_regions_noLD works", {

  # finemap a single region with no-LD version
  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_noLD_res.RDS", package = "ctwas"))
  expected_finemap_res <- ctwas_res$finemap_res
  expected_susie_alpha_res <- ctwas_res$susie_alpha_res

  screened_region_data <- ctwas_res$screen_res$screened_region_data
  param <- ctwas_res$param
  group_prior <- param$group_prior
  group_prior_var <- param$group_prior_var
  rm(ctwas_res)

  region_id <- sample(names(screened_region_data),1)
  expected_finemap_res <- expected_finemap_res[expected_finemap_res$region_id == region_id,]
  rownames(expected_finemap_res) <- NULL

  expected_susie_alpha_res <- expected_susie_alpha_res[expected_susie_alpha_res$region_id == region_id,]
  rownames(expected_susie_alpha_res) <- NULL

  capture.output({
    res <- finemap_regions_noLD(region_data = screened_region_data[region_id],
                                group_prior = group_prior,
                                group_prior_var = group_prior_var)
    finemap_res <- res$finemap_res
    susie_alpha_res <- res$susie_alpha_res
  })

  expect_equal(finemap_res$susie_pip, expected_finemap_res$susie_pip)
  expect_equal(susie_alpha_res$susie_alpha, expected_susie_alpha_res$susie_alpha)

})

test_that("finemap_regions works", {

  # finemap a single region with LD
  LD_map <- readRDS(system.file("extdata/sample_data", "LDL_example.LD_map.RDS", package = "ctwas"))
  skip_if_no_LD_file(LD_map$LD_file)

  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))

  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_res.RDS", package = "ctwas"))
  expected_finemap_res <- ctwas_res$finemap_res
  expected_susie_alpha_res <- ctwas_res$susie_alpha_res
  screened_region_data <- ctwas_res$screen_res$screened_region_data
  screened_region_L <- ctwas_res$screen_res$screened_region_L
  param <- ctwas_res$param
  group_prior <- param$group_prior
  group_prior_var <- param$group_prior_var
  rm(ctwas_res)

  region_id <- sample(names(screened_region_data),1)
  expected_finemap_res <- expected_finemap_res[expected_finemap_res$region_id %in% region_id,]
  rownames(expected_finemap_res) <- NULL

  expected_susie_alpha_res <- expected_susie_alpha_res[expected_susie_alpha_res$region_id %in% region_id,]
  rownames(expected_susie_alpha_res) <- NULL

  capture.output({
    res <- finemap_regions(region_data = screened_region_data[region_id],
                           LD_map = LD_map,
                           weights = weights,
                           group_prior = group_prior,
                           group_prior_var = group_prior_var,
                           L = screened_region_L[region_id])
    finemap_res <- res$finemap_res
    susie_alpha_res <- res$susie_alpha_res
  })

  expect_equal(finemap_res$susie_pip, expected_finemap_res$susie_pip)
  expect_equal(susie_alpha_res$susie_alpha, expected_susie_alpha_res$susie_alpha)

})

