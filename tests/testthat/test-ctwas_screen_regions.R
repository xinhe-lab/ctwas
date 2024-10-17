test_that("screen_regions_noLD works", {

  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))

  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_noLD_res.RDS", package = "ctwas"))
  region_data <- ctwas_res$region_data
  param <- ctwas_res$param
  group_prior <- param$group_prior
  group_prior_var <- param$group_prior_var
  expected_screen_res <- ctwas_res$screen_res
  rm(ctwas_res)

  capture.output({
    screen_res <- screen_regions_noLD(region_data,
                                      group_prior = group_prior,
                                      group_prior_var = group_prior_var,
                                      min_nonSNP_PIP = 0.5,
                                      ncore = 2)
    screened_region_data <- screen_res$screened_region_data

    # expand selected regions with all SNPs
    screened_region_data <- expand_region_data(screened_region_data,
                                               snp_map,
                                               z_snp,
                                               maxSNP = 20000)
    screen_res$screened_region_data <- screened_region_data
  })

  expect_equal(screen_res$screened_region_data, expected_screen_res$screened_region_data)
  expect_equal(screen_res$screen_summary, expected_screen_res$screen_summary)

})

test_that("screen_regions (filter by L) works", {

  LD_map <- readRDS(system.file("extdata/sample_data", "LDL_example.LD_map.RDS", package = "ctwas"))
  skip_if_no_LD_file(LD_map$LD_file)

  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))

  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))

  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_res.RDS", package = "ctwas"))
  region_data <- ctwas_res$region_data
  param <- ctwas_res$param
  group_prior <- param$group_prior
  group_prior_var <- param$group_prior_var
  expected_screen_res <- ctwas_res$screen_res
  rm(ctwas_res)

  capture.output({
    screen_res <- screen_regions(region_data,
                                 LD_map = LD_map,
                                 weights = weights,
                                 group_prior = group_prior,
                                 group_prior_var = group_prior_var,
                                 filter_L = TRUE,
                                 filter_nonSNP_PIP = FALSE,
                                 ncore = 2)
    screened_region_data <- screen_res$screened_region_data

    # expand selected regions with all SNPs
    screened_region_data <- expand_region_data(screened_region_data,
                                               snp_map,
                                               z_snp,
                                               maxSNP = 20000)
    screen_res$screened_region_data <- screened_region_data
  })

  expect_equal(screen_res$screened_region_data, expected_screen_res$screened_region_data)
  expect_equal(screen_res$screened_region_L, expected_screen_res$screened_region_L)
  expect_equal(screen_res$screen_summary, expected_screen_res$screen_summary)

})



test_that("screen_regions (filter by non-SNP PIPs) works", {

  LD_map <- readRDS(system.file("extdata/sample_data", "LDL_example.LD_map.RDS", package = "ctwas"))
  skip_if_no_LD_file(LD_map$LD_file)

  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))

  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))

  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_nonSNP_PIP_res.RDS", package = "ctwas"))
  region_data <- ctwas_res$region_data
  param <- ctwas_res$param
  group_prior <- param$group_prior
  group_prior_var <- param$group_prior_var
  expected_screen_res <- ctwas_res$screen_res
  rm(ctwas_res)

  capture.output({
    screen_res <- screen_regions(region_data,
                                 LD_map = LD_map,
                                 weights = weights,
                                 group_prior = group_prior,
                                 group_prior_var = group_prior_var,
                                 filter_L = FALSE,
                                 filter_nonSNP_PIP = TRUE,
                                 min_nonSNP_PIP = 0.5,
                                 ncore = 2)
    screened_region_data <- screen_res$screened_region_data

    # expand selected regions with all SNPs
    screened_region_data <- expand_region_data(screened_region_data,
                                               snp_map,
                                               z_snp,
                                               maxSNP = 20000)
    screen_res$screened_region_data <- screened_region_data
  })

  expect_equal(screen_res$screened_region_data, expected_screen_res$screened_region_data)
  expect_equal(screen_res$screened_region_L, expected_screen_res$screened_region_L)
  expect_equal(screen_res$screen_summary, expected_screen_res$screen_summary)

})
