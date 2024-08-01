test_that("screen_regions (no-LD version) work", {

  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))

  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_noLD_res.RDS", package = "ctwas"))
  region_data <- ctwas_res$region_data
  param <- ctwas_res$param
  group_prior <- param$group_prior
  group_prior_var <- param$group_prior_var
  expected_screen_regions_res <- ctwas_res$screen_regions_res

  capture.output({
    screen_regions_res <- screen_regions(region_data,
                                         use_LD = FALSE,
                                         snp_map = snp_map,
                                         z_snp = z_snp,
                                         group_prior = group_prior,
                                         group_prior_var = group_prior_var,
                                         min_nonSNP_PIP = 0.5,
                                         expand = TRUE,
                                         maxSNP = 20000)
  })

  expect_equal(screen_regions_res, expected_screen_regions_res)

})

test_that("screen_regions (LD version, filter by L) work", {

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
  expected_screen_regions_res <- ctwas_res$screen_regions_res

  capture.output({
    screen_regions_res <- screen_regions(region_data,
                                         use_LD = TRUE,
                                         LD_map = LD_map,
                                         snp_map = snp_map,
                                         weights = weights,
                                         z_snp = z_snp,
                                         group_prior = group_prior,
                                         group_prior_var = group_prior_var,
                                         filter_L = TRUE,
                                         filter_nonSNP_PIP = FALSE,
                                         expand = TRUE,
                                         maxSNP = 20000,
                                         ncore = 2)
  })

  expect_equal(screen_regions_res, expected_screen_regions_res)

})
