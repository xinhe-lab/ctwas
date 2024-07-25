test_that("screen_regions (no-LD version) work", {

  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))
  z_gene <- readRDS(system.file("extdata/sample_data", "LDL_example.z_gene.RDS", package = "ctwas"))

  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_noLD_res.RDS", package = "ctwas"))
  region_data <- ctwas_res$region_data
  param <- ctwas_res$param
  group_prior <- param$group_prior
  group_prior_var <- param$group_prior_var
  expected_screened_region_data <- ctwas_res$screened_region_data

  capture.output({
    screen_regions_res <- screen_regions(region_data,
                                         use_LD = FALSE,
                                         group_prior = group_prior,
                                         group_prior_var = group_prior_var,
                                         min_nonSNP_PIP = 0.5)
    screened_region_data <- screen_regions_res$screened_region_data

    # Expand screened region_data with all SNPs in the regions
    screened_region_data <- expand_region_data(screened_region_data,
                                               snp_map,
                                               z_snp,
                                               z_gene,
                                               trim_by = "z",
                                               maxSNP = 20000)
  })

  expect_equal(screened_region_data, expected_screened_region_data)

})

test_that("screen_regions (LD version, filter by L) work", {

  LD_map <- readRDS(system.file("extdata/sample_data", "LDL_example.LD_map.RDS", package = "ctwas"))
  skip_if_no_LD_file(LD_map$LD_file)

  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))

  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))
  z_gene <- readRDS(system.file("extdata/sample_data", "LDL_example.z_gene.RDS", package = "ctwas"))

  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_res.RDS", package = "ctwas"))
  region_data <- ctwas_res$region_data
  param <- ctwas_res$param
  group_prior <- param$group_prior
  group_prior_var <- param$group_prior_var

  expected_screened_region_data <- ctwas_res$screened_region_data

  capture.output({
    screen_regions_res <- screen_regions(region_data,
                                         use_LD = TRUE,
                                         LD_map = LD_map,
                                         snp_map = snp_map,
                                         weights = weights,
                                         group_prior = group_prior,
                                         group_prior_var = group_prior_var,
                                         filter_L = TRUE,
                                         filter_nonSNP_PIP = FALSE,
                                         ncore = 6)
    screened_region_data <- screen_regions_res$screened_region_data

    # Expand screened region_data with all SNPs in the regions
    screened_region_data <- expand_region_data(screened_region_data,
                                               snp_map,
                                               z_snp,
                                               z_gene,
                                               trim_by = "z",
                                               maxSNP = 20000)
  })

  expect_equal(screened_region_data, expected_screened_region_data)

})
