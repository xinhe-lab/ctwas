test_that("screen_regions works", {

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
    thin <- min(sapply(region_data, "[[", "thin"))
    if (thin < 1){
      region_data <- expand_region_data(region_data,
                                            snp_map,
                                            z_snp,
                                            maxSNP = 20000,
                                            ncore = 2)
    } else {
      region_data <- region_data
    }

    screen_res <- screen_regions(region_data,
                                 group_prior = group_prior,
                                 group_prior_var = group_prior_var,
                                 min_nonSNP_PIP = 0.5,
                                 null_method = "ctwas",
                                 ncore = 2)
  })

  expect_equal(screen_res$screened_region_data, expected_screen_res$screened_region_data)
  expect_equal(screen_res$screen_summary$n_gids, expected_screen_res$screen_summary$n_gids)

})

