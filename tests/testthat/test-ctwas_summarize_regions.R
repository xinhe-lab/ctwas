test_that("get_regions_minP works", {

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

  region_summary <- summarize_region_signals(region_data)
  selected_region_ids <- region_summary[which(region_summary$min_gene_p < 5e-8 | region_summary$min_snp_p < 5e-8), "region_id"]

})

