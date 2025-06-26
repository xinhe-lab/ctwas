test_that("postprocess_merge_regions works", {

  LD_map <- readRDS(system.file("extdata/sample_data", "LDL_example.LD_map.RDS", package = "ctwas"))
  skip_if_no_LD_file(LD_map$LD_file)

  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))
  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))
  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  mapping_table <- readRDS(system.file("extdata/sample_data", "mapping_table.RDS", package = "ctwas"))
  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_res.RDS", package = "ctwas"))
  z_gene <- ctwas_res$z_gene
  region_data <- ctwas_res$region_data
  group_prior <- ctwas_res$param$group_prior
  group_prior_var <- ctwas_res$param$group_prior_var
  finemap_res <- ctwas_res$finemap_res
  susie_alpha_res <- ctwas_res$susie_alpha_res

  expected_postprocess_merge_region_res <- readRDS("LDL_example.postprocess_merge_region_res.RDS")

  capture.output({
    postprocess_merge_region_res <- postprocess_merge_regions(region_info,
                                                              region_data,
                                                              z_snp,
                                                              z_gene,
                                                              weights,
                                                              LD_map,
                                                              snp_map,
                                                              finemap_res = finemap_res,
                                                              susie_alpha_res = susie_alpha_res,
                                                              mapping_table = mapping_table,
                                                              L = 5,
                                                              group_prior = group_prior,
                                                              group_prior_var = group_prior_var,
                                                              pip_thresh = 0.5,
                                                              filter_cs = TRUE,
                                                              maxSNP = 20000,
                                                              save_cor = FALSE,
                                                              cor_dir = NULL,
                                                              verbose = FALSE,
                                                              ncore = 1)
  })

  # saveRDS(postprocess_merge_region_res, "LDL_example.postprocess_merge_region_res.RDS")

  expect_equal(postprocess_merge_region_res, expected_postprocess_merge_region_res)

})

test_that("postprocess_merge_regions_noLD works", {

  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))
  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))
  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_res.RDS", package = "ctwas"))
  z_gene <- ctwas_res$z_gene
  region_data <- ctwas_res$region_data
  group_prior <- ctwas_res$param$group_prior
  group_prior_var <- ctwas_res$param$group_prior_var
  finemap_res <- ctwas_res$finemap_res
  susie_alpha_res <- ctwas_res$susie_alpha_res
  mapping_table <- readRDS(system.file("extdata/sample_data", "mapping_table.RDS", package = "ctwas"))

  expected_postprocess_merge_region_noLD_res <- readRDS("LDL_example.postprocess_merge_region_noLD_res.RDS")

  capture.output({
    postprocess_merge_region_noLD_res <- postprocess_merge_regions_noLD(region_info,
                                                                        region_data,
                                                                        z_snp,
                                                                        z_gene,
                                                                        weights,
                                                                        snp_map,
                                                                        finemap_res = finemap_res,
                                                                        susie_alpha_res = susie_alpha_res,
                                                                        mapping_table = mapping_table,
                                                                        group_prior = group_prior,
                                                                        group_prior_var = group_prior_var,
                                                                        pip_thresh = 0.5,
                                                                        filter_cs = TRUE,
                                                                        maxSNP = 20000,
                                                                        verbose = FALSE,
                                                                        ncore = 1)
  })

  # saveRDS(postprocess_merge_region_noLD_res, "LDL_example.postprocess_merge_region_noLD_res.RDS")

  expect_equal(postprocess_merge_region_noLD_res, expected_postprocess_merge_region_noLD_res)

})
