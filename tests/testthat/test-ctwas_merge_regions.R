test_that("merge_region_data_noLD works", {

  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))
  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))
  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_noLD_res.RDS", package = "ctwas"))
  z_gene <- ctwas_res$z_gene
  boundary_genes <- ctwas_res$boundary_genes
  region_data <- ctwas_res$region_data
  rm(ctwas_res)

  expected_merge_region_res <- readRDS("LDL_example.merge_region_noLD_res.RDS")

  capture.output({
    merge_region_res <- merge_region_data_noLD(boundary_genes,
                                               region_data,
                                               region_info = region_info,
                                               snp_map = snp_map,
                                               z_snp = z_snp,
                                               z_gene = z_gene,
                                               expand = TRUE,
                                               maxSNP = 20000,
                                               ncore = 2)
  })

  # expect_equal(merge_region_res, expected_merge_region_res)
  expect_equal(merge_region_res$merged_region_data, expected_merge_region_res$merged_region_data)
  expect_equal(merge_region_res$merged_region_info, expected_merge_region_res$merged_region_info)
  expect_equal(merge_region_res$merged_snp_map, expected_merge_region_res$merged_snp_map)
  expect_equal(merge_region_res$merged_region_id_map, expected_merge_region_res$merged_region_id_map)

})

test_that("merge_region_data works", {

  LD_map <- readRDS(system.file("extdata/sample_data", "LDL_example.LD_map.RDS", package = "ctwas"))
  skip_if_no_LD_file(LD_map$LD_file)

  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))
  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))
  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_res.RDS", package = "ctwas"))
  z_gene <- ctwas_res$z_gene
  boundary_genes <- ctwas_res$boundary_genes
  region_data <- ctwas_res$region_data
  group_prior <- ctwas_res$param$group_prior
  group_prior_var <- ctwas_res$param$group_prior_var
  rm(ctwas_res)

  expected_merge_region_res <- readRDS("LDL_example.merge_region_res.RDS")

  capture.output({
    merge_region_res <- merge_region_data(boundary_genes,
                                          region_data,
                                          region_info = region_info,
                                          snp_map = snp_map,
                                          LD_map = LD_map,
                                          weights = weights,
                                          z_snp = z_snp,
                                          z_gene = z_gene,
                                          expand = TRUE,
                                          maxSNP = 20000,
                                          ncore = 2)
  })

  # expect_equal(merge_region_res, expected_merge_region_res)
  expect_equal(merge_region_res$merged_region_data, expected_merge_region_res$merged_region_data)
  expect_equal(merge_region_res$merged_region_info, expected_merge_region_res$merged_region_info)
  expect_equal(merge_region_res$merged_snp_map, expected_merge_region_res$merged_snp_map)
  expect_equal(merge_region_res$merged_region_id_map, expected_merge_region_res$merged_region_id_map)

})
