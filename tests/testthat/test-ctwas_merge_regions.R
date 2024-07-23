test_that("merge_region_data works", {

  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))
  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))
  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_noLD_res.RDS", package = "ctwas"))
  z_gene <- ctwas_res$z_gene
  boundary_genes <- ctwas_res$boundary_genes
  region_data <- ctwas_res$region_data

  precomputed_merged_region_data <- readRDS(system.file("tests/testthat", "LDL_example.merged_region_data.RDS", package = "ctwas"))

  capture.output({
    res <- merge_region_data(boundary_genes,
                             region_data,
                             region_info = region_info,
                             snp_map = snp_map,
                             z_snp = z_snp,
                             z_gene = z_gene,
                             use_LD = FALSE,
                             maxSNP = 20000)
    merged_region_data <- res$merged_region_data
    merged_region_info <- res$merged_region_info
    merged_LD_map <- res$merged_LD_map
    merged_snp_map <- res$merged_snp_map
    merged_region_id_map <- res$merged_region_id_map
    merged_region_L <- res$merged_region_L
  })

  expect_equal(merged_region_data, precomputed_merged_region_data)

})
test_that("merge_region_data and finemapping without LD works", {

  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))
  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))
  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_noLD_res.RDS", package = "ctwas"))
  z_gene <- ctwas_res$z_gene
  boundary_genes <- ctwas_res$boundary_genes
  region_data <- ctwas_res$region_data
  group_prior <- ctwas_res$param$group_prior
  group_prior_var <- ctwas_res$param$group_prior_var

  precomputed_merged_region_data <- readRDS(system.file("tests/testthat", "LDL_example.merged_region_data.RDS", package = "ctwas"))
  precomputed_finemap_merged_regions_res <- readRDS(system.file("tests/testthat", "LDL_example.finemap_merged_regions_noLD_res.RDS", package = "ctwas"))

  capture.output({
    res <- merge_region_data(boundary_genes,
                             region_data,
                             region_info = region_info,
                             snp_map = snp_map,
                             z_snp = z_snp,
                             z_gene = z_gene,
                             use_LD = FALSE,
                             maxSNP = 20000)
    merged_region_data <- res$merged_region_data
    merged_region_info <- res$merged_region_info
    merged_snp_map <- res$merged_snp_map
    merged_region_id_map <- res$merged_region_id_map

    finemap_merged_regions_res <- finemap_regions(merged_region_data,
                                                  use_LD = FALSE,
                                                  group_prior = group_prior,
                                                  group_prior_var = group_prior_var)
  })

  expect_equal(merged_region_data, precomputed_merged_region_data)
  expect_equal(finemap_merged_regions_res, precomputed_finemap_merged_regions_res)

})

test_that("merge_region_data and finemapping with LD works", {

  LD_map <- readRDS(system.file("extdata/sample_data", "LDL_example.LD_map.RDS", package = "ctwas"))
  skip_if_no_LD_matrix(LD_map$LD_matrix)

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

  precomputed_merged_region_data <- readRDS(system.file("tests/testthat", "LDL_example.merged_region_data.RDS", package = "ctwas"))
  precomputed_finemap_merged_regions_res <- readRDS(system.file("tests/testthat", "LDL_example.finemap_merged_regions_res.RDS", package = "ctwas"))

  capture.output({
    res <- merge_region_data(boundary_genes,
                             region_data,
                             region_info = region_info,
                             snp_map = snp_map,
                             LD_map = LD_map,
                             z_snp = z_snp,
                             z_gene = z_gene,
                             use_LD = TRUE,
                             estimate_L = TRUE,
                             maxSNP = 20000)
    merged_region_data <- res$merged_region_data
    merged_region_info <- res$merged_region_info
    merged_LD_map <- res$merged_LD_map
    merged_snp_map <- res$merged_snp_map
    merged_region_id_map <- res$merged_region_id_map
    merged_region_L <- res$merged_region_L

    finemap_merged_regions_res <- finemap_regions(merged_region_data,
                                                  use_LD = TRUE,
                                                  LD_map = merged_LD_map,
                                                  snp_map = merged_snp_map,
                                                  weights = weights,
                                                  group_prior = group_prior,
                                                  group_prior_var = group_prior_var,
                                                  L = merged_region_L)
  })

  expect_equal(merged_region_data, precomputed_merged_region_data)
  expect_equal(finemap_merged_regions_res, precomputed_finemap_merged_regions_res)

})
