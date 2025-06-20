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
                                               maxSNP = 20000,
                                               ncore = 2)
  })

  # saveRDS(merge_region_res, "LDL_example.merge_region_noLD_res.RDS")

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
                                          maxSNP = 20000,
                                          ncore = 2)
  })

  # saveRDS(merge_region_res, "LDL_example.merge_region_res.RDS")

  # expect_equal(merge_region_res, expected_merge_region_res)
  expect_equal(merge_region_res$merged_region_data, expected_merge_region_res$merged_region_data)
  expect_equal(merge_region_res$merged_region_info, expected_merge_region_res$merged_region_info)
  expect_equal(merge_region_res$merged_snp_map, expected_merge_region_res$merged_snp_map)
  expect_equal(merge_region_res$merged_region_id_map, expected_merge_region_res$merged_region_id_map)

})

test_that("create_merged_snp_map works", {

  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_noLD_res.RDS", package = "ctwas"))
  boundary_genes <- ctwas_res$boundary_genes

  expected_merge_region_res <- readRDS("LDL_example.merge_region_noLD_res.RDS")

  capture.output({
    res <- create_merged_snp_map(boundary_genes, region_info, snp_map)
    merged_region_info <- res$merged_region_info
    merged_snp_map <- res$merged_snp_map
    merged_region_id_map <- res$merged_region_id_map
  })

  expect_equal(merged_region_info, expected_merge_region_res$merged_region_info)
  expect_equal(merged_snp_map, expected_merge_region_res$merged_snp_map)
  expect_equal(merged_region_id_map, expected_merge_region_res$merged_region_id_map)

})

test_that("create_merged_snp_LD_map works", {

  LD_map <- readRDS(system.file("extdata/sample_data", "LDL_example.LD_map.RDS", package = "ctwas"))
  skip_if_no_LD_file(LD_map$LD_file)

  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_noLD_res.RDS", package = "ctwas"))
  boundary_genes <- ctwas_res$boundary_genes

  expected_merge_region_res <- readRDS("LDL_example.merge_region_res.RDS")

  capture.output({
    res <- create_merged_snp_LD_map(boundary_genes, region_info, snp_map, LD_map)
    merged_region_info <- res$merged_region_info
    merged_LD_map <- res$merged_LD_map
    merged_snp_map <- res$merged_snp_map
    merged_region_id_map <- res$merged_region_id_map
  })

  expect_equal(merged_region_info, expected_merge_region_res$merged_region_info)
  expect_equal(merged_LD_map, expected_merge_region_res$merged_LD_map)
  expect_equal(merged_snp_map, expected_merge_region_res$merged_snp_map)
  expect_equal(merged_region_id_map, expected_merge_region_res$merged_region_id_map)

})


test_that("label_overlapping_regions works", {

  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))
  z_gene <- readRDS(system.file("extdata/sample_data", "LDL_example.z_gene.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  mapping_table <- readRDS(system.file("extdata/sample_data", "mapping_real_data.RDS", package = "ctwas"))

  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_noLD_res.RDS", package = "ctwas"))
  boundary_genes <- ctwas_res$boundary_genes

  expected_labeled_boundary_genes <- readRDS("LDL_example.labeled_boundary_genes.RDS")

  capture.output({
    labeled_boundary_genes <- label_overlapping_regions(boundary_genes)
  })
  # saveRDS(labeled_boundary_genes, "LDL_example.labeled_boundary_genes.RDS")

  expect_equal(labeled_boundary_genes, expected_labeled_boundary_genes)

})
