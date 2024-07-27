test_that("ctwas_sumstats works", {

  LD_map <- readRDS(system.file("extdata/sample_data", "LDL_example.LD_map.RDS", package = "ctwas"))
  skip_if_no_LD_file(LD_map$LD_file)

  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))

  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))
  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))

  precomputed_ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_res.RDS", package = "ctwas"))

  set.seed(99)

  capture.output({
    ctwas_res <- ctwas_sumstats(z_snp = z_snp,
                                weights = weights,
                                region_info = region_info,
                                snp_map = snp_map,
                                LD_map = LD_map,
                                thin = 0.1,
                                maxSNP = 20000,
                                filter_L = TRUE,
                                filter_nonSNP_PIP = FALSE,
                                ncore = 4)
  })

  expect_equal(ctwas_res, precomputed_ctwas_res)
  # expect_equal(ctwas_res$z_gene, precomputed_ctwas_res$z_gene)
  # expect_equal(ctwas_res$param, precomputed_ctwas_res$param)
  # expect_equal(ctwas_res$finemap_res, precomputed_ctwas_res$finemap_res)
  # expect_equal(ctwas_res$boundary_genes, precomputed_ctwas_res$boundary_genes)

})

test_that("ctwas_sumstats with nonSNP_PIP filtering works", {

  LD_map <- readRDS(system.file("extdata/sample_data", "LDL_example.LD_map.RDS", package = "ctwas"))
  skip_if_no_LD_file(LD_map$LD_file)

  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))

  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))
  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))

  precomputed_ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_nonSNP_PIP_res.RDS", package = "ctwas"))

  set.seed(99)

  capture.output({
    ctwas_res <- ctwas_sumstats(z_snp = z_snp,
                                weights = weights,
                                region_info = region_info,
                                snp_map = snp_map,
                                LD_map = LD_map,
                                thin = 0.1,
                                maxSNP = 20000,
                                filter_L = FALSE,
                                filter_nonSNP_PIP = TRUE,
                                min_nonSNP_PIP = 0.5,
                                ncore = 4)
  })

  expect_equal(ctwas_res, precomputed_ctwas_res)
  # expect_equal(ctwas_res$z_gene, precomputed_ctwas_res$z_gene)
  # expect_equal(ctwas_res$param, precomputed_ctwas_res$param)
  # expect_equal(ctwas_res$finemap_res, precomputed_ctwas_res$finemap_res)
  # expect_equal(ctwas_res$boundary_genes, precomputed_ctwas_res$boundary_genes)

})
