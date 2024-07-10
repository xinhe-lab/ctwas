test_that("get_region_cor correctly loads correlation matrices", {

  region_id <- "16_71020125_72901251"
  R_gene <- readRDS(system.file("extdata/sample_data/cor_matrix", paste0("region.", region_id, ".R_gene.RDS"), package = "ctwas"))
  R_snp_gene <- readRDS(system.file("extdata/sample_data/cor_matrix", paste0("region.", region_id, ".R_snp_gene.RDS"), package = "ctwas"))
  R_snp <- readRDS(system.file("extdata/sample_data/cor_matrix", paste0("region.", region_id, ".R_snp.RDS"), package = "ctwas"))
  precomputed_cor_res <- list("R_snp" = R_snp,
                              "R_snp_gene" = R_snp_gene,
                              "R_gene" = R_gene)

  # test loading the correlation matrices
  cor_dir <- system.file("extdata/sample_data", "cor_matrix", package = "ctwas")
  cor_res <- get_region_cor(region_id, cor_dir = cor_dir)

  expect_equal(cor_res, precomputed_cor_res)
})

test_that("get_region_cor correctly computes correlation matrices", {
  region_id <- "16_71020125_72901251"

  LD_info <- readRDS(system.file("extdata/sample_data", "LDL_example.LD_info.RDS", package = "ctwas"))
  skip_if_no_LD_matrix(LD_info$LD_matrix)

  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_info <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_info.RDS", package = "ctwas"))
  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))
  ctwas_res <- readRDS("LDL_example.ctwas_sumstats_res.RDS")
  screened_region_data <- ctwas_res$screened_region_data
  # test computing the correlation matrices
  cor_res <- get_region_cor(region_id,
                            region_data = screened_region_data,
                            LD_info = LD_info,
                            snp_info = snp_info,
                            weights = weights)

  R_gene <- readRDS(system.file("extdata/sample_data/cor_matrix", paste0("region.", region_id, ".R_gene.RDS"), package = "ctwas"))
  R_snp_gene <- readRDS(system.file("extdata/sample_data/cor_matrix", paste0("region.", region_id, ".R_snp_gene.RDS"), package = "ctwas"))
  R_snp <- readRDS(system.file("extdata/sample_data/cor_matrix", paste0("region.", region_id, ".R_snp.RDS"), package = "ctwas"))
  precomputed_cor_res <- list("R_snp" = R_snp,
                              "R_snp_gene" = R_snp_gene,
                              "R_gene" = R_gene)

  expect_equal(cor_res, precomputed_cor_res)

})
