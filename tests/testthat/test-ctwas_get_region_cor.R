test_that("get_region_cor correctly computes correlation matrices", {
  region_id <- "16_71020125_72901251"

  LD_map <- readRDS(system.file("extdata/sample_data", "LDL_example.LD_map.RDS", package = "ctwas"))
  skip_if_no_LD_file(LD_map$LD_file)

  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))
  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_res.RDS", package = "ctwas"))

  screened_region_data <- ctwas_res$screen_res$screened_region_data

  # test computing the correlation matrices
  regiondata <- extract_region_data(screened_region_data, region_id)
  gids <- regiondata$gid
  sids <- regiondata$sid

  cor_res <- get_region_cor(region_id,
                            sids = sids,
                            gids = gids,
                            LD_map = LD_map,
                            weights = weights)

  R_gene <- readRDS(system.file("extdata/sample_data/cor_matrix", paste0("region.", region_id, ".R_gene.RDS"), package = "ctwas"))
  R_snp_gene <- readRDS(system.file("extdata/sample_data/cor_matrix", paste0("region.", region_id, ".R_snp_gene.RDS"), package = "ctwas"))
  R_snp <- readRDS(system.file("extdata/sample_data/cor_matrix", paste0("region.", region_id, ".R_snp.RDS"), package = "ctwas"))
  expected_cor_res <- list("R_snp" = R_snp,
                           "R_snp_gene" = R_snp_gene,
                           "R_gene" = R_gene)

  expect_equal(cor_res, expected_cor_res)

})

test_that("get_region_cor correctly loads correlation matrices", {
  region_id <- "16_71020125_72901251"
  R_gene <- readRDS(system.file("extdata/sample_data/cor_matrix", paste0("region.", region_id, ".R_gene.RDS"), package = "ctwas"))
  R_snp_gene <- readRDS(system.file("extdata/sample_data/cor_matrix", paste0("region.", region_id, ".R_snp_gene.RDS"), package = "ctwas"))
  R_snp <- readRDS(system.file("extdata/sample_data/cor_matrix", paste0("region.", region_id, ".R_snp.RDS"), package = "ctwas"))
  expected_cor_res <- list("R_snp" = R_snp,
                           "R_snp_gene" = R_snp_gene,
                           "R_gene" = R_gene)

  # test loading the correlation matrices
  cor_dir <- system.file("extdata/sample_data", "cor_matrix", package = "ctwas")
  cor_res <- get_region_cor(region_id, cor_dir = cor_dir)

  expect_equal(cor_res, expected_cor_res)
})


test_that("load_region_cor correctly loads correlation matrices", {

  region_id <- "16_71020125_72901251"
  R_gene <- readRDS(system.file("extdata/sample_data/cor_matrix", paste0("region.", region_id, ".R_gene.RDS"), package = "ctwas"))
  R_snp_gene <- readRDS(system.file("extdata/sample_data/cor_matrix", paste0("region.", region_id, ".R_snp_gene.RDS"), package = "ctwas"))
  R_snp <- readRDS(system.file("extdata/sample_data/cor_matrix", paste0("region.", region_id, ".R_snp.RDS"), package = "ctwas"))
  expected_cor_res <- list("R_snp" = R_snp,
                           "R_snp_gene" = R_snp_gene,
                           "R_gene" = R_gene)

  # test loading the correlation matrices
  cor_dir <- system.file("extdata/sample_data", "cor_matrix", package = "ctwas")
  cor_res <- load_region_cor(region_id, cor_dir = cor_dir)
  expect_equal(cor_res, expected_cor_res)

  R_sg_file <- file.path(cor_dir, paste0("region.", region_id, ".R_snp_gene.RDS"))
  R_g_file <- file.path(cor_dir, paste0("region.", region_id, ".R_gene.RDS"))
  R_s_file <- file.path(cor_dir, paste0("region.", region_id, ".R_snp.RDS"))
  cor_res <- load_region_cor(R_sg_file = R_sg_file,
                             R_g_file = R_g_file,
                             R_s_file = R_s_file)
  expect_equal(cor_res, expected_cor_res)
})
