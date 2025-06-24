test_that("postprocess_LD_mismatch works", {

  LD_map <- readRDS(system.file("extdata/sample_data", "LDL_example.LD_map.RDS", package = "ctwas"))
  skip_if_no_LD_file(LD_map$LD_file)

  gwas_n <- 343621
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))
  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))
  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_res.RDS", package = "ctwas"))
  z_gene <- ctwas_res$z_gene
  region_data <- ctwas_res$region_data
  group_prior <- ctwas_res$param$group_prior
  group_prior_var <- ctwas_res$param$group_prior_var
  finemap_res <- ctwas_res$finemap_res
  susie_alpha_res <- ctwas_res$susie_alpha_res

  expected_LD_diagnosis_res <- readRDS("LDL_example.LD_diagnosis_res.RDS")

  capture.output({
    set.seed(99)
    postprocess_LD_mismatch_res <- postprocess_LD_mismatch(region_data,
                                                           z_snp,
                                                           weights,
                                                           LD_map,
                                                           snp_map,
                                                           gwas_n,
                                                           finemap_res,
                                                           susie_alpha_res,
                                                           group_prior = group_prior,
                                                           group_prior_var = group_prior_var,
                                                           min_nonSNP_PIP = 0.8,
                                                           filter_cs = TRUE,
                                                           p_diff_thresh = 5e-6,
                                                           pip_thresh = 0.5,
                                                           plot = FALSE,
                                                           maxSNP = 20000,
                                                           ncore = 1)
  })

  problematic_snps <- postprocess_LD_mismatch_res$problematic_snps
  flipped_snps <- postprocess_LD_mismatch_res$flipped_snps
  condz_stats <- postprocess_LD_mismatch_res$condz_stats

  # saveRDS(postprocess_LD_mismatch_res, "LDL_example.postprocess_LD_mismatch_res.RDS")

  expect_equal(postprocess_LD_mismatch_res$problematic_snps, expected_LD_diagnosis_res$problematic_snps)
  expect_equal(postprocess_LD_mismatch_res$condz_stats$p_diff, expected_LD_diagnosis_res$condz_stats$p_diff)
  expect_equal(postprocess_LD_mismatch_res$flipped_snps, expected_LD_diagnosis_res$flipped_snps)

})

