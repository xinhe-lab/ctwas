test_that("diagnose_LD_mismatch_susie works", {

  LD_map <- readRDS(system.file("extdata/sample_data", "LDL_example.LD_map.RDS", package = "ctwas"))
  skip_if_no_LD_file(LD_map$LD_file)

  gwas_n <- 343621
  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))
  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_res.RDS", package = "ctwas"))
  finemap_res <- ctwas_res$finemap_res

  expected_problematic_snps <- readRDS("LDL_example.problematic_snps.RDS")
  expected_LD_diagnosis_res <- readRDS("LDL_example.LD_diagnosis_res.RDS")

  capture.output({
    set.seed(99)
    nonSNP_PIPs <- compute_region_nonSNP_PIPs(finemap_res, filter_cs = TRUE)
    selected_region_ids <- names(nonSNP_PIPs[nonSNP_PIPs > 0.8])

    LD_diagnosis_res <- diagnose_LD_mismatch_susie(selected_region_ids,
                                                   z_snp,
                                                   LD_map,
                                                   gwas_n,
                                                   ncore = 2)

    problematic_snps <- LD_diagnosis_res$problematic_snps
    flipped_snps <- LD_diagnosis_res$flipped_snps
    condz_stats <- LD_diagnosis_res$condz_stats
  })

  expect_equal(problematic_snps, expected_problematic_snps)
  expect_equal(LD_diagnosis_res$condz_stats$p_diff, expected_LD_diagnosis_res$condz_stats$p_diff)

})

