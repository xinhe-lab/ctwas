test_that("assemble_region_data works", {

  # Load region info
  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  # Load SNP info
  snp_info <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_info.RDS", package = "ctwas"))
  # Load weights
  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))
  # Load z_snp
  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))
  # Load gene z-scores
  z_gene <- readRDS(system.file("extdata/sample_data", "LDL_example.z_gene.RDS", package = "ctwas"))

  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_noLD_res.RDS", package = "ctwas"))

  preassembled_region_data <- ctwas_res$region_data
  preassembled_boundary_genes <- ctwas_res$boundary_genes

  capture.output({
    res <- assemble_region_data(region_info = region_info,
                                z_snp = z_snp,
                                z_gene = z_gene,
                                weights = weights,
                                snp_info = snp_info)
  })

  region_data <- res$region_data
  boundary_genes <- res$boundary_genes

  expect_equal(region_data, preassembled_region_data)
  expect_equal(boundary_genes, preassembled_boundary_genes)

})
