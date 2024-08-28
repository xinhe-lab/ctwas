test_that("compute_gene_z works", {

  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))
  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))
  precomputed_z_gene <- readRDS(system.file("extdata/sample_data", "LDL_example.z_gene.RDS", package = "ctwas"))

  capture.output({
    z_gene <- compute_gene_z(z_snp, weights, ncore = 2)
  })

  expect_equal(z_gene, precomputed_z_gene)

})
