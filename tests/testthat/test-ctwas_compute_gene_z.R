test_that("compute_gene_z works", {

  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))
  z_snp <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.z_snp.RDS", package = "ctwas"))
  precomputed_z_gene <- readRDS(system.file("extdata/sample_data", "LDL_example.z_gene.RDS", package = "ctwas"))

  capture.output({
    z_gene <- compute_gene_z(z_snp, weights, ncore = 2)
  })

  expect_equal(z_gene, precomputed_z_gene)

})

test_that("get_gene_info works", {

  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))
  expected_gene_info <- readRDS(system.file("extdata/sample_data", "LDL_example.gene_info.RDS", package = "ctwas"))

  capture.output({
    gene_info <- get_gene_info(weights)
  })

  # saveRDS(gene_info, "inst/extdata/sample_data/LDL_example.gene_info.RDS")

  expect_equal(gene_info, expected_gene_info)

})

test_that("get_boundary_genes works", {

  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))
  z_gene <- readRDS(system.file("extdata/sample_data", "LDL_example.z_gene.RDS", package = "ctwas"))
  expected_boundary_genes <- readRDS(system.file("extdata/sample_data", "LDL_example.boundary_genes.RDS", package = "ctwas"))

  capture.output({
    boundary_genes <- get_boundary_genes(region_info,
                                         weights,
                                         gene_ids = z_gene$id,
                                         ncore = 2)
  })

  # saveRDS(boundary_genes, "inst/extdata/sample_data/LDL_example.boundary_genes.RDS")

  expect_equal(boundary_genes, expected_boundary_genes)

})

test_that("get_boundary_genes with mapping_table works", {

  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))
  z_gene <- readRDS(system.file("extdata/sample_data", "LDL_example.z_gene.RDS", package = "ctwas"))
  mapping_table <- readRDS(system.file("extdata/sample_data", "mapping_table.RDS", package = "ctwas"))
  expected_boundary_genes <- readRDS(system.file("extdata/sample_data", "LDL_example.combined_boundary_genes.RDS", package = "ctwas"))

  capture.output({
    boundary_genes <- get_boundary_genes(region_info,
                                         weights,
                                         gene_ids = z_gene$id,
                                         mapping_table = mapping_table,
                                         show_mapping = TRUE,
                                         ncore = 2)
  })

  # saveRDS(boundary_genes, "inst/extdata/sample_data/LDL_example.combined_boundary_genes.RDS")

  expect_equal(boundary_genes, expected_boundary_genes)

})


test_that("select_boundary_genes works", {

  region_info <- readRDS(system.file("extdata/sample_data", "LDL_example.region_info.RDS", package = "ctwas"))
  weights <- readRDS(system.file("extdata/sample_data", "LDL_example.preprocessed.weights.RDS", package = "ctwas"))
  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_res.RDS", package = "ctwas"))
  z_gene <- ctwas_res$z_gene
  finemap_res <- ctwas_res$finemap_res
  susie_alpha_res <- ctwas_res$susie_alpha_res

  combined_pip_res <- readRDS("LDL_example.combined_pip_res.RDS")
  combined_pip_res <- combined_pip_res$by_context
  mapping_table <- readRDS(system.file("extdata/sample_data", "mapping_table.RDS", package = "ctwas"))

  expected_postprocess_merge_region_res <- readRDS("LDL_example.postprocess_merge_region_res.RDS")

  capture.output({
    res <- select_boundary_genes(region_info,
                                 weights,
                                 gene_ids = z_gene$id,
                                 finemap_res = finemap_res,
                                 susie_alpha_res = susie_alpha_res,
                                 combine_PIPs = TRUE,
                                 mapping_table = mapping_table,
                                 pip_thresh = 0.5,
                                 filter_cs = TRUE,
                                 ncore = 1)
    boundary_genes <- res$boundary_genes
    selected_boundary_genes <- res$selected_boundary_genes
  })

  expect_equal(boundary_genes, expected_postprocess_merge_region_res$boundary_genes)
  expect_equal(selected_boundary_genes, expected_postprocess_merge_region_res$selected_boundary_genes)

})
