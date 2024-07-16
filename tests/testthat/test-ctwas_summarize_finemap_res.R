test_that("anno_finemap_res works", {

  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_noLD_res.RDS", package = "ctwas"))
  finemap_res <- ctwas_res$finemap_res
  gene_annot <- readRDS(system.file("extdata/sample_data", "LDL_example.gene_annot.RDS", package = "ctwas"))

  expected_annotated_finemap_res <- readRDS("LDL_example.annotated_finemap_res.RDS")

  capture.output({
    annotated_finemap_res <- anno_finemap_res(finemap_res,
                                              snp_map,
                                              gene_annot,
                                              use_gene_pos = "mid")
  })

  expect_equal(annotated_finemap_res, expected_annotated_finemap_res)

})
