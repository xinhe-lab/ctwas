test_that("anno_finemap_res works", {

  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_res.RDS", package = "ctwas"))
  finemap_res <- ctwas_res$finemap_res
  gene_annot <- readRDS(system.file("extdata/sample_data", "LDL_example.gene_annot.RDS", package = "ctwas"))
  colnames(gene_annot)[colnames(gene_annot) == "gene_id"] <- "molecular_id"

  expected_annotated_finemap_res <- readRDS("LDL_example.annotated_finemap_res.RDS")

  capture.output({
    finemap_res$molecular_id <- get_molecular_ids(finemap_res)
    annotated_finemap_res <- anno_finemap_res(finemap_res,
                                              snp_map = snp_map,
                                              mapping_table = gene_annot,
                                              add_gene_annot = TRUE,
                                              map_by = "molecular_id",
                                              drop_unmapped = TRUE,
                                              add_position = TRUE,
                                              use_gene_pos = "mid")
  })

  expect_equal(annotated_finemap_res$gene_name, expected_annotated_finemap_res$gene_name)
  expect_equal(annotated_finemap_res$gene_type, expected_annotated_finemap_res$gene_type)
  expect_equal(annotated_finemap_res$pos, expected_annotated_finemap_res$pos)

})

test_that("compute_combined_pips works", {

  finemap_gene_res <- data.frame(molecular_id = c("a1", "a2", "a3", "a4", "b1", "b2"),
                                 gene_name = c(rep("A", 4), rep("B", 2)),
                                 region_id = c(rep("region_1", 4), rep("region_2", 2)),
                                 susie_pip = c(0.5, 0.2, 0.3, 0.1, 0.4, 0.5),
                                 cs_index = c(1, 1, 2, 0, 1, 2))

  capture.output({
    combined_gene_pips <- compute_combined_pips(finemap_gene_res,
                                                group_by = "gene_name",
                                                method = "combine_cs",
                                                filter_cs = TRUE)
  })

  expect_equal(combined_gene_pips$combined_pip, c(0.79, 0.70))

})
