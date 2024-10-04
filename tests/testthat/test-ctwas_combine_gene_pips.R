test_that("anno_susie_alpha_res works", {

  snp_map <- readRDS(system.file("extdata/sample_data", "LDL_example.snp_map.RDS", package = "ctwas"))
  ctwas_res <- readRDS(system.file("extdata/sample_data", "LDL_example.ctwas_sumstats_res.RDS", package = "ctwas"))
  finemap_res <- ctwas_res$finemap_res
  susie_alpha_res <- ctwas_res$susie_alpha_res
  gene_annot <- readRDS(system.file("extdata/sample_data", "LDL_example.gene_annot.RDS", package = "ctwas"))
  colnames(gene_annot)[colnames(gene_annot) == "gene_id"] <- "molecular_id"

  expected_annotated_susie_alpha_res <- readRDS("LDL_example.annotated_susie_alpha_res.RDS")

  capture.output({
    susie_alpha_res$molecular_id <- get_molecular_ids(susie_alpha_res)
    annotated_susie_alpha_res <- anno_susie_alpha_res(susie_alpha_res,
                                                      mapping_table = gene_annot,
                                                      map_by = "molecular_id",
                                                      drop_unmapped = TRUE)
  })

  expect_equal(annotated_susie_alpha_res$id, expected_annotated_susie_alpha_res$id)
  expect_equal(annotated_susie_alpha_res$gene_name, expected_annotated_susie_alpha_res$gene_name)

})

test_that("combine_gene_pips works", {

  annotated_susie_alpha_res <- readRDS("LDL_example.annotated_susie_alpha_res.RDS")
  expected_combined_pip_res <- readRDS("LDL_example.combined_pip_res.RDS")

  capture.output({
    combined_pip_by_context <- combine_gene_pips(annotated_susie_alpha_res,
                                                 group_by = "gene_name",
                                                 by = "context",
                                                 method = "combine_cs",
                                                 filter_cs = FALSE,
                                                 include_cs_id = FALSE,
                                                 include_set_id = FALSE)

    combined_pip_by_type <- combine_gene_pips(annotated_susie_alpha_res,
                                              group_by = "gene_name",
                                              by = "type",
                                              method = "combine_cs",
                                              filter_cs = TRUE,
                                              include_cs_id = FALSE,
                                              include_set_id = FALSE)

    combined_pip_by_group <- combine_gene_pips(annotated_susie_alpha_res,
                                               group_by = "gene_name",
                                               by = "group",
                                               method = "combine_cs",
                                               filter_cs = TRUE,
                                               include_cs_id = FALSE,
                                               include_set_id = FALSE)

    # saveRDS(list("by_context" = combined_pip_by_context,
    #              "by_type" = combined_pip_by_type,
    #              "by_group" = combined_pip_by_group), "LDL_example.combined_pip_res.RDS")

  })

  expect_equal(combined_pip_by_context, expected_combined_pip_res$by_context)
  expect_equal(combined_pip_by_type, expected_combined_pip_res$by_type)
  expect_equal(combined_pip_by_group, expected_combined_pip_res$by_group)

})
