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
                                              filter_cs = FALSE,
                                              include_cs_id = FALSE,
                                              include_set_id = FALSE)

    combined_pip_by_group <- combine_gene_pips(annotated_susie_alpha_res,
                                               group_by = "gene_name",
                                               by = "group",
                                               method = "combine_cs",
                                               filter_cs = FALSE,
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
