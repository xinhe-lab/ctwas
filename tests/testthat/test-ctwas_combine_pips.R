
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
