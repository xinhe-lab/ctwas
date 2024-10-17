test_that("get_weight_LD_from_predictdb correctly computes weight LD", {

  load("ENSG00000103152.11.R_wgt.Rd")
  R_wgt_res <- get_weight_LD_from_predictdb(g.cov_table, g.weight_table, convert_cov_to_cor = TRUE)
  expect_equal(R_wgt_res, R_wgt)

  load("ENSG00000183549.10.R_wgt.Rd")
  R_wgt_res <- get_weight_LD_from_predictdb(g.cov_table, g.weight_table, convert_cov_to_cor = TRUE)
  expect_equal(R_wgt_res, R_wgt)
})

