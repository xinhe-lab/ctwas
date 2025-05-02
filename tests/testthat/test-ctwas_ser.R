test_that("ctwas_ser_rss works", {
  n = 1000
  p = 100
  # prior_weights = rep(1/p,p)
  prior_weights = rep(0.005,p)
  prior_variance_beta = 0.05
  prior_variance_z = prior_variance_beta * n
  V = matrix(rep(prior_variance_z, p), nrow = 1)
  simulated_data <- simulate_single_group(n, p,
                                          prior_weights,
                                          prior_variance = V,
                                          effect_scale = "z",
                                          compute_R = TRUE,
                                          seed = 1)
  z <- simulated_data$z
  R <- simulated_data$R

  # Fit SER model (summary statistics version) to z-scores
  ser_rss_res <- ctwas_ser_rss(z,
                               prior_variance = V,
                               prior_weights = prior_weights,
                               residual_variance = 1,
                               null_method = "susie")

  # Fit susie_rss with L = 1
  null_weight = max(0, 1 - sum(prior_weights))
  scaled_prior_weights = prior_weights/(1-null_weight)

  susie_rss_res <- susie_rss(z,
                             R,
                             prior_weights = scaled_prior_weights,
                             prior_variance = V,
                             estimate_prior_variance = FALSE,
                             L = 1,
                             max_iter = 1,
                             null_weight = null_weight,
                             warn_converge_fail = FALSE)

  expect_equal(drop(susie_rss_res$alpha)[1:p], drop(ser_rss_res$alpha))
  expect_equal(drop(susie_rss_res$mu2[1:p]), drop(ser_rss_res$mu2))
})

test_that("ctwas_ser_rss works under the null", {
  n = 1000
  p = 100
  # prior_weights = rep(1/p,p)
  prior_weights = rep(0.005,p)
  prior_variance_beta = 0.05
  prior_variance_z = prior_variance_beta * n
  V <- matrix(rep(prior_variance_z, p), nrow = 1)

  simulated_data <- simulate_single_group(n, p,
                                          prior_weights,
                                          prior_variance = V,
                                          effect_scale = "z",
                                          simulate_null_model = TRUE,
                                          compute_R = TRUE,
                                          seed = 1)
  z <- simulated_data$z
  R <- simulated_data$R

  ser_rss_susie_res <- ctwas_ser_rss(z,
                                     prior_variance = V,
                                     prior_weights = prior_weights,
                                     residual_variance = 1,
                                     null_method = "susie")

  null_weight = max(0, 1 - sum(prior_weights))
  scaled_prior_weights = prior_weights/(1-null_weight)

  susie_rss_res <- susie_rss(z,
                             R,
                             prior_weights = scaled_prior_weights,
                             prior_variance = V,
                             estimate_prior_variance = FALSE,
                             L = 1,
                             max_iter = 1,
                             null_weight = null_weight,
                             warn_converge_fail = FALSE)

  expect_equal(drop(susie_rss_res$alpha)[1:p], drop(ser_rss_susie_res$alpha))
  expect_equal(drop(susie_rss_res$mu2)[1:p], drop(ser_rss_susie_res$mu2))
})


# test_that("ctwas_ser works", {
#   n = 1000
#   p = 100
#   # prior_weights = rep(1/p,p)
#   prior_weights = rep(0.005,p)
#   prior_variance_beta = 0.05
#   V_beta <- matrix(rep(prior_variance_beta, p), nrow = 1)
#
#   simulated_data <- simulate_single_group(n, p,
#                                           prior_weights,
#                                           prior_variance = V_beta,
#                                           effect_scale = "beta",
#                                           compute_R = FALSE,
#                                           seed = 1)
#   X <- simulated_data$X
#   y <- simulated_data$y
#
#   ser_res <- ctwas_ser(X, y,
#                        scaled_prior_variance = V_beta,
#                        prior_weights = prior_weights,
#                        residual_variance = 1,
#                        null_method = "susie")
#
#   # Fit susie with L = 1
#   null_weight = max(0, 1 - sum(prior_weights))
#   scaled_prior_weights = prior_weights/(1-null_weight)
#
#   capture.output({
#     susie_res <- susie(X, y,
#                        scaled_prior_variance = V_beta,
#                        L = 1,
#                        max_iter = 1,
#                        residual_variance = 1,
#                        prior_weights = scaled_prior_weights,
#                        estimate_residual_variance = FALSE,
#                        estimate_prior_variance = FALSE,
#                        null_weight = null_weight,
#                        standardize = TRUE,
#                        intercept = TRUE,
#                        warn_converge_fail = FALSE)
#   })
#
#   expect_equal(drop(susie_res$alpha)[1:p], drop(ser_res$alpha))
#   expect_equal(drop(susie_res$mu2[1:p]), drop(ser_res$mu2))
# })
