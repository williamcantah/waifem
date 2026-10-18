test_that("data loads", {
  d <- waifem_data("Country A Monthly")
  expect_equal(nrow(d), 240)
  expect_s3_class(d$date, "Date")
  expect_true("Country A Quarterly" %in% waifem_sheets())
})

test_that("ffbs returns a path of the right length", {
  set.seed(1)
  pi_m <- waifem_data("Country A Monthly")$inflation_mom
  y <- pi_m[-1]; x <- pi_m[-240]
  out <- ffbs(y - 0.59, x, r = 0.115, q = 0.0017)
  expect_length(out$draw, length(y) + 1)
  expect_length(ffbs_tv(y - 0.59, x, r = rep(0.115, length(y)), q = 0.0017), length(y) + 1)
})

test_that("short samplers run and give sensible numbers", {
  set.seed(1)
  pi_m <- waifem_data("Country A Monthly")$inflation_mom
  y <- pi_m[-1]; x <- pi_m[-240]
  g <- gibbs_tvp(y, x, n_iter = 600, n_burn = 300)
  expect_equal(dim(g$rho), c(300, length(y)))
  s <- gibbs_tvp_sv(y, x, n_iter = 300, n_burn = 100)
  expect_equal(colnames(s$par), c("alpha", "q", "mu", "phi", "sigma"))
})

test_that("bayes_reg recovers a simple regression", {
  set.seed(2)
  X <- cbind(const = 1, x = rnorm(200)); y <- X %*% c(1, 2) + rnorm(200, sd = 0.5)
  f <- bayes_reg(y, X, b0 = c(0, 0), sd0 = c(100, 100), n_iter = 1500, n_burn = 500)
  expect_equal(unname(apply(f$B, 2, median)), c(1, 2), tolerance = 0.1)
  expect_equal(colnames(post_table(f$B)), c("q05", "median", "q95", "P_pos"))
})

test_that("small helpers work", {
  expect_named(arch_lm(rnorm(100)), c("LM", "p_value"))
  expect_length(garch_vol(c(0.5, 0.1, 0.8), rnorm(50)), 50)
  expect_true(log_bf_g(0.9, 50, 2, 50) > 0)
  expect_equal(nrow(summarise_path(matrix(rnorm(300), 100, 3))), 3)
  d <- waifem_data("Country A Quarterly")
  expect_equal(dim(scenario_X(d[nrow(d), ], g_ngdp = 18)), c(4, 5))
})
