# =====================================================================
# L-24: Bayesian model evaluation
# =====================================================================

#' Fit the week's four inflation models
#'
#' AR1 (constant persistence and variance), AR1_SV (stochastic volatility),
#' TVP (time-varying persistence), TVP_SV (both). Uses [bayes_reg()],
#' [stochvol::svsample()], [gibbs_tvp()] and [gibbs_tvp_sv()].
#'
#' @param y Inflation this month.
#' @param x Inflation last month.
#' @param seed Random seed.
#' @return A named list; each element holds `mu` and `sd`, the predicted mean
#'   and standard deviation of \eqn{\pi_t} (draws x months).
#' @export
fit_four_models <- function(y, x, seed = 1) {
  T <- length(y); X1 <- cbind(1, x); xr <- function(n) matrix(x, n, T, byrow = TRUE)
  set.seed(seed); f1 <- bayes_reg(y, X1, b0 = c(0, 0), sd0 = c(10, 10),
                                  n_iter = 3000, n_burn = 1000)
  set.seed(seed); f2 <- stochvol::svsample(c(x[1], y), designmatrix = "ar1",
                                 priorphi = c(20, 1.5), draws = 4000,
                                 burnin = 2000, thin = 2, quiet = TRUE)
  set.seed(seed); f3 <- gibbs_tvp(y, x, n_iter = 9000, n_burn = 5000)
  set.seed(seed); f4 <- gibbs_tvp_sv(y, x, n_iter = 9000, n_burn = 5000)
  list(
    AR1    = list(mu = f1$B %*% t(X1),
                  sd = matrix(sqrt(f1$s2), nrow(f1$B), T)),
    AR1_SV = list(mu = as.matrix(f2$beta) %*% t(X1),
                  sd = exp(as.matrix(stochvol::latent(f2)) / 2)),
    TVP    = list(mu = f3$par[, "alpha"] + f3$rho * xr(nrow(f3$rho)),
                  sd = matrix(sqrt(f3$par[, "r"]), nrow(f3$rho), T)),
    TVP_SV = list(mu = f4$par[, "alpha"] + f4$rho * xr(nrow(f4$rho)),
                  sd = exp(f4$h / 2)))
}

#' Put two chains into the format of the posterior package
#'
#' @param a,b Matrices of parameter draws (draws x parameters), one per chain.
#' @return A `draws_array` (draws x chains x parameters).
#' @export
two_chains <- function(a, b) {
  arr <- array(NA, c(nrow(a), 2, ncol(a)), dimnames = list(NULL, NULL, colnames(a)))
  arr[, 1, ] <- a; arr[, 2, ] <- b
  posterior::as_draws_array(arr)
}

#' Posterior predictive check statistics and p-values
#'
#' `ppc_stats()` computes three statistics with economic meaning: fat tails
#' (kurtosis), volatility clustering (autocorrelation of squared deviations)
#' and the worst (highest) month. `ppc_pvalues()` simulates a fake series for
#' `n` draws and counts how often fake >= observed. p near 0 or 1 means the
#' model cannot reproduce that feature of the data.
#'
#' @param v A numeric series.
#' @param fit One element of [fit_four_models()] (with `mu` and `sd`).
#' @param y Observed series.
#' @param n Number of simulated series.
#' @return `ppc_stats()`: named vector. `ppc_pvalues()`: list with `p`,
#'   `sims` and `obs`.
#' @export
ppc_stats <- function(v) {
  e <- v - mean(v)
  c(fat_tails  = mean(e^4) / mean(e^2)^2,             # kurtosis
    clustering = acf(e^2, 1, plot = FALSE)$acf[2],     # volatility clustering
    worst_month = max(v))                               # highest inflation month
}

#' @rdname ppc_stats
#' @export
ppc_pvalues <- function(fit, y, n = 1000) {
  rows <- round(seq(1, nrow(fit$mu), length.out = n))
  sims <- t(sapply(rows, function(i)
    ppc_stats(rnorm(length(y), fit$mu[i, ], fit$sd[i, ]))))
  obs  <- ppc_stats(y)
  list(p = colMeans(sims >= matrix(obs, n, 3, byrow = TRUE)), sims = sims, obs = obs)
}

#' Standardised residuals at posterior medians
#'
#' @param fit Element of [fit_four_models()].
#' @param y Observed series.
#' @return \eqn{(y - \hat\mu) / \hat\sigma}.
#' @export
std_resid <- function(fit, y)
  (y - apply(fit$mu, 2, median)) / apply(fit$sd, 2, median)

#' Probability integral transform
#'
#' \eqn{u_t = P(\pi_t^{pred} \le \pi_t | data)}, averaged over posterior
#' draws. For a calibrated model \eqn{u_t \sim U(0, 1)}.
#'
#' @inheritParams std_resid
#' @return Vector of PIT values, one per month.
#' @export
pit <- function(fit, y)
  colMeans(pnorm(matrix(y, nrow(fit$mu), length(y), byrow = TRUE), fit$mu, fit$sd))

#' Pointwise log-likelihood matrix for WAIC and LOO
#'
#' @inheritParams std_resid
#' @param n Number of (spread-out) draws to use.
#' @return Matrix draws x months, ready for [loo::loo()] and [loo::waic()].
#' @export
loglik <- function(fit, y, n = 2000) {
  rows <- round(seq(1, nrow(fit$mu), length.out = n))
  dnorm(matrix(y, n, length(y), byrow = TRUE), fit$mu[rows, ], fit$sd[rows, ], log = TRUE)
}

#' Log Bayes factor under Zellner's g-prior
#'
#' Log Bayes factor of a regression against the intercept-only model
#' (Liang et al., 2008):
#' \eqn{\log BF = \frac{n-1-k}{2}\log(1+g) - \frac{n-1}{2}\log(1 + g(1-R^2))}.
#'
#' @param R2 R-squared of the regression.
#' @param n Number of observations.
#' @param k Number of regressors (excluding the intercept).
#' @param g Prior scale (`g = n` is the "unit information" prior).
#' @return Log Bayes factor.
#' @examples
#' log_bf_g(R2 = c(0.93, 0.96), n = 77, k = c(2, 4), g = 77)
#' @export
log_bf_g <- function(R2, n, k, g)
  ((n - 1 - k) / 2) * log(1 + g) - ((n - 1) / 2) * log(1 + g * (1 - R2))

#' Exact leave-future-out validation: AR1 versus AR1-SV
#'
#' For each of the last `H` months: refit on data up to \eqn{t-1}, predict
#' month \eqn{t}, and record the log predictive density.
#'
#' @param y Inflation this month.
#' @param x Inflation last month.
#' @param H Number of final months to predict.
#' @return Matrix H x 2 of log predictive densities (`AR1`, `AR1_SV`).
#' @export
lfo_ar1_vs_sv <- function(y, x, H = 72) {
  T <- length(y); out <- matrix(NA, H, 2, dimnames = list(NULL, c("AR1", "AR1_SV")))
  for (j in 1:H) {
    t <- T - H + j; yy <- y[1:(t - 1)]; xx <- x[1:(t - 1)]
    set.seed(j)
    f1 <- bayes_reg(yy, cbind(1, xx), b0 = c(0, 0), sd0 = c(10, 10),
                    n_iter = 2500, n_burn = 500)
    out[j, 1] <- log(mean(dnorm(y[t], f1$B %*% c(1, x[t]), sqrt(f1$s2))))
    f2 <- stochvol::svsample(c(x[1], yy), designmatrix = "ar1", priorphi = c(20, 1.5),
                             draws = 2000, burnin = 1000, quiet = TRUE)
    pr <- suppressWarnings(predict(f2, steps = 1, newdata = matrix(c(1, x[t]), 1)))
    h  <- as.numeric(stochvol::predlatent(pr)[[1]])     # predicted log-volatility
    out[j, 2] <- log(mean(dnorm(y[t], as.matrix(f2$beta) %*% c(1, x[t]), exp(h / 2))))
  }
  out
}

#' Convergence table for two chains
#'
#' R-hat, bulk-ESS and tail-ESS (Vehtari et al., 2021). Rules of thumb:
#' R-hat < 1.01 and ESS > 400.
#'
#' @param a,b Matrices of parameter draws, one per chain.
#' @return Data frame with `rhat`, `ess_bulk`, `ess_tail`.
#' @export
conv_table <- function(a, b) {
  s <- posterior::summarise_draws(two_chains(a, b), "rhat", "ess_bulk", "ess_tail")
  data.frame(rhat = round(s$rhat, 3), ess_bulk = round(s$ess_bulk),
             ess_tail = round(s$ess_tail), row.names = s$variable)
}
