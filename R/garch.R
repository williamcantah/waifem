# =====================================================================
# Day 8 L-21: ARCH / GARCH tools
# =====================================================================

#' Engle's (1982) ARCH-LM test
#'
#' Regresses \eqn{e_t^2} on \eqn{e_{t-1}^2, \dots, e_{t-q}^2};
#' \eqn{LM = T \times R^2 \sim \chi^2_q}. H0: no ARCH (constant variance).
#' A small p-value indicates volatility clustering. Used in Day 8 L-21,
#' Day 8 L-22 and L-24.
#'
#' @param e Residuals (surprises) from the mean equation.
#' @param q Number of lags of squared residuals.
#' @return Named vector: `LM` statistic and `p_value`.
#' @examples
#' arch_lm(rnorm(200))
#' @export
arch_lm <- function(e, q = 6) {
  E    <- embed(e^2, q + 1)            # column 1: e_t^2; others: its lags
  fit  <- lm(E[, 1] ~ E[, -1])         # auxiliary regression
  LM   <- nrow(E) * summary(fit)$r.squared
  c(LM = LM, p_value = pchisq(LM, df = q, lower.tail = FALSE))
}

#' Conditional volatility path of a GARCH(1,1)
#'
#' \eqn{h_t = \omega + \alpha e_{t-1}^2 + \beta h_{t-1}}, started at the
#' sample variance.
#'
#' @param par `c(omega, alpha, beta)`.
#' @param e Residuals (surprises).
#' @return The volatility path \eqn{\sqrt{h_t}}.
#' @examples
#' garch_vol(c(0.7, 0.18, 0.76), rnorm(100, sd = 3))
#' @export
garch_vol <- function(par, e) {
  h <- numeric(length(e))
  h[1] <- var(e)                                   # start at sample variance
  for (t in 2:length(e))
    h[t] <- par[1] + par[2] * e[t - 1]^2 + par[3] * h[t - 1]
  sqrt(h)
}

#' Bayesian GARCH volatility path with bands
#'
#' Computes the volatility path for many posterior draws of a Bayesian
#' GARCH(1,1) (e.g. from `bayesGARCH`) and summarises it month by month.
#'
#' @param draws Matrix with columns `alpha0`, `alpha1`, `beta`.
#' @param e Residuals.
#' @param n Number of (spread-out) draws to use.
#' @return Data frame with `lower` (5\%), `median`, `upper` (95\%) volatility.
#' @examples
#' d <- cbind(alpha0 = runif(50, 0.5, 1), alpha1 = 0.2, beta = 0.7)
#' head(bayes_vol_band(d, rnorm(100, sd = 3), n = 50))
#' @export
bayes_vol_band <- function(draws, e, n = 500) {
  rows <- round(seq(1, nrow(draws), length.out = n))   # spread-out draws
  V <- sapply(rows, function(i)
         garch_vol(draws[i, c("alpha0", "alpha1", "beta")], e))
  data.frame(lower  = apply(V, 1, quantile, 0.05),
             median = apply(V, 1, median),
             upper  = apply(V, 1, quantile, 0.95))
}
