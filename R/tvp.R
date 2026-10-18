# =====================================================================
# L-23: Forward Filtering Backward Sampling and the TVP-AR(1) model
#   measurement:  pi_t  = alpha + rho_t * pi_{t-1} + e_t,   e_t ~ N(0, r)
#   state:        rho_t = rho_{t-1} + eta_t,                eta_t ~ N(0, q)
# =====================================================================

#' Forward Filtering, Backward Sampling for the TVP-AR(1) persistence
#'
#' Runs the Kalman filter forward (month 1 to T) and draws one path of
#' \eqn{\rho_0, \dots, \rho_T} backward (Carter and Kohn, 1994;
#' Fruhwirth-Schnatter, 1994). Used in session L-23.
#'
#' @param ys Inflation minus alpha, \eqn{\pi_t - \alpha}.
#' @param z Last month's inflation, \eqn{\pi_{t-1}}.
#' @param r Measurement-shock variance.
#' @param q State-innovation variance.
#' @param m0,C0 Prior mean and variance of the starting value \eqn{\rho_0}.
#' @return A list: `draw` = one sampled path \eqn{(\rho_0, \dots, \rho_T)};
#'   `filt` = filtered means \eqn{m_t} (real-time estimates). Position
#'   `t + 1` holds month `t`; position 1 holds \eqn{\rho_0}.
#' @examples
#' pi_m <- waifem_data("Country A Monthly")$inflation_mom
#' y <- pi_m[-1]; x <- pi_m[-length(pi_m)]
#' out <- ffbs(y - 0.59, x, r = 0.115, q = 0.0017)
#' head(out$filt)
#' @export
ffbs <- function(ys, z, r, q, m0 = 0.5, C0 = 1) {
  T <- length(ys)
  m <- numeric(T + 1)       # filtered means
  C <- numeric(T + 1)       # filtered variances
  m[1] <- m0; C[1] <- C0    # start: prior for rho_0

  # ---- FORWARD: Kalman filter, month 1 to month T ----
  for (t in 1:T) {
    R <- C[t] + q                         # predicted variance of rho_t
    Q <- z[t]^2 * R + r                   # variance of the forecast of y_t
    K <- R * z[t] / Q                     # Kalman gain
    surprise <- ys[t] - z[t] * m[t]       # forecast error this month
    m[t + 1] <- m[t] + K * surprise       # update the mean
    C[t + 1] <- R - K * z[t] * R          # update the variance
  }

  # ---- BACKWARD: sample a path, month T back to month 0 ----
  rho <- numeric(T + 1)
  rho[T + 1] <- rnorm(1, m[T + 1], sqrt(C[T + 1]))   # last month
  for (t in T:1) {
    B <- C[t] / (C[t] + q)                    # weight on next month's draw
    h <- m[t] + B * (rho[t + 1] - m[t])       # corrected mean
    H <- C[t] * (1 - B)                       # corrected variance
    rho[t] <- rnorm(1, h, sqrt(H))            # draw this month
  }
  list(draw = rho, filt = m)
}

#' The three simple Gibbs steps of the TVP-AR(1) sampler
#'
#' Each is one draw from a known distribution. Priors:
#' \eqn{\alpha \sim N(0, 10)}, \eqn{r \sim IG(2, 0.1)}, \eqn{q \sim IG(3, b_q)}.
#' An inverse-gamma draw is `1 / rgamma(1, a, b)`. Used in L-23.
#'
#' @param y Inflation this month.
#' @param x Inflation last month.
#' @param rho A sampled path \eqn{(\rho_0, \dots, \rho_T)} from [ffbs()].
#' @param r Measurement-shock variance.
#' @param alpha The constant.
#' @param b_q Prior scale of \eqn{q}.
#' @return One draw of \eqn{\alpha}, \eqn{r} or \eqn{q}.
#' @name gibbs_steps
NULL

#' @rdname gibbs_steps
#' @export
draw_alpha <- function(y, x, rho, r) {           # Normal (regression)
  e1 <- y - rho[-1] * x                          # y minus the rho part
  v  <- 1 / (length(y) / r + 1 / 10)             # posterior variance
  rnorm(1, v * sum(e1) / r, sqrt(v))             # posterior draw
}

#' @rdname gibbs_steps
#' @export
draw_r <- function(y, x, rho, alpha) {           # inverse-gamma
  e <- y - alpha - rho[-1] * x                   # measurement errors
  1 / rgamma(1, 2 + length(y) / 2, 0.1 + sum(e^2) / 2)
}

#' @rdname gibbs_steps
#' @export
draw_q <- function(rho, b_q = 0.004) {           # inverse-gamma
  moves <- diff(rho)                             # rho_t - rho_{t-1}
  1 / rgamma(1, 3 + length(moves) / 2, b_q + sum(moves^2) / 2)
}

#' Gibbs sampler for the TVP-AR(1) model
#'
#' Repeats: FFBS for the path of \eqn{\rho_t}, then draws of \eqn{\alpha},
#' \eqn{r} and \eqn{q}. Used in L-23 (and L-24).
#'
#' @param y Inflation this month.
#' @param x Inflation last month.
#' @param n_iter Total draws.
#' @param n_burn Early draws thrown away (burn-in).
#' @param b_q Prior scale of \eqn{q}.
#' @return A list: `rho` = matrix of kept paths (one row per draw, one column
#'   per month); `par` = matrix of kept `alpha`, `r`, `q`.
#' @examples
#' \donttest{
#' pi_m <- waifem_data("Country A Monthly")$inflation_mom
#' y <- pi_m[-1]; x <- pi_m[-length(pi_m)]
#' set.seed(1)
#' g <- gibbs_tvp(y, x, n_iter = 2000, n_burn = 1000)
#' apply(g$par, 2, median)
#' }
#' @export
gibbs_tvp <- function(y, x, n_iter = 15000, n_burn = 5000, b_q = 0.004) {
  alpha <- 0; r <- var(y); q <- 0.002          # starting values
  keep <- n_iter - n_burn
  rho_store <- matrix(NA, keep, length(y))
  par_store <- matrix(NA, keep, 3, dimnames = list(NULL, c("alpha", "r", "q")))
  for (it in 1:n_iter) {
    rho   <- ffbs(y - alpha, x, r, q)$draw      # 1. states (FFBS)
    alpha <- draw_alpha(y, x, rho, r)           # 2. alpha
    r     <- draw_r(y, x, rho, alpha)           # 3. r
    q     <- draw_q(rho, b_q)                   # 4. q
    if (it > n_burn) {                          # 5. keep after burn-in
      rho_store[it - n_burn, ] <- rho[-1]
      par_store[it - n_burn, ] <- c(alpha, r, q)
    }
  }
  list(rho = rho_store, par = par_store)
}

#' Median and 90\% band of a path, month by month
#'
#' @param draws Matrix of draws: one row per draw, one column per month.
#' @return A data frame with columns `lower` (5\%), `median` and `upper` (95\%).
#' @examples
#' summarise_path(matrix(rnorm(3000), 1000, 3))
#' @export
summarise_path <- function(draws) {
  data.frame(lower  = apply(draws, 2, quantile, 0.05),
             median = apply(draws, 2, median),
             upper  = apply(draws, 2, quantile, 0.95))
}
