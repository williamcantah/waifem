# =====================================================================
# Day 8 L-22: stochastic volatility and the TVP-AR(1)-SV model
#   measurement: pi_t  = alpha + rho_t * pi_{t-1} + exp(h_t / 2) * eps_t
#   state 1:     rho_t = rho_{t-1} + eta_t,                 eta_t ~ N(0, q)
#   state 2:     h_t   = mu + phi * (h_{t-1} - mu) + sigma * u_t
# =====================================================================

#' FFBS with a time-varying measurement variance
#'
#' The L-23 [ffbs()] for \eqn{\rho_t}, but the measurement variance may
#' differ month by month (`r` is a vector, \eqn{r_t = e^{h_t}}). The only
#' change from L-23 is `r[t]` instead of `r` in the forecast variance.
#'
#' @inheritParams ffbs
#' @param r Vector of measurement variances, one per month.
#' @return One sampled path \eqn{(\rho_0, \dots, \rho_T)}.
#' @examples
#' pi_m <- waifem_data("Country A Monthly")$inflation_mom
#' y <- pi_m[-1]; x <- pi_m[-length(pi_m)]
#' head(ffbs_tv(y - 0.59, x, r = rep(0.115, length(y)), q = 0.0017))
#' @export
ffbs_tv <- function(ys, z, r, q, m0 = 0.5, C0 = 1) {
  T <- length(ys); m <- C <- numeric(T + 1)
  m[1] <- m0; C[1] <- C0                       # prior for rho_0
  for (t in 1:T) {                             # ---- forward filter
    R <- C[t] + q                              # predicted variance
    Q <- z[t]^2 * R + r[t]                     # forecast variance (r_t!)
    K <- R * z[t] / Q                          # Kalman gain
    m[t + 1] <- m[t] + K * (ys[t] - z[t] * m[t])
    C[t + 1] <- R - K * z[t] * R
  }
  rho <- numeric(T + 1)                        # ---- backward sampler
  rho[T + 1] <- rnorm(1, m[T + 1], sqrt(C[T + 1]))
  for (t in T:1) {
    B <- C[t] / (C[t] + q)
    rho[t] <- rnorm(1, m[t] + B * (rho[t + 1] - m[t]), sqrt(C[t] * (1 - B)))
  }
  rho
}

#' Priors for the stochastic-volatility part
#'
#' \eqn{\mu \sim N(0, 10)}, \eqn{(\phi + 1)/2 \sim Beta(20, 1.5)},
#' \eqn{\sigma^2 \sim Gamma(0.5, 0.5)} (Nakajima, 2011; Kastner, 2016),
#' in the format of [stochvol::specify_priors()].
#' @format A `stochvol` prior specification.
#' @export
sv_priors <- stochvol::specify_priors(
  mu     = stochvol::sv_normal(mean = 0, sd = 10),
  phi    = stochvol::sv_beta(shape1 = 20, shape2 = 1.5),
  sigma2 = stochvol::sv_gamma(shape = 0.5, rate = 0.5))

#' Gibbs sampler for the TVP-AR(1)-SV model
#'
#' The L-23 sampler plus one stochastic-volatility step: (1) path of
#' \eqn{\rho_t} by [ffbs_tv()] with \eqn{r_t = e^{h_t}}; (2) \eqn{\alpha} by
#' weighted regression; (3) path of \eqn{h_t} and \eqn{\mu, \phi, \sigma}
#' with [stochvol::svsample_fast_cpp()] (Kim, Shephard and Chib mixture
#' sampler with ASIS); (4) \eqn{q} (inverse-gamma). Priors for
#' \eqn{\alpha, q, \rho_0} as in L-23. Used in Day 8 L-22 and L-24.
#'
#' @inheritParams gibbs_tvp
#' @return A list: `rho` (paths of \eqn{\rho_t}), `h` (paths of \eqn{h_t}),
#'   `par` (matrix of `alpha`, `q`, `mu`, `phi`, `sigma`).
#' @examples
#' \donttest{
#' pi_m <- waifem_data("Country A Monthly")$inflation_mom
#' y <- pi_m[-1]; x <- pi_m[-length(pi_m)]
#' set.seed(1)
#' g <- gibbs_tvp_sv(y, x, n_iter = 1500, n_burn = 500)
#' apply(g$par, 2, median)
#' }
#' @export
gibbs_tvp_sv <- function(y, x, n_iter = 15000, n_burn = 5000, b_q = 0.004) {
  T <- length(y); keep <- n_iter - n_burn
  alpha <- 0; q <- 0.002                               # starting values
  h <- rep(log(var(y)), T)                             # start: constant vol
  sv <- list(mu = log(var(y)), phi = 0.9, sigma = 0.3, nu = Inf,
             rho = 0, beta = NA, latent0 = log(var(y)))
  rho_store <- h_store <- matrix(NA, keep, T)
  par_store <- matrix(NA, keep, 5,
                      dimnames = list(NULL, c("alpha", "q", "mu", "phi", "sigma")))
  for (it in 1:n_iter) {
    r   <- exp(h)                                      # r_t = exp(h_t)
    # 1. persistence path rho_t: FFBS with time-varying r_t
    rho <- ffbs_tv(y - alpha, x, r, q)
    # 2. alpha: weighted regression (calm months get more weight)
    e1  <- y - rho[-1] * x
    v   <- 1 / (sum(1 / r) + 1 / 10)
    alpha <- rnorm(1, v * sum(e1 / r), sqrt(v))
    # 3. SV step: log-volatility path h_t and mu, phi, sigma (stochvol)
    e   <- e1 - alpha                                  # inflation surprises
    s   <- stochvol::svsample_fast_cpp(e, startpara = sv, startlatent = h,
                                       priorspec = sv_priors)
    h   <- as.numeric(s$latent)
    sv$mu <- s$para[1, "mu"]; sv$phi <- s$para[1, "phi"]
    sv$sigma <- s$para[1, "sigma"]; sv$latent0 <- as.numeric(s$latent0)
    # 4. q: how much rho_t moves (inverse-gamma, as in L-23)
    q   <- 1 / rgamma(1, 3 + T / 2, b_q + sum(diff(rho)^2) / 2)
    if (it > n_burn) {                                 # 5. store
      k <- it - n_burn
      rho_store[k, ] <- rho[-1]; h_store[k, ] <- h
      par_store[k, ] <- c(alpha, q, sv$mu, sv$phi, sv$sigma)
    }
  }
  list(rho = rho_store, h = h_store, par = par_store)
}
