# =====================================================================
# W-05: Bayesian fiscal regression and revenue / balance scenarios
# =====================================================================

#' Bayesian linear regression (two-block Gibbs sampler)
#'
#' Model \eqn{y = Xb + e}, \eqn{e \sim N(0, s^2)}. Independent priors in
#' natural units: \eqn{b_j \sim N(b0_j, sd0_j^2)}, \eqn{s^2 \sim IG(a0, d0)}.
#' Gibbs steps: (1) \eqn{b | s^2} is Normal, a precision-weighted average of
#' prior and data; (2) \eqn{s^2 | b} is inverse-gamma from the squared
#' residuals. Used in W-05 and L-24.
#'
#' @param y Dependent variable.
#' @param X Regressor matrix (include a column of 1s for the constant).
#' @param b0 Prior means, one per column of `X`.
#' @param sd0 Prior standard deviations, one per column of `X`.
#' @param a0,d0 Inverse-gamma prior for the error variance (vague).
#' @param n_iter Total draws.
#' @param n_burn Burn-in draws discarded.
#' @return A list: `B` = draws of the coefficients; `s2` = draws of the
#'   error variance.
#' @examples
#' X <- cbind(const = 1, x = rnorm(50)); y <- X %*% c(1, 2) + rnorm(50)
#' fit <- bayes_reg(y, X, b0 = c(0, 0), sd0 = c(100, 100), n_iter = 1500, n_burn = 500)
#' post_table(fit$B)
#' @export
bayes_reg <- function(y, X, b0, sd0, a0 = 2, d0 = 0.001,
                      n_iter = 6000, n_burn = 1000) {
  y   <- as.numeric(y)                     # accept a vector or 1-column matrix
  k   <- ncol(X)
  V0i <- diag(1 / sd0^2, k)                # prior precision (1 / variance)
  XtX <- crossprod(X); Xty <- crossprod(X, y)
  s2  <- var(y)                            # starting value
  B   <- matrix(NA, n_iter - n_burn, k, dimnames = list(NULL, colnames(X)))
  S2  <- numeric(n_iter - n_burn)
  for (it in 1:n_iter) {
    # 1. coefficients: posterior precision = prior + data precision
    Vn <- solve(V0i + XtX / s2)
    bn <- Vn %*% (V0i %*% b0 + Xty / s2)   # weighted average
    b  <- as.numeric(bn + t(chol(Vn)) %*% rnorm(k))
    # 2. error variance: inverse-gamma from the residuals
    e  <- y - X %*% b
    s2 <- 1 / rgamma(1, a0 + length(y) / 2, d0 + sum(e^2) / 2)
    if (it > n_burn) { B[it - n_burn, ] <- b; S2[it - n_burn] <- s2 }
  }
  list(B = B, s2 = S2)
}

#' Posterior summary table
#'
#' @param B Matrix of draws (one column per coefficient).
#' @return For each column: 5\%, median, 95\% and P(coefficient > 0).
#' @examples
#' post_table(cbind(a = rnorm(1000, 1), b = rnorm(1000)))
#' @export
post_table <- function(B) {
  round(t(apply(B, 2, function(b)
    c(q05 = quantile(b, 0.05, names = FALSE), median = median(b),
      q95 = quantile(b, 0.95, names = FALSE), P_pos = mean(b > 0)))), 3)
}

#' Predictive draws of revenue for future quarters
#'
#' Each draw uses one set of coefficients and (optionally) a new random
#' shock, so the result contains parameter and shock uncertainty.
#'
#' @param fit Output of [bayes_reg()] for the log-revenue equation.
#' @param Xnew One row per future quarter, same columns as `X`
#'   (e.g. from [scenario_X()]).
#' @param shocks Add shock uncertainty? `FALSE` keeps parameter uncertainty only.
#' @return Matrix of revenue (LCU bn): draws x future quarters.
#' @export
predict_revenue <- function(fit, Xnew, shocks = TRUE) {
  n  <- nrow(fit$B)
  mu <- fit$B %*% t(Xnew)                                  # log revenue
  if (shocks) mu <- mu + matrix(rnorm(length(mu)), n) * sqrt(fit$s2)
  exp(mu)                                                  # back to LCU bn
}

#' Future regressor rows for the tax equation under a scenario
#'
#' Builds four quarters ahead of the columns `const`, `ly` (log nominal GDP),
#' `lm` (log import share), `lc` (log commodity index) and `tax_policy`.
#'
#' @param last Last row of the quarterly data (a data frame row with
#'   `nominal_gdp_bn`, `imports_bn`, `commodity_index`, `tax_policy`).
#' @param g_ngdp Nominal GDP growth, \% per year.
#' @param com_chg Change in the commodity index over the year, \%.
#' @param imp_chg Change in the import share of GDP over the year, \%.
#' @param policy Number of tax-policy measures in force in each quarter.
#' @return A 4 x 5 matrix.
#' @examples
#' d <- waifem_data("Country A Quarterly")
#' scenario_X(d[nrow(d), ], g_ngdp = 18)
#' @export
scenario_X <- function(last, g_ngdp, com_chg = 0, imp_chg = 0, policy = rep(last$tax_policy, 4)) {
  h <- 1:4                                                  # 4 quarters ahead
  cbind(const      = 1,
        ly         = log(last$nominal_gdp_bn) + h * log(1 + g_ngdp / 100) / 4,
        lm         = log(last$imports_bn / last$nominal_gdp_bn) + h * log(1 + imp_chg / 100) / 4,
        lc         = log(last$commodity_index) + h * log(1 + com_chg / 100) / 4,
        tax_policy = policy)
}

#' Primary balance next year under a scenario
#'
#' Combines revenue draws with the expenditure equation to give draws of
#' annual revenue, expenditure and primary balance (\% of GDP). Debt is held
#' at its last value; debt dynamics are modelled in W-06.
#'
#' @param rev_draws Draws x 4 quarters of revenue (from [predict_revenue()]).
#' @param ngdp Four quarters of nominal GDP in the scenario (LCU bn).
#' @param efit [bayes_reg()] fit of the expenditure equation with columns
#'   const, revenue_gdp, output_gap, election, fiscal_rule, debt_l.
#' @param gap Output gap in each of the 4 quarters.
#' @param election,rule Election and fiscal-rule dummies.
#' @param debt Debt/GDP (held fixed).
#' @return Data frame of draws: `revenue`, `expenditure`, `primary_bal`.
#' @export
scenario_balance <- function(rev_draws, ngdp, efit, gap, election = 0,
                             rule = 1, debt) {
  n <- nrow(rev_draws)
  rev_gdp <- 100 * rev_draws / matrix(ngdp, n, 4, byrow = TRUE)  # % of GDP
  exp_gdp <- matrix(NA, n, 4)
  for (h in 1:4) {
    Xe <- cbind(1, rev_gdp[, h], gap[h], election, rule, debt)  # one row per draw
    exp_gdp[, h] <- rowSums(efit$B * Xe) + rnorm(n) * sqrt(efit$s2)
  }
  w <- ngdp / sum(ngdp)                                          # GDP weights
  data.frame(revenue     = rev_gdp %*% w,
             expenditure = exp_gdp %*% w,
             primary_bal = (rev_gdp - exp_gdp) %*% w)
}
