#' waifem: Bayesian tools for the WAIFEM Module III training
#'
#' Helper functions used in the WAIFEM Module III slides and exercises.
#' Every function has the same name as on the slides, so code from the
#' slides runs unchanged after `library(waifem)`.
#'
#' Sessions and functions:
#' \itemize{
#'   \item L-23 (FFBS, TVP-AR(1)): [ffbs()], [draw_alpha()], [draw_r()],
#'     [draw_q()], [gibbs_tvp()], [summarise_path()]
#'   \item Day 8 L-21 (ARCH/GARCH): [arch_lm()], [garch_vol()], [bayes_vol_band()]
#'   \item Day 8 L-22 (SV, TVP-AR(1)-SV): [ffbs_tv()], [sv_priors],
#'     [gibbs_tvp_sv()]
#'   \item W-05 (fiscal): [bayes_reg()], [post_table()], [predict_revenue()],
#'     [scenario_X()], [scenario_balance()]
#'   \item L-24 (model evaluation): [fit_four_models()], [two_chains()],
#'     [ppc_stats()], [ppc_pvalues()], [std_resid()], [pit()], [loglik()],
#'     [log_bf_g()], [lfo_ar1_vs_sv()], [conv_table()]
#'   \item Data: [waifem_data()], [waifem_file()]; setup: [waifem_check()]
#' }
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom stochvol svsample
#' @importFrom posterior as_draws_array
#' @importFrom loo loo
#' @importFrom stats acf dnorm embed lm median pchisq pnorm predict quantile rgamma rnorm var
NULL
