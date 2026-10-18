# waifem

Bayesian time-series and fiscal tools for the **WAIFEM Module III** training
(Lagos, November 2026): *Bayesian macroeconometrics for central banks and
ministries of finance in English-speaking West Africa*.

Every function has **the same name as on the training slides**, so code from
the slides runs unchanged after `library(waifem)`. The simulated teaching data
(`WAIFEM_Data.xlsx`) ship with the package, so no file paths are needed.

## Installation

```r
# install.packages("remotes")
remotes::install_github("williamcantah/waifem")
```

Offline (e.g. from a USB stick), install the downloaded file:

```r
install.packages("waifem_0.1.0.tar.gz", repos = NULL, type = "source")
```

The package is written in plain R, so no compiler (Rtools) is needed.
Then check that every package used in the training is installed:

```r
library(waifem)
waifem_check()                 # lists anything missing
waifem_check(install = TRUE)   # installs missing packages from CRAN
```

## Data

```r
waifem_sheets()                          # all sheets
a_m <- waifem_data("Country A Monthly")  # monthly data, date as Date
a_q <- waifem_data("Country A Quarterly")
```

Slides that use `read_excel("data/WAIFEM_Data.xlsx", sheet = ...)` work
unchanged; `waifem_data(sheet)` is the path-free alternative.

## Functions by session

| Session | Functions |
|:--|:--|
| L-23: FFBS, TVP-AR(1) | `ffbs()`, `draw_alpha()`, `draw_r()`, `draw_q()`, `gibbs_tvp()`, `summarise_path()` |
| Day 8 L-21: ARCH/GARCH | `arch_lm()`, `garch_vol()`, `bayes_vol_band()` |
| Day 8 L-22: SV, TVP-AR(1)-SV | `ffbs_tv()`, `sv_priors`, `gibbs_tvp_sv()` |
| W-05: Bayesian fiscal modelling | `bayes_reg()`, `post_table()`, `predict_revenue()`, `scenario_X()`, `scenario_balance()` |
| L-24: model evaluation | `fit_four_models()`, `two_chains()`, `ppc_stats()`, `ppc_pvalues()`, `std_resid()`, `pit()`, `loglik()`, `log_bf_g()`, `lfo_ar1_vs_sv()`, `conv_table()` |
| Data and setup | `waifem_data()`, `waifem_sheets()`, `waifem_file()`, `waifem_check()` |

The code is deliberately simple and fully commented. Type a function name
without brackets (e.g. `gibbs_tvp`) to read its source, or `?gibbs_tvp` for help.

## Example: TVP-AR(1) by FFBS and Gibbs sampling

```r
library(waifem)
pi_m <- waifem_data("Country A Monthly")$inflation_mom
y <- pi_m[-1]                 # inflation this month
x <- pi_m[-length(pi_m)]      # inflation last month
set.seed(1)
g <- gibbs_tvp(y, x)          # 15,000 draws, 5,000 burn-in
apply(g$par, 2, median)       # alpha, r, q
path <- summarise_path(g$rho) # median and 90% band of rho_t
```

## Data note

All data are **simulated** for teaching: Country A to Country E behave like
West African economies but are not official statistics.

## Author

William G. Cantah, Department of Data Science and Economic Policy,
University of Cape Coast. Licence: MIT.
