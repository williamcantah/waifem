#' Check that everything needed for the training is installed
#'
#' Lists the packages used in the WAIFEM Module III slides and exercises and
#' reports which are missing, with the command to install them.
#'
#' @param install If `TRUE`, install missing CRAN packages.
#' @return Invisibly, a logical vector (installed or not) per package.
#' @examples
#' waifem_check()
#' @export
waifem_check <- function(install = FALSE) {
  pkgs <- c("readxl", "vars", "BVAR", "tseries", "dlm", "shrinkTVP", "coda",
            "strucchange", "fGarch", "bayesGARCH", "stochvol", "lmtest",
            "posterior", "loo")
  ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  for (p in pkgs) message(sprintf("  %-12s %s", p, if (ok[p]) "OK" else "MISSING"))
  if (all(ok)) {
    message("All packages for the WAIFEM training are installed.")
  } else {
    miss <- pkgs[!ok]
    cmd <- sprintf('install.packages(c(%s))', paste0('"', miss, '"', collapse = ", "))
    if (install) utils::install.packages(miss) else message("Install with:\n  ", cmd)
  }
  invisible(ok)
}
