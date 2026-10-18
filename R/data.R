#' Path to the WAIFEM teaching data file
#'
#' The simulated data set `WAIFEM_Data.xlsx` ships with the package.
#' @return Full path to the Excel file.
#' @examples
#' waifem_file()
#' @export
waifem_file <- function()
  system.file("extdata", "WAIFEM_Data.xlsx", package = "waifem", mustWork = TRUE)

#' Read one sheet of the WAIFEM teaching data
#'
#' Reads a sheet of the simulated data shipped with the package, so no file
#' path is needed. Sheets: `"Country A Quarterly"` ... `"Country E Quarterly"`,
#' `"Country A Monthly"` ... `"Country E Monthly"`, the matching `"Q Truth"`
#' and `"M Truth"` sheets, `"Codebook"` and `"README"`. The `date` column is
#' converted to class `Date`.
#'
#' @param sheet Sheet name, e.g. `"Country A Monthly"`.
#' @return A data frame.
#' @examples
#' head(waifem_data("Country A Quarterly"))
#' waifem_sheets()
#' @export
waifem_data <- function(sheet) {
  d <- as.data.frame(readxl::read_excel(waifem_file(), sheet = sheet))
  if ("date" %in% names(d)) d$date <- as.Date(d$date)
  d
}

#' @rdname waifem_data
#' @export
waifem_sheets <- function() readxl::excel_sheets(waifem_file())
