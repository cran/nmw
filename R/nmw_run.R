#' Run All NONMEM Post-Processing Reports
#'
#' Generates all diagnostic PDF reports for a NONMEM run in sequence:
#' S1-OFV, S2-Parameters, S3-Predictions, S4-Residuals, S5-EBE,
#' SA-Input, SB-Output, SC-IndiPK.
#'
#' @param run_dir character, path to the NONMEM run directory (default: current directory)
#' @param reports character vector, which reports to generate.
#'   Options: "ofv", "param", "pred", "resid", "ebe", "input", "output", "indipk".
#'   Default is all reports.
#' @export
#' @examples
#' \dontrun{
#' # Generate all reports
#' nmw_run("/path/to/nonmem/run/directory")
#'
#' # Generate only OFV and parameter reports
#' nmw_run("/path/to/run", reports = c("ofv", "param"))
#' }
nmw_run <- function(run_dir = getwd(),
                    reports = c("ofv", "param", "pred", "resid", "ebe",
                                "input", "output", "indipk")) {

  report_funcs <- list(
    ofv    = nmw_report_ofv,
    param  = nmw_report_param,
    pred   = nmw_report_pred,
    resid  = nmw_report_resid,
    ebe    = nmw_report_ebe,
    input  = nmw_report_input,
    output = nmw_report_output,
    indipk = nmw_report_indipk
  )

  report_names <- list(
    ofv    = "S1-OFV",
    param  = "S2-Parameters",
    pred   = "S3-Predictions",
    resid  = "S4-Residuals",
    ebe    = "S5-EBE",
    input  = "SA-Input",
    output = "SB-Output",
    indipk = "SC-IndiPK"
  )

  reports <- match.arg(reports, names(report_funcs), several.ok = TRUE)

  message("Starting NONMEM post-processing reports...")
  message(paste("Run directory:", run_dir))
  message(paste("Reports to generate:", paste(reports, collapse = ", ")))
  message("")

  for (rpt in reports) {
    message(paste0("Generating ", report_names[[rpt]], "..."))
    tryCatch({
      report_funcs[[rpt]](run_dir = run_dir)
    }, error = function(e) {
      warning(paste0("Error generating ", report_names[[rpt]], ": ", e$message))
    })
  }

  message("")
  message("All reports completed.")
}
