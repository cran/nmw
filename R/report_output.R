#' Generate NONMEM Output Report (SB-Output.PDF)
#'
#' Includes the NONMEM PRINT.OUT file content in a PDF report.
#'
#' @param run_dir character, path to the NONMEM run directory (default: current directory)
#' @export
nmw_report_output <- function(run_dir = getwd()) {
  owd <- setwd(run_dir)
  on.exit(setwd(owd))

  PrepPDF("SB-Output.PDF")
  Print.Out <- readLines("PRINT.OUT")
  PrinMTxt(Print.Out, Cex = 0.6)
  ClosePDF()
  message("SB-Output.PDF generated.")
}
