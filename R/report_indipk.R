#' Generate Individual PK Parameter Report (SC-IndiPK.PDF)
#'
#' Generates a PDF report with individual PK parameter distributions
#' including summary statistics, histograms, and QQ plots.
#'
#' @param run_dir character, path to the NONMEM run directory (default: current directory)
#' @export
nmw_report_indipk <- function(run_dir = getwd()) {
  owd <- setwd(run_dir)
  on.exit(setwd(owd))

  defpar <- par(no.readonly = TRUE)
  defpar$new <- NULL
  on.exit(par(defpar), add = TRUE)

  PrepPDF("SC-IndiPK.PDF")
  IndiPK <- unique(read.table("patab", skip = 1, header = TRUE))

  PKParas <- setdiff(colnames(IndiPK), "ID")
  nPKPara <- length(PKParas)

  if ("ETA1" %in% PKParas) IndiPK <- IndiPK[IndiPK$ETA1 != 0, ]

  RowPerPage <- 3
  for (i in 1:nPKPara) {
    if (i %% RowPerPage == 1) {
      par(oma = c(1, 1, 3, 1), mfrow = c(RowPerPage, 3), lty = 1)
    }
    var.data <- IndiPK[, PKParas[i]]

    PlotDistribution(var.data, PKParas[i], show_cv = TRUE)

    if (i %% RowPerPage == 0 | i == nPKPara) {
      mtext(outer = TRUE, side = 3, "PK Parameter distribution")
    }
  }

  PrinMTxt(capture.output(IndiPK), Cex = 0.8)
  ClosePDF()
  message("SC-IndiPK.PDF generated.")
}
