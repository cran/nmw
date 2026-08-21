#' Generate Individual PK Parameter Report (SC-IndiPK.PDF)
#'
#' Generates a PDF report with individual PK parameter distributions
#' including summary statistics, histograms, and QQ plots.
#'
#' @param run_dir character, path to the NONMEM run directory (default: current directory)
#' @param file character, output PDF filename
#' @param model accepted for a uniform report-function interface; unused because
#'   this report always reads \file{patab}
#' @return Invisibly, a list with the output \code{file}, number of \code{pages},
#'   parameter count and names, and subject count.
#' @export
nmw_report_indipk <- function(run_dir = getwd(), file = "SC-IndiPK.PDF", model = NULL) {
  owd <- setwd(run_dir); on.exit(setwd(owd))

  defpar <- par(no.readonly = TRUE)
  defpar$new <- NULL
  on.exit(par(defpar), add = TRUE)

  ## ---- patab reading + parameter stats: reused verbatim from nmw -----------
  IndiPK <- unique(read.table("patab", skip = 1, header = TRUE))

  PKParas <- setdiff(colnames(IndiPK), "ID")
  nPKPara <- length(PKParas)

  if ("ETA1" %in% PKParas) IndiPK <- IndiPK[IndiPK$ETA1 != 0, ]

  ## ---- PDF via simPDF ------------------------------------------------------
  doc <- sp_new(file, paper = "letter", family = "Courier", size = 10)
  frame_set(doc, top = doc$H - 45, bottom = 45, left = 45, right = doc$W - 45)

  ## Per-parameter distribution plots: 3 parameters per page, each parameter is
  ## one row of 3 panels (summary / histogram+density / QQ) drawn by
  ## PlotDistribution().  Replaces the old par(mfrow)/AddPage device juggling.
  RowPerPage <- 3
  if (nPKPara > 0) {
    nGroup <- ceiling(nPKPara / RowPerPage)
    for (g in seq_len(nGroup)) {
      idx <- ((g - 1) * RowPerPage + 1):min(g * RowPerPage, nPKPara)
      sp_figure_page(doc, {
        for (i in idx) {
          var.data <- IndiPK[, PKParas[i]]
          PlotDistribution(var.data, PKParas[i], show_cv = TRUE)
        }
        ## pad any short final page so exactly prod(mfrow) panels are drawn
        npad <- (RowPerPage - length(idx)) * 3
        if (npad > 0) for (k in seq_len(npad)) plot.new()
        mtext(outer = TRUE, side = 3, "PK Parameter distribution")
      }, mfrow = c(RowPerPage, 3), oma = c(1, 1, 3, 1))
    }
  }

  ## Individual parameter table (replaces PrinMTxt(capture.output(IndiPK))).
  IndiPKtab <- IndiPK
  rownames(IndiPKtab) <- NULL

  flow_run(doc, list(
    block_para("Individual PK Parameters", size = 13, font = 2),
    block_rule(),
    block_table(IndiPKtab, size = 9, family = "mono")),
    footer = block_para("CONFIDENTIAL", size = 8, align = "center"))

  np <- doc$page_no
  sp_close(doc)
  message(paste(file, "generated."))
  invisible(list(file = file, pages = np, nPKPara = nPKPara,
                 PKParas = PKParas, nSubject = nrow(IndiPK)))
}
