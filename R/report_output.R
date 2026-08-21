#' Generate NONMEM Output Report (SB-Output.PDF)
#'
#' Includes the NONMEM PRINT.OUT file content in a PDF report.
#'
#' @param run_dir character, path to the NONMEM run directory (default: current directory)
#' @param file character, output PDF filename
#' @param model accepted for a uniform report-function interface; unused because
#'   this report always reads \file{PRINT.OUT}
#' @return Invisibly, a list with the output \code{file}, number of \code{pages},
#'   and number of source \code{lines}.
#' @export
nmw_report_output <- function(run_dir = getwd(), file = "SB-Output.PDF", model = NULL) {
  owd <- setwd(run_dir)
  on.exit(setwd(owd))

  ## ---- parsing: reused verbatim from nmw ----------------------------------
  if (!file.exists("PRINT.OUT"))
    stop("PRINT.OUT not found in ", run_dir)
  Print.Out <- readLines("PRINT.OUT")

  ## ---- PDF via simPDF (replaces PrepPDF/PrinMTxt/ClosePDF) -----------------
  ## PrinMTxt(Print.Out, Cex = 0.6)  ->  block_pre(Print.Out, size = 7)
  content <- if (length(Print.Out) == 0L) {
    block_para("(PRINT.OUT is empty)", size = 10, font = 3)
  } else {
    block_pre(Print.Out, size = 7)
  }

  blocks <- list(
    block_keep(list(
      block_para("NONMEM Output (PRINT.OUT)", size = 16, font = 2),
      block_rule())),
    content
  )

  doc <- sp_new(file, paper = "letter", family = "Courier", size = 10)
  frame_set(doc, top = doc$H - 45, bottom = 45, left = 45, right = doc$W - 45)
  flow_run(doc, Filter(Negate(is.null), blocks),
           footer = block_para("CONFIDENTIAL", size = 8, align = "center"))
  np <- doc$page_no
  sp_close(doc)

  message(file, " generated.")
  invisible(list(file = file, pages = np, lines = length(Print.Out)))
}
