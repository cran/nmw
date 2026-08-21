#' Generate Input Data Summary Report (SA-Input.PDF)
#'
#' Generates a PDF report with input data validation including
#' individuals with no dosing/DV, duplicate records, and summary statistics.
#'
#' @param run_dir character, path to the NONMEM run directory (default: current directory)
#' @param file character, output PDF filename
#' @param model character, run/model name; NULL auto-detects via GetCurModelName()
#' @return Invisibly, a list with the output \code{file}, number of \code{pages},
#'   record/duplicate counts, and counts of subjects without dose or DV records.
#' @export
nmw_report_input <- function(run_dir = getwd(), file = "SA-Input.PDF", model = NULL) {
  owd <- setwd(run_dir); on.exit(setwd(owd))

  ## ---- parsing + statistics: reused verbatim from nmw ----------------------
  CtlName <- if (is.null(model)) GetCurModelName() else model

  XML <- readLines(paste0(CtlName, ".xml"))

  nRec <- GetNRecFromXML(XML)

  FDATA <- ReadFDATA()
  VarStat <- NMVarStat(FDATA)
  IDStat <- NMIDStat(FDATA)

  # Individuals with no dosing / no DV (subsets of the Individuals data.frame)
  NoDose <- IDStat[[1]][IDStat[[1]][, "nAMT"] == 0, , drop = FALSE]
  NoDV   <- IDStat[[1]][IDStat[[1]][, "nDV"] == 0, , drop = FALSE]

  # Duplicate detection (verbatim)
  DuplData <- NULL
  i <- 1
  while (i < nRec) {
    cID <- FDATA[i, "ID"]
    cTIME <- FDATA[i, "TIME"]
    cMDV <- FDATA[i, "MDV"]
    tDATA <- FDATA[FDATA[, "ID"] == cID & FDATA[, "TIME"] == cTIME & FDATA[, "MDV"] == cMDV, ]
    if (length(tDATA[, 1]) > 1) {
      DuplData <- rbind(DuplData, tDATA)
    }
    i <- i + length(tDATA[, 1])
  }

  # DV == 0 -> NA (verbatim; must precede summary(FDATA), follow the stats above)
  FDATA[!is.na(FDATA$DV) & FDATA$DV == 0, "DV"] <- NA

  # capture.output dumps are rendered as block_pre; keep NONMEM's width so the
  # summary() column groups wrap exactly as in the original report.
  oldWidth <- getOption("width"); on.exit(options(width = oldWidth), add = TRUE)
  options(width = 90)
  SummaryDump <- capture.output(summary(FDATA))
  IDStatDump  <- capture.output(IDStat)

  ## ---- PDF via simPDF (replaces PrepPDF/PrinMTxt/AddPage/ClosePDF) ----------
  blocks <- list(
    block_para(paste0("Summary A - Input Data Check   [", CtlName, "]"),
               size = 16, font = 2),
    block_para(sprintf("No. of data records in data set : %d", nRec), size = 10),
    block_rule(),

    # Individuals with no dosing
    if (nrow(NoDose) > 0) block_keep(list(
      block_para("Individuals with no dosing", size = 13, font = 2),
      block_table(round(NoDose, 5), size = 9, family = "mono"))),

    # Individuals with no DV
    if (nrow(NoDV) > 0) block_keep(list(
      block_para("Individuals with no DV", size = 13, font = 2),
      block_table(round(NoDV, 5), size = 9, family = "mono"))),

    # Potentially harmful (duplicate) records
    if (!is.null(DuplData)) block_keep(list(
      block_para("Potentially Harmful Records (duplicates)", size = 13, font = 2),
      block_para("Records sharing the same ID / TIME / MDV.", size = 9, font = 3))),
    if (!is.null(DuplData)) block_table(DuplData, size = 8, family = "mono"),
    if (!is.null(DuplData)) block_spacer(6),

    # Summary of the (DV==0 -> NA adjusted) input data
    block_keep(list(
      block_para("Summary of Input Data", size = 13, font = 2),
      block_para("(DV values of 0 are treated as missing.)", size = 9, font = 3))),
    block_pre(SummaryDump, size = 9),
    block_spacer(6),

    # Variable-level statistics
    block_keep(list(
      block_para("Variable Statistics", size = 13, font = 2),
      block_para("nNA/nZero/nPos counts, uniqueness, range, type and functional dependence on ID/TIME/MDV/EVID.",
                 size = 9, font = 3))),
    block_table(round(VarStat, 5), size = 8, family = "mono"),
    block_spacer(6),

    # Per-individual statistics
    block_keep(list(
      block_para("Individual Statistics", size = 13, font = 2),
      block_para(sprintf("ID Sorted : %s     Time Sorted : %s",
                         IDStat[["ID Sorted"]], IDStat[["Time Sorted"]]),
                 size = 10))),
    block_table(round(IDStat[[1]], 5), size = 8, family = "mono"),
    block_spacer(6),
    block_para("Full list dump (record positions and sort flags):", size = 9, font = 3),
    block_pre(IDStatDump, size = 8))

  doc <- sp_new(file, paper = "letter", family = "Courier", size = 10)
  frame_set(doc, top = doc$H - 45, bottom = 45, left = 45, right = doc$W - 45)
  flow_run(doc, Filter(Negate(is.null), blocks),
           footer = block_para("CONFIDENTIAL", size = 8, align = "center"))
  np <- doc$page_no
  sp_close(doc)

  message(paste0(file, " generated."))
  invisible(list(file = file, pages = np,
                 nRec = nRec, nDupl = if (is.null(DuplData)) 0L else nrow(DuplData),
                 nNoDose = nrow(NoDose), nNoDV = nrow(NoDV)))
}
