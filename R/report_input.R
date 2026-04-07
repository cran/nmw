#' Generate Input Data Summary Report (SA-Input.PDF)
#'
#' Generates a PDF report with input data validation including
#' individuals with no dosing/DV, duplicate records, and summary statistics.
#'
#' @param run_dir character, path to the NONMEM run directory (default: current directory)
#' @export
nmw_report_input <- function(run_dir = getwd()) {
  owd <- setwd(run_dir)
  on.exit(setwd(owd))

  defpar <- par(no.readonly = TRUE)
  defpar$new <- NULL
  on.exit(par(defpar), add = TRUE)

  CtlName <- GetCurModelName()

  XML <- readLines(paste0(CtlName, ".xml"))

  TagList <- c(" NO. OF DATA RECS IN DATA SET:", " NO. OF DATA ITEMS IN DATA SET:")
  TagIndex <- rep(NA, length(TagList))
  for (i in 1:length(TagList)) {
    TagIndex[i] <- ParseOut(TagList[i], XML)
  }
  nRec <- as.integer(TagIndex[1])

  FDATA <- ReadFDATA()
  VarStat <- NMVarStat(FDATA)
  IDStat <- NMIDStat(FDATA)

  PrepPDF("SA-Input.PDF")

  if (length(IDStat[[1]][IDStat[[1]][, "nAMT"] == 0, 1]) > 0) {
    PrinMTxt(capture.output(IDStat[[1]][IDStat[[1]][, "nAMT"] == 0, ]),
             Header1 = "Individuals with no dosing")
  }
  if (length(IDStat[[1]][IDStat[[1]][, "nDV"] == 0, 1]) > 0) {
    PrinMTxt(capture.output(IDStat[[1]][IDStat[[1]][, "nDV"] == 0, ]),
             Header1 = "Individuals with no DV")
  }

  # Duplicate detection
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

  options(width = 128)
  if (!is.null(DuplData)) {
    PrinMTxt(capture.output(DuplData), Header1 = "Potentially Harmful Records", Cex = 0.6)
  }
  FDATA[!is.na(FDATA$DV) & FDATA$DV == 0, "DV"] <- NA
  options(width = 90)
  PrinMTxt(capture.output(summary(FDATA)))
  PrinMTxt(capture.output(VarStat))
  PrinMTxt(capture.output(IDStat))

  ClosePDF()
  message("SA-Input.PDF generated.")
}
