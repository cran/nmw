#' Generate OFV Diagnostic Report (S1-OFV.PDF)
#'
#' Generates a PDF report with objective function value diagnostics including
#' total OFV, AICc, BIC, individual OFV per DV summary, and gradient analysis.
#'
#' @param run_dir character, path to the NONMEM run directory (default: current directory)
#' @param file character, output PDF filename
#' @param model character, run/model name; NULL auto-detects via GetCurModelName()
#' @return Invisibly, a list with the output \code{file}, number of \code{pages},
#'   subject and observation counts, and the objective function value.
#' @export
nmw_report_ofv <- function(run_dir = getwd(), file = "S1-OFV.PDF", model = NULL) {
  owd <- setwd(run_dir)
  on.exit(setwd(owd))

  ## ---- parsing + statistics: reused verbatim from nmw report_ofv.R ----------
  CtlName <- if (is.null(model)) GetCurModelName() else model

  XML <- readLines(paste0(CtlName, ".xml"))
  EXT <- ReadLastTable(paste0(CtlName, ".ext"))   # chained-$EST safe (last table)
  GRD <- ReadLastTable(paste0(CtlName, ".grd"))
  PHI <- ReadLastTable(paste0(CtlName, ".phi"))

  # Final OFV from the last $EST's final-estimate row (ITERATION==-1000000000), read BEFORE
  # the >=0 filter and RmvFixed. Chained-$EST safe: when the last $EST converges in very few
  # iterations the recorded iteration rows are (near-)identical, so RmvFixed drops the constant
  # OBJ column and the old "OFV at max iteration" lookup returned empty -> blank OFV/AICc/BIC.
  OFV <- EXT[EXT[, "ITERATION"] == -1000000000, "OBJ"]

  EXT <- EXT[EXT[, "ITERATION"] >= 0, ]
  EXT <- RmvFixed(EXT)

  nPara <- dim(GRD)[2] - 1

  IterCnt <- max(EXT[, "ITERATION"])
  Grad <- GRD[GRD[, "ITERATION"] == IterCnt, ]

  TagList <- c(
    " NO. OF DATA RECS IN DATA SET:",
    " NO. OF DATA ITEMS IN DATA SET:",
    " TOT. NO. OF OBS RECS:",
    " TOT. NO. OF INDIVIDUALS:",
    "0LENGTH OF THETA:",
    " NO. OF FUNCTION EVALUATIONS USED:",
    " NO. OF SIG. DIGITS IN FINAL EST.:"
  )

  TagIndex <- rep(NA, length(TagList))
  for (i in 1:length(TagList)) {
    TagIndex[i] <- ParseOut(TagList[i], XML)
  }

  nRec <- as.integer(TagIndex[1])
  nItem <- as.integer(TagIndex[2])
  nObs <- as.integer(TagIndex[3])
  nID <- as.integer(TagIndex[4])
  nTheta <- as.integer(TagIndex[5])
  nFuncEval <- as.integer(TagIndex[6])
  nSigDigit <- as.double(TagIndex[7])
  nDV <- nObs

  FDATA <- ReadFDATA()

  merged <- MergeIDStatOFV(FDATA, PHI)
  IDStat3 <- merged$IDStat3
  IDStat3$FRecD <- ifelse(IDStat3$FRec == IDStat3$FDRec, "T", "F")
  IDStat4 <- IDStat3[, -(11:14)]
  IDStat4$OFVpDV <- IDStat4[, "iOFV"] / IDStat4[, "nDV"]
  IDStat5 <- IDStat4[order(IDStat4[, "OFVpDV"], decreasing = TRUE), ]

  AICc <- OFV + 2 * nPara + 2 * nPara * (nPara + 1) / (nDV - nPara - 1)
  SBIC <- OFV + nPara * log(nDV)

  sIOFVpDV <- summary(IDStat5[, "OFVpDV"])
  sdIOFVpDV <- sd(IDStat5[, "OFVpDV"])

  ## ---- PDF via simPDF (replaces PrepPDF/split.screen/AddPage/PrinTxt) -------

  # Page-1 summary text (aligned key:value lists -> block_pre, monospace)
  probLine <- paste("PROBLEM :", XML[grep(" PROBLEM NO.:", XML) + 1])

  ofvLines <- c(
    paste("Number of Total Records :", nRec),
    paste("Number of DV Records    :", nObs),
    paste("Number of Items(Columns):", nItem),
    paste("Number of Parameters    :", nPara),
    paste("Objective Function Value:", OFV),
    paste("OFV per DV              :", format(OFV / nDV, digits = 5)),
    paste("Corrected AIC Value     :", format(AICc, digits = 10)),
    paste("Schwartz Criterion(BIC) :", format(SBIC, digits = 10)),
    paste("# of Gradients Over |1| :", sum(abs(Grad[-1]) > 1)),
    paste("Number of Sig Digits    :", nSigDigit)
  )

  iofvSumLines <- c(
    paste("Minimum :", sIOFVpDV[1]),
    paste("1st Qu. :", sIOFVpDV[2]),
    paste("Median  :", sIOFVpDV[3]),
    paste("Mean    :", sIOFVpDV[4]),
    paste("3rd Qu. :", sIOFVpDV[5]),
    paste("Maximum :", sIOFVpDV[6]),
    paste("Std Dev :", format(sdIOFVpDV, digits = 4)),
    paste("Coe Var :", format(sdIOFVpDV / sIOFVpDV[4], digits = 4)),
    paste("S-W test:", format(shapiro.test(IDStat5[, "OFVpDV"])$p.value, digits = 4))
  )

  headerInfoLines <- c(
    "ID     : Subject ID",
    "iOFV   : Individual Objective Function Value",
    "nRec   : Number of Records",
    "nDV    : Number of DV(Dependent Variable) Records",
    "nMDV   : Number of Missing DV Records",
    "nAMT   : Number of Dosing(AMT) Records",
    "nEVIDx : Number of Records with EVID == x",
    "FRecD  : Is the first record is a dosing record?",
    "OFVpDV : Objective Function Value per DV",
    "*Table is ordered by decreasing OFVpDV"
  )

  abbrevLines <- c(
    "PRED   : Typical Prediction",
    "IPRE   : Individual Prediction",
    "WRES   : Weighted Residual",
    "CWRE   : Conditional Weighted Residual",
    "IWRE   : Individual Weighted Residual",
    "LL     : Lower Limit of Confidence Interval",
    "UL     : Upper Limit of Confidence Interval",
    "RSE    : Relative Standard Error (SE / Point Estimate)",
    "SHR    : Shrinkage (Observed SD / Estimated SD)",
    "ZERO   : Does the confidence interval include 0?",
    "ONE    : Does the confidence interval include 1?"
  )

  # Per-subject iOFV table (was capture.output(IDStat5) hand-painted via split.screen).
  # Round only for display; computation above is untouched.
  iOFVtab <- IDStat5
  rownames(iOFVtab) <- NULL              # clean sequential rank ordering
  iOFVtab$iOFV   <- round(iOFVtab$iOFV, 4)
  iOFVtab$OFVpDV <- round(iOFVtab$OFVpDV, 4)

  blocks <- list(
    block_para("Summary 1 - Objective Function Values", size = 16, font = 2),
    block_rule(),
    block_para(probLine, size = 10),
    block_spacer(4),
    block_keep(list(
      block_para("Objective Function Values", size = 13, font = 2),
      block_pre(ofvLines, size = 9))),
    block_spacer(6),
    block_keep(list(
      block_para("Summary of Individual OFV per DV", size = 13, font = 2),
      block_pre(iofvSumLines, size = 9))),
    block_spacer(6),

    # --- 3-panel iOFV diagnostic (boxplot-style / histogram / qq) -----------
    block_keep(list(
      block_para("Distribution of Individual OFV per DV", size = 13, font = 2),
      block_plot({
        x <- IDStat5[, "OFVpDV"]
        y <- jitter(rep(0, length(x)))
        plot(x, y, ylim = c(-0.2, 0.2), xaxt = "n", yaxt = "n", xlab = "", ylab = "",
             bty = "n", cex = 0.4, main = "iOFV per DV (box summary)")
        lines(c(sIOFVpDV[2], sIOFVpDV[2]), c(-0.05, +0.05))
        lines(c(sIOFVpDV[3], sIOFVpDV[3]), c(-0.05, +0.05))
        lines(c(sIOFVpDV[4], sIOFVpDV[4]), c(-0.07, +0.07))
        lines(c(sIOFVpDV[5], sIOFVpDV[5]), c(-0.05, +0.05))
        lines(c(sIOFVpDV[2], sIOFVpDV[5]), c(-0.05, -0.05))
        lines(c(sIOFVpDV[2], sIOFVpDV[5]), c(+0.05, +0.05))
        lines(c(sIOFVpDV[4] - 2 * sdIOFVpDV, sIOFVpDV[4] - 2 * sdIOFVpDV), c(-0.03, +0.03))
        lines(c(sIOFVpDV[4] + 2 * sdIOFVpDV, sIOFVpDV[4] + 2 * sdIOFVpDV), c(-0.03, +0.03))
      }, height = 150, mar = c(3, 2, 2, 1)))),
    block_plot({
      var.data <- IDStat5[, "OFVpDV"]
      h.res <- hist(var.data, plot = FALSE)
      d.res <- density(var.data, na.rm = TRUE)
      h.rat <- max(h.res$counts) / max(h.res$density)
      xrange <- range(d.res$x)
      yrange <- c(0, max(h.res$counts, max(d.res$y) * h.rat))
      plot(h.res, xlim = xrange, ylim = yrange, xlab = "IOFV per DV", main = "")
      lines(d.res$x, h.rat * d.res$y)
    }, height = 190, mar = c(4, 4, 2, 1)),
    block_plot({
      qqnorm(IDStat5[, "OFVpDV"], datax = TRUE, main = "", ylab = "IOFV per DV", cex = 0.4)
    }, height = 190, mar = c(4, 4, 2, 1)),
    block_spacer(6),

    block_keep(list(
      block_para("Header information for the next table", size = 13, font = 2),
      block_pre(headerInfoLines, size = 9))),
    block_spacer(6),
    block_keep(list(
      block_para("Abbreviations for Tables", size = 13, font = 2),
      block_pre(abbrevLines, size = 9))),
    block_spacer(6),

    block_keep(list(
      block_para("Table of Individual Objective Function Values and Records Summary",
                 size = 13, font = 2),
      block_para("*ordered by decreasing OFVpDV", size = 8))),
    block_table(iOFVtab, size = 8, family = "mono")
  )

  doc <- sp_new(file, paper = "letter", family = "Courier", size = 10)
  frame_set(doc, top = doc$H - 45, bottom = 45, left = 45, right = doc$W - 45)
  flow_run(doc, Filter(Negate(is.null), blocks),
           footer = block_para("CONFIDENTIAL", size = 8, align = "center"))
  np <- doc$page_no
  sp_close(doc)

  message(paste0(file, " generated."))
  invisible(list(file = file, pages = np, nID = nID, nDV = nDV, OFV = OFV))
}
