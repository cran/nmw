#' Generate OFV Diagnostic Report (S1-OFV.PDF)
#'
#' Generates a PDF report with objective function value diagnostics including
#' total OFV, AICc, BIC, individual OFV per DV summary, and gradient analysis.
#'
#' @param run_dir character, path to the NONMEM run directory (default: current directory)
#' @export
nmw_report_ofv <- function(run_dir = getwd()) {
  owd <- setwd(run_dir)
  on.exit(setwd(owd))

  defpar <- par(no.readonly = TRUE)
  defpar$new <- NULL
  on.exit(par(defpar), add = TRUE)

  CtlName <- GetCurModelName()

  XML <- readLines(paste0(CtlName, ".xml"))
  EXT <- read.table(paste0(CtlName, ".ext"), skip = 1, header = TRUE)
  GRD <- read.table(paste0(CtlName, ".grd"), skip = 1, header = TRUE)
  PHI <- read.table(paste0(CtlName, ".phi"), skip = 1, header = TRUE)

  EXT <- EXT[EXT[, "ITERATION"] >= 0, ]
  EXT <- RmvFixed(EXT)

  IDs <- PHI[, "ID"]
  nID <- length(IDs)

  params <- CountEXTParams(EXT)
  nTheta <- params$nTheta
  nEta <- params$nEta
  nEps <- params$nEps
  nPara <- dim(GRD)[2] - 1

  IterCnt <- max(EXT[, "ITERATION"])
  OFV <- EXT[EXT[, "ITERATION"] == IterCnt, "OBJ"]
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
  IDStat6 <- capture.output(IDStat5)

  AICc <- OFV + 2 * nPara + 2 * nPara * (nPara + 1) / (nDV - nPara - 1)
  SBIC <- OFV + nPara * log(nDV)

  sIOFVpDV <- summary(IDStat5[, "OFVpDV"])
  sdIOFVpDV <- sd(IDStat5[, "OFVpDV"])

  # --- PDF Generation ---
  PrepPDF("S1-OFV.PDF")

  split.screen(OFV_SCREEN_LAYOUT)
  options(width = 90)

  AddPage()

  screen(1)
  PrinTxt(1, 1, "Summary 1 - Objective Function Values", Cex = 1.2)
  PrinTxt(3, 1, paste("PROBLEM :", XML[grep(" PROBLEM NO.:", XML) + 1]), Cex = 0.9)
  PrinTxt(5, 3, paste("Number of Total Records :", nRec))
  PrinTxt(6, 3, paste("Number of DV Records    :", nObs))
  PrinTxt(7, 3, paste("Number of Items(Columns):", nItem))
  PrinTxt(8, 3, paste("Number of Parameters    :", nPara))
  PrinTxt(9, 3, paste("Objective Function Value:", OFV))
  PrinTxt(10, 3, paste("OFV per DV              :", format(OFV / nDV, digits = 5)))
  PrinTxt(11, 3, paste("Corrected AIC Value     :", format(AICc, digits = 10)))
  PrinTxt(12, 3, paste("Schwartz Criterion(BIC) :", format(SBIC, digits = 10)))
  PrinTxt(13, 3, paste("# of Gradients Over |1| :", sum(abs(Grad[-1]) > 1)))
  PrinTxt(14, 3, paste("Number of Sig Digits    :", nSigDigit))

  PrinTxt(16, 1, "Summary of Individual OFV per DV", Cex = 0.9)
  PrinTxt(18, 3, paste("Minimum :", sIOFVpDV[1]))
  PrinTxt(19, 3, paste("1st Qu. :", sIOFVpDV[2]))
  PrinTxt(20, 3, paste("Median  :", sIOFVpDV[3]))
  PrinTxt(21, 3, paste("Mean    :", sIOFVpDV[4]))
  PrinTxt(22, 3, paste("3rd Qu. :", sIOFVpDV[5]))
  PrinTxt(23, 3, paste("Maximum :", sIOFVpDV[6]))
  PrinTxt(24, 3, paste("Std Dev :", format(sdIOFVpDV, digits = 4)))
  PrinTxt(25, 3, paste("Coe Var :", format(sdIOFVpDV / sIOFVpDV[4], digits = 4)))
  PrinTxt(26, 3, paste("S-W test:", format(shapiro.test(IDStat5[, "OFVpDV"])$p.value, digits = 4)))

  PrinTxt(28, 1, "Header information for the next table", Cex = 0.9)
  PrinTxt(30, 3, "ID     : Subject ID")
  PrinTxt(31, 3, "iOFV   : Individual Objective Function Value")
  PrinTxt(32, 3, "nRec   : Number of Records")
  PrinTxt(33, 3, "nDV    : Number of DV(Dependent Variable) Records")
  PrinTxt(34, 3, "nMDV   : Number of Missing DV Records")
  PrinTxt(35, 3, "nAMT   : Number of Dosing(AMT) Recods")
  PrinTxt(36, 3, "nEVIDx : Number of Records with EVID == x")
  PrinTxt(37, 3, "FRecD  : Is the first record is a dosing record?")
  PrinTxt(38, 3, "OFVpDV : Objective Function Value per DV")
  PrinTxt(39, 3, "*Table is ordered by decreasing OFVpDV", Cex = 0.7)

  PrinTxt(41, 1, "Abbreviations for Tables", Cex = 0.9)
  PrinTxt(43, 3, "PRED   : Typical Prediction")
  PrinTxt(44, 3, "IPRE   : Individual Prediction")
  PrinTxt(45, 3, "WRES   : Weighted Residual")
  PrinTxt(46, 3, "CWRE   : Conditional Weighted Residual")
  PrinTxt(47, 3, "IWRE   : Individual Weighted Residual")
  PrinTxt(48, 3, "LL     : Lower Limit of Confidence Interval")
  PrinTxt(49, 3, "UL     : Upper Limit of Confidence Interval")
  PrinTxt(50, 3, "RSE    : Relative Standard Error (SE / Point Estimate)")
  PrinTxt(51, 3, "SHR    : Shrinkage (Observed SD / Estimated SD)")
  PrinTxt(52, 3, "ZERO   : Does the confidence interval include 0?")
  PrinTxt(53, 3, "ONE    : Does the confidence interval include 1?")

  par(defpar)

  screen(2)
  x <- IDStat5[, "OFVpDV"]
  y <- jitter(rep(0, length(x)))
  plot(x, y, ylim = c(-0.2, 0.2), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", cex = 0.4)
  lines(c(sIOFVpDV[2], sIOFVpDV[2]), c(-0.05, +0.05))
  lines(c(sIOFVpDV[3], sIOFVpDV[3]), c(-0.05, +0.05))
  lines(c(sIOFVpDV[4], sIOFVpDV[4]), c(-0.07, +0.07))
  lines(c(sIOFVpDV[5], sIOFVpDV[5]), c(-0.05, +0.05))
  lines(c(sIOFVpDV[2], sIOFVpDV[5]), c(-0.05, -0.05))
  lines(c(sIOFVpDV[2], sIOFVpDV[5]), c(+0.05, +0.05))
  lines(c(sIOFVpDV[4] - 2 * sdIOFVpDV, sIOFVpDV[4] - 2 * sdIOFVpDV), c(-0.03, +0.03))
  lines(c(sIOFVpDV[4] + 2 * sdIOFVpDV, sIOFVpDV[4] + 2 * sdIOFVpDV), c(-0.03, +0.03))

  par(defpar)

  screen(3)
  var.data <- IDStat5[, "OFVpDV"]
  h.res <- hist(var.data, plot = FALSE)
  d.res <- density(var.data, na.rm = TRUE)
  h.rat <- max(h.res$counts) / max(h.res$density)
  xrange <- range(d.res$x)
  yrange <- c(0, max(h.res$counts, max(d.res$y) * h.rat))
  plot(h.res, xlim = xrange, ylim = yrange, xlab = "IOFV per DV", main = "")
  lines(d.res$x, h.rat * d.res$y)

  screen(4)
  qqnorm(IDStat5[, "OFVpDV"], datax = TRUE, main = "", ylab = "IOFV per DV", cex = 0.4)

  close.screen(all.screens = TRUE)

  nRow <- 55
  AddPage(Cex = 0.8)
  PrinTxt(1, 1, "Table of Individual Objective Function Values and Records Summary", Cex = 1)
  PrinTxt(2, 1, IDStat6[1])

  for (i in 2:length(IDStat6)) {
    if (i %% (nRow - 1) == 1) {
      AddPage(Cex = 0.8)
      PrinTxt(1, 1, IDStat6[1])
    }
    PrinTxt((i - 1) %% (nRow - 1) + 2, 1, IDStat6[i])
  }

  ClosePDF()
  message("S1-OFV.PDF generated.")
}
