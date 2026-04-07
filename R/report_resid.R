#' Generate Residual Diagnostics Report (S4-Residuals.PDF)
#'
#' Generates a PDF report with WRES, CWRES, IWRES diagnostics including
#' scatter plots, histograms, normality tests, binomial count tests,
#' and individual residual curves with run tests.
#'
#' @param run_dir character, path to the NONMEM run directory (default: current directory)
#' @export
nmw_report_resid <- function(run_dir = getwd()) {
  owd <- setwd(run_dir)
  on.exit(setwd(owd))

  defpar <- par(no.readonly = TRUE)
  defpar$new <- NULL
  on.exit(par(defpar), add = TRUE)

  RunNumber <- GetCurModelName()

  PHI <- read.table(paste0(RunNumber, ".phi"), skip = 1, header = TRUE)
  XML <- readLines(paste0(RunNumber, ".xml"))

  TagList <- c(" NO. OF DATA RECS IN DATA SET:", " NO. OF DATA ITEMS IN DATA SET:")
  TagIndex <- rep(NA, length(TagList))
  for (i in 1:length(TagList)) {
    TagIndex[i] <- ParseOut(TagList[i], XML)
  }
  nRec <- as.integer(TagIndex[1])

  sdtab <- read.table("sdtab", header = TRUE, skip = 1, na.strings = "***********")
  FDATA <- ReadFDATA()
  FDATA[is.na(FDATA)] <- 0
  sdtab[, "ID"] <- FDATA[, "ID"]

  merged <- MergeIDStatOFV(FDATA, PHI)
  VarStat <- merged$VarStat
  IDStat3 <- merged$IDStat3
  nDV <- VarStat["MDV", "nZero"]

  IDStat4 <- IDStat3[, -10]
  IDStat4$OFVpDV <- IDStat4[, "iOFV"] / IDStat4[, "nDV"]
  IDStat5 <- IDStat4[order(IDStat4[, "OFVpDV"], decreasing = TRUE), ]

  CWRES <- sdtab[sdtab[, "MDV"] == 0, c("ID", "TIME", "DV", "PRED", "IPRE", "WRES", "CWRES", "IWRE")]
  colnames(CWRES) <- c("ID", "TIME", "DV", "PRED", "IPRE", "WRES", "CWRE", "IWRE")

  # Calculate TaLD
  if ("AMT" %in% colnames(FDATA)) {
    tald_result <- CalcTaLDForReport(FDATA, nRec)
    TaLD <- tald_result$TaLD
    NewID <- tald_result$NewID

    FITDATA3 <- cbind(NewID[sdtab[, "MDV"] == 0], TaLD[sdtab[, "MDV"] == 0],
                      CWRES[, c("PRED", "IPRE", "WRES", "CWRE", "IWRE")])
    colnames(FITDATA3) <- c("ID", "TALD", "PRED", "IPRE", "WRES", "CWRE", "IWRE")
  } else {
    FITDATA3 <- CWRES
    FITDATA3$TALD <- FITDATA3$TIME
  }

  resRun <- vector()
  resRun2 <- vector()

  PrepPDF("S4-Residuals.PDF")

  # Residual scatter plots
  par(defpar)
  par(oma = c(1, 1, 3, 1), mfrow = c(3, 2), lty = 1)
  DxPlotPost(FITDATA3$PRED, FITDATA3$WRES, mat = FITDATA3, xlbl = "PRED",
             ylbl = paste("WRES, SD =", format(sd(FITDATA3$WRES, na.rm = TRUE), digits = 4)), smooth = "T")
  DxPlotPost(FITDATA3$TALD, FITDATA3$WRES, mat = FITDATA3, xlbl = "Time after latest dose", ylbl = "WRES", smooth = "T")
  DxPlotPost(FITDATA3$IPRE, FITDATA3$CWRE, mat = FITDATA3, xlbl = "IPRE",
             ylbl = paste("CWRES, SD =", format(sd(FITDATA3$CWRE, na.rm = TRUE), digits = 4)), smooth = "T")
  DxPlotPost(FITDATA3$TALD, FITDATA3$CWRE, mat = FITDATA3, xlbl = "Time after latest dose", ylbl = "CWRES", smooth = "T")
  DxPlotPost(FITDATA3$IPRE, FITDATA3$IWRE, mat = FITDATA3, xlbl = "IPRE",
             ylbl = paste("IWRES, SD =", format(sd(FITDATA3$IWRE, na.rm = TRUE), digits = 4)), smooth = "T")
  DxPlotPost(FITDATA3$TALD, FITDATA3$IWRE, mat = FITDATA3, xlbl = "Time after latest dose", ylbl = "IWRES", smooth = "T")
  mtext(outer = TRUE, side = 3, paste("Residuals of Model", toupper(RunNumber)))

  # Absolute residuals
  par(oma = c(1, 1, 3, 1), mfrow = c(3, 2), lty = 1)
  DxPlotPost(FITDATA3$PRED, abs(FITDATA3$WRES), mat = FITDATA3, xlbl = "PRED", ylbl = "|WRES|", smooth = "T")
  DxPlotPost(FITDATA3$TALD, abs(FITDATA3$WRES), mat = FITDATA3, xlbl = "Time after latest dose", ylbl = "|WRES|", smooth = "T")
  DxPlotPost(FITDATA3$IPRE, abs(FITDATA3$CWRE), mat = FITDATA3, xlbl = "IPRE", ylbl = "|CWRES|", smooth = "T")
  DxPlotPost(FITDATA3$TALD, abs(FITDATA3$CWRE), mat = FITDATA3, xlbl = "Time after latest dose", ylbl = "|CWRES|", smooth = "T")
  DxPlotPost(FITDATA3$IPRE, abs(FITDATA3$IWRE), mat = FITDATA3, xlbl = "IPRE", ylbl = "|IWRES|", smooth = "T")
  DxPlotPost(FITDATA3$TALD, abs(FITDATA3$IWRE), mat = FITDATA3, xlbl = "Time after latest dose", ylbl = "|IWRES|", smooth = "T")
  mtext(outer = TRUE, side = 3, paste("Absolute Values of Residuals of Model", toupper(RunNumber)))

  # Histograms and QQ plots
  par(oma = c(1, 1, 3, 1), mfrow = c(3, 2), lty = 1)

  for (vname in c("WRES", "CWRE", "IWRE")) {
    var.data <- FITDATA3[, vname]
    h.res <- hist(var.data, plot = FALSE)
    d.res <- density(var.data, na.rm = TRUE)
    h.rat <- max(h.res$counts) / max(h.res$density)
    xrange <- range(d.res$x)
    yrange <- c(0, max(h.res$counts, max(d.res$y) * h.rat))
    plot(h.res, xlim = xrange, ylim = yrange,
         xlab = paste(vname, "SD =", format(sd(var.data), digits = 4)), main = "")
    lines(d.res$x, h.rat * d.res$y)
    if (length(var.data) < 5000) {
      p.value <- format(shapiro.test(var.data)$p.value, digits = 4)
      mtext(paste("Shapiro-Wilk test p-value =", p.value), cex = 0.7, line = 0, adj = 0)
    }
    qqnorm(var.data, main = "", ylab = vname, cex = 0.4)
  }
  mtext(outer = TRUE, side = 3, paste("Distribution of Residuals of Model", toupper(RunNumber)))

  # Binomial test
  BTest.Res <- ResTest(CWRES)
  options(width = 110)
  PrinMTxt(capture.output(BTest.Res), Header1 = "Test of residual counts using binomial distribution", Cex = 0.8)

  # Count graphs
  BTest.Res2 <- BTest.Res[BTest.Res[, "Expected"] > 0, ]
  n.points <- length(BTest.Res2[, 1])
  x <- BTest.Res2[, "Z-value"]
  y11 <- BTest.Res2[, "CWRE Cnt"] / BTest.Res2[, "Expected"]
  y12 <- BTest.Res2[, "IWRE Cnt"] / BTest.Res2[, "Expected"]
  y13 <- BTest.Res2[, "WRES Cnt"] / BTest.Res2[, "Expected"]
  y2 <- BTest.Res[, "LB"] / BTest.Res[, "Expected"]
  y3 <- BTest.Res[, "UB"] / BTest.Res[, "Expected"]

  par(defpar)

  for (info in list(list(y13, "W", "w", "WRES"), list(y11, "C", "c", "CWRES"), list(y12, "I", "i", "IWRES"))) {
    yvals <- info[[1]]; ch_out <- info[[2]]; ch_in <- info[[3]]; lbl <- info[[4]]
    plot(0, 0, type = "n", xlim = c(-3, 3), ylim = c(0, max(yvals)),
         xlab = "Z-value", ylab = "Observed Count/Expected Count",
         main = paste("Ratio of", lbl, "Count to Expected Count"))
    for (i in 1:n.points) {
      if (yvals[i] < y2[i] | yvals[i] > y3[i]) text(x[i], yvals[i], ch_out, col = "red")
      else text(x[i], yvals[i], ch_in)
    }
    lines(x, yvals, lty = 3)
    lines(BTest.Res[, "Z-value"], y2, lty = 1)
    lines(BTest.Res[, "Z-value"], y3, lty = 1)
    abline(h = 1, lty = 2)
  }

  options(width = 95)
  PrinMTxt(capture.output(CWRES[abs(CWRES[, "WRES"]) > 3, c("ID", "TIME", "DV", "PRED", "IPRE", "WRES", "CWRE", "IWRE")]),
           Header1 = "Extreme WRES values larger than 3")
  PrinMTxt(capture.output(CWRES[abs(CWRES[, "WRES"]) > 2, c("ID", "TIME", "DV", "PRED", "IPRE", "WRES", "CWRE", "IWRE")]),
           Header1 = "Extreme WRES values larger than 2")
  PrinMTxt(capture.output(CWRES[abs(CWRES[, "CWRE"]) > 2, c("ID", "TIME", "DV", "PRED", "IPRE", "WRES", "CWRE", "IWRE")]),
           Header1 = "Extreme CWRE values larger than 2")
  PrinMTxt(capture.output(CWRES[abs(CWRES[, "IWRE"]) > 2, c("ID", "TIME", "DV", "PRED", "IPRE", "WRES", "CWRE", "IWRE")]),
           Header1 = "Extreme IWRE values larger than 2")

  # Individual residual curves
  ID2s <- unique(CWRES[, "ID"])
  nID2 <- length(ID2s)
  for (i in 1:nID2) {
    ID <- ID2s[i]
    par(defpar)
    par(oma = c(1, 1, 3, 1), mfrow = c(3, 2), lty = 1)

    x1 <- CWRES[CWRES[, "ID"] == ID, "PRED"]
    y1 <- CWRES[CWRES[, "ID"] == ID, "WRES"]
    ord1 <- order(x1)
    plot(x1[ord1], y1[ord1], type = "b", cex = 0.5, xlab = "PRED", ylab = "WRES")
    abline(h = 0, lty = 3)

    FITi <- FITDATA3[floor(FITDATA3[, "ID"]) == ID, ]
    x4 <- FITi[, "TALD"]
    y4 <- FITi[, "WRES"]
    DxPlotPost(x4, y4, mat = FITi, xlbl = "Time after Latest Dose", ylbl = "WRES", smooth = "T")

    x2 <- CWRES[CWRES[, "ID"] == ID, "IPRE"]
    y2 <- CWRES[CWRES[, "ID"] == ID, "CWRE"]
    ord2 <- order(x2)
    plot(x2[ord2], y2[ord2], type = "b", cex = 0.5, xlab = "IPRE", ylab = "CWRES")
    abline(h = 0, lty = 3)

    y5 <- FITi[, "CWRE"]
    DxPlotPost(x4, y5, mat = FITi, xlbl = "Time after Latest Dose", ylbl = "CWRES", smooth = "T")

    y3 <- CWRES[CWRES[, "ID"] == ID, "IWRE"]
    plot(x2[ord2], y3[ord2], type = "b", cex = 0.5, xlab = "IPRE", ylab = "IWRES")
    abline(h = 0, lty = 3)

    y6 <- FITi[, "IWRE"]
    DxPlotPost(x4, y6, mat = FITi, xlbl = "Time after Latest Dose", ylbl = "IWRES", smooth = "T")
    mtext(outer = FALSE, side = 1, line = 4, cex = 0.7, "*Fractional number means dosing occasion.")

    mtext(outer = TRUE, side = 3, paste("Individiual Residual Curve, ID =", ID2s[i]))

    cCWRE <- CWRES[CWRES[, "ID"] == ID, "CWRE"]
    resRun[i] <- run.test.nm(cCWRE)
    Col <- "black"
    if (resRun[i] < 0.05) Col <- "red"

    mtext(outer = TRUE, side = 3, line = -1, col = Col, cex = 0.8,
          paste("IOFV =", format(IDStat5[IDStat5[, "ID"] == ID, "iOFV"], digits = 6),
                ", OFV per DV =", format(IDStat5[IDStat5[, "ID"] == ID, "OFVpDV"], digits = 5),
                ", Run Test(CWRES) p =", format(resRun[i], digits = 3)))

    nCWRE <- length(cCWRE)
    if (nCWRE > 1) {
      cCWRE2 <- as.numeric(cCWRE[2:nCWRE] >= cCWRE[1:(nCWRE - 1)]) - 0.5
    } else {
      cCWRE2 <- 0.5
    }
    resRun2[i] <- run.test.nm(cCWRE2)
  }

  # Bad run test summary
  AddPage(Cex = 0.8)
  resRun <- cbind(ID2s, resRun)
  colnames(resRun) <- c("ID", "p.value")
  PrinTxt(1, 1, "Subjects with Bad Run Test Result on CWRES", Cex = 1)
  nBadRun <- sum(resRun[, "p.value"] < 0.05)

  if (nBadRun > 0) {
    BadRun <- resRun[resRun[, "p.value"] < 0.05, ]
    sBadRun <- capture.output(BadRun)
    for (i in 1:length(sBadRun)) {
      PrinTxt(3 + i, 3, sBadRun[i], Cex = 0.8)
    }
  } else {
    PrinTxt(4, 3, "None")
  }

  AddPage(Cex = 0.8)
  resRun2 <- cbind(ID2s, resRun2)
  colnames(resRun2) <- c("ID", "p.value")
  PrinTxt(1, 1, "Subjects with Bad Run Test 2 Result on CWRES", Cex = 1)
  nBadRun2 <- sum(resRun2[, "p.value"] < 0.05)

  if (nBadRun2 > 0) {
    BadRun2 <- resRun2[resRun2[, "p.value"] < 0.05, ]
    sBadRun2 <- capture.output(BadRun2)
    for (i in 1:length(sBadRun2)) {
      PrinTxt(3 + i, 3, sBadRun2[i], Cex = 0.8)
    }
  } else {
    PrinTxt(4, 3, "None")
  }

  ClosePDF()
  message("S4-Residuals.PDF generated.")
}
