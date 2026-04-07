#' Generate Prediction Diagnostics Report (S3-Predictions.PDF)
#'
#' Generates a PDF report with spaghetti plots, DV vs PRED/IPRE,
#' and individual fitting curves.
#'
#' @param run_dir character, path to the NONMEM run directory (default: current directory)
#' @export
nmw_report_pred <- function(run_dir = getwd()) {
  owd <- setwd(run_dir)
  on.exit(setwd(owd))

  defpar <- par(no.readonly = TRUE)
  defpar$new <- NULL
  on.exit(par(defpar), add = TRUE)

  RunNumber <- GetCurModelName()

  fLogDV <- tryCatch(LogDV2(), error = function(e) FALSE)

  sdtab <- read.table("sdtab", header = TRUE, skip = 1, na.strings = "***********")
  nRec <- length(sdtab[, 1])
  FDATA <- ReadFDATA()
  FDATA[is.na(FDATA)] <- 0

  PrepPDF("S3-Predictions.PDF")

  if ("AMT" %in% colnames(FDATA)) {
    tald_result <- CalcTaLDForReport(FDATA, nRec)
    TaLD <- tald_result$TaLD
    NewID <- tald_result$NewID

    if ("CMT" %in% colnames(FDATA)) {
      FITDATA <- cbind(NewID, TaLD, sdtab[, c("CMT", "MDV", "DV", "PRED", "IPRE")])
      colnames(FITDATA) <- c("ID", "TALD", "CMT", "MDV", "DV", "PRED", "IPRE")
    } else {
      FITDATA <- cbind(NewID, TaLD, sdtab[, c("MDV", "DV", "PRED", "IPRE")])
      colnames(FITDATA) <- c("ID", "TALD", "MDV", "DV", "PRED", "IPRE")
    }
  } else {
    FITDATA <- sdtab
    FITDATA$TALD <- FITDATA$TIME
  }

  if (!is.null(FITDATA$MDV)) {
    FITDATA <- FITDATA[FITDATA$MDV == 0, ]
  }

  n.CMT <- 1
  if (!is.null(FITDATA$CMT)) {
    CMTs <- sort(unique(FITDATA$CMT))
    n.CMT <- length(CMTs)
  }

  .plot_spaghetti <- function(FITDATA2, RunNumber, fLogDV, cmt_label = "") {
    par(defpar)
    par(oma = c(1, 1, 3, 1), mfrow = c(2, 1), lty = 1)
    DxPlotPost(FITDATA2$PRED, FITDATA2$DV, mat = FITDATA2, xlbl = "PRED", ylbl = "DV", smooth = "F")
    DxPlotPost(FITDATA2$IPRE, FITDATA2$DV, mat = FITDATA2, xlbl = "IPRE", ylbl = "DV", smooth = "F")
    mtext(outer = TRUE, side = 3, paste("Spaghetti Plot of Model", toupper(RunNumber), cmt_label))

    par(oma = c(1, 1, 3, 1), mfrow = c(3, 1), lty = 1)
    DxPlotPost(FITDATA2$TALD, FITDATA2$DV, mat = FITDATA2, xlbl = "Time after Latest Dose", ylbl = "DV", smooth = "F")
    DxPlotPost(FITDATA2$TALD, FITDATA2$PRED, mat = FITDATA2, xlbl = "Time after Latest Dose", ylbl = "PRED", smooth = "F")
    DxPlotPost(FITDATA2$TALD, FITDATA2$IPRE, mat = FITDATA2, xlbl = "Time after Latest Dose", ylbl = "IPRE", smooth = "F")
    mtext(outer = TRUE, side = 3, paste("Spaghetti Plot of Model", toupper(RunNumber), cmt_label))

    par(oma = c(1, 1, 3, 1), mfrow = c(3, 1), lty = 1)
    if (fLogDV == FALSE) {
      DxPlotPost(FITDATA2$TALD, FITDATA2$DV, mat = FITDATA2, xlbl = "Time after Latest Dose", ylbl = "DV", smooth = "F", Log = "y")
      DxPlotPost(FITDATA2$TALD, FITDATA2$PRED, mat = FITDATA2, xlbl = "Time after Latest Dose", ylbl = "PRED", smooth = "F", Log = "y")
      DxPlotPost(FITDATA2$TALD, FITDATA2$IPRE, mat = FITDATA2, xlbl = "Time after Latest Dose", ylbl = "IPRE", smooth = "F", Log = "y")
      mtext(outer = TRUE, side = 3, paste("Spaghetti Plot of Model", toupper(RunNumber), "(log scale on Y)", cmt_label))
    } else {
      DxPlotPost(FITDATA2$TALD, exp(FITDATA2$DV), mat = FITDATA2, xlbl = "Time after Latest Dose", ylbl = "DV", smooth = "F")
      DxPlotPost(FITDATA2$TALD, exp(FITDATA2$PRED), mat = FITDATA2, xlbl = "Time after Latest Dose", ylbl = "PRED", smooth = "F")
      DxPlotPost(FITDATA2$TALD, exp(FITDATA2$IPRE), mat = FITDATA2, xlbl = "Time after Latest Dose", ylbl = "IPRE", smooth = "F")
      mtext(outer = TRUE, side = 3, paste("Spaghetti Plot of Model", toupper(RunNumber), "(original scale on Y)", cmt_label))
    }
  }

  if (n.CMT == 1) {
    .plot_spaghetti(FITDATA, RunNumber, fLogDV)
  } else {
    for (i in 1:n.CMT) {
      FITDATA2 <- FITDATA[FITDATA$CMT == CMTs[i], ]
      .plot_spaghetti(FITDATA2, RunNumber, fLogDV, cmt_label = paste(", CMT=", CMTs[i], sep = ""))
    }
  }

  # Individual fitting curves
  sdtab[, "ID"] <- FDATA[, "ID"]
  IDs <- unique(sdtab[sdtab$DV > 0, "ID"])
  nID <- length(IDs)

  par(defpar)

  if (!exists("CMTs")) CMTs <- 1

  for (i in 1:nID) {
    cID <- IDs[i]

    for (m in 1:n.CMT) {
      if ("CMT" %in% colnames(sdtab)) {
        x <- sdtab[sdtab[, "ID"] == cID & sdtab[, "MDV"] == 0 & sdtab[, "CMT"] == CMTs[m], "DV"]
        y <- sdtab[sdtab[, "ID"] == cID & sdtab[, "MDV"] == 0 & sdtab[, "CMT"] == CMTs[m], "PRED"]
        z <- sdtab[sdtab[, "ID"] == cID & sdtab[, "MDV"] == 0 & sdtab[, "CMT"] == CMTs[m], "IPRE"]
        w <- sdtab[sdtab[, "ID"] == cID & sdtab[, "MDV"] == 0 & sdtab[, "CMT"] == CMTs[m], "TIME"]
        tx <- sdtab[sdtab[, "MDV"] == 0 & sdtab[, "CMT"] == CMTs[m], "DV"]
        ty <- sdtab[sdtab[, "MDV"] == 0 & sdtab[, "CMT"] == CMTs[m], "PRED"]
        tz <- sdtab[sdtab[, "MDV"] == 0 & sdtab[, "CMT"] == CMTs[m], "IPRE"]
        maxy <- max(tx, ty, tz)
      } else {
        x <- sdtab[sdtab[, "ID"] == cID & sdtab[, "MDV"] == 0, "DV"]
        y <- sdtab[sdtab[, "ID"] == cID & sdtab[, "MDV"] == 0, "PRED"]
        z <- sdtab[sdtab[, "ID"] == cID & sdtab[, "MDV"] == 0, "IPRE"]
        w <- sdtab[sdtab[, "ID"] == cID & sdtab[, "MDV"] == 0, "TIME"]
        tx <- sdtab[sdtab[, "MDV"] == 0, "DV"]
        ty <- sdtab[sdtab[, "MDV"] == 0, "PRED"]
        tz <- sdtab[sdtab[, "MDV"] == 0, "IPRE"]
        maxy <- max(tx, ty, tz)
      }

      if (fLogDV == FALSE) {
        w1 <- w[x > 0]; x1 <- x[x > 0]; y1 <- y[x > 0]; z1 <- z[x > 0]
        w2 <- w[x > 0 & y > 0 & z > 0]; x2 <- x[x > 0 & y > 0 & z > 0]
        y2 <- y[x > 0 & y > 0 & z > 0]; z2 <- z[x > 0 & y > 0 & z > 0]
      } else {
        w2 <- w1 <- w; x2 <- x1 <- x; y2 <- y1 <- y; z2 <- z1 <- z
      }

      if (length(x1) == 0 || length(x2) == 0) next

      maxx <- 1.05 * max(w)
      maxy1 <- max(x1, y1, z1)
      miny2 <- min(x2, y2, z2)

      par(oma = c(1, 1, 3, 1), mfrow = c(2, 1), lty = 1)

      if (fLogDV == FALSE) {
        plot(w1, x1, type = "p", cex = 0.7, xlim = c(0, maxx), ylim = c(0, maxy1),
             xlab = "Time (Circle=DV, Dotted Line=PRED, Line=IPRE)", ylab = "DV,PRED,IPRE")
      } else {
        plot(w1, x1, type = "p", cex = 0.7, xlim = c(0, maxx), ylim = c(miny2, maxy1),
             xlab = "Time (Circle=DV, Dotted Line=PRED, Line=IPRE)", ylab = "DV,PRED,IPRE (log scale)")
      }
      lines(w1, y1, type = "b", cex = 0.3, col = "blue", lty = 3)
      lines(w1, z1, type = "b", cex = 0.5, col = "red")

      if (length(intersect(names(FDATA), "ADDL")) == 1) {
        DosingHist <- FDATA[FDATA[, "ID"] == cID & FDATA[, "AMT"] > 0 & FDATA[, "MDV"] == 1,
                            c("AMT", "TIME", "II", "ADDL")]
        nDoseRec <- length(DosingHist[, "AMT"])
        for (j in 1:nDoseRec) {
          cAMT <- DosingHist[j, "AMT"]; cTIME <- DosingHist[j, "TIME"]
          cII <- DosingHist[j, "II"]; cADDL <- DosingHist[j, "ADDL"]
          for (k in 0:cADDL) text(cTIME + k * cII, 0, cAMT, adj = c(0, 1), cex = 0.7)
        }
      } else if ("AMT" %in% colnames(FDATA)) {
        a <- FDATA[FDATA[, "ID"] == cID & FDATA[, "MDV"] == 1 & FDATA[, "AMT"] > 0, "AMT"]
        d <- FDATA[FDATA[, "ID"] == cID & FDATA[, "MDV"] == 1 & FDATA[, "AMT"] > 0, "TIME"]
        for (j in 1:length(d)) text(d[j], 0, a[j], adj = c(0, 1), cex = 0.7)
      }

      mtext(outer = FALSE, side = 1, line = 4, cex = 0.8, "*Numbers within graph are dosing amounts.")

      if (fLogDV == FALSE) {
        plot(w2, x2, type = "p", cex = 0.7, xlim = c(0, maxx), ylim = c(miny2, maxy1),
             xlab = "Time (Circle=DV, Dotted Line=PRED, Line=IPRE)", ylab = "DV,PRED,IPRE (log interval)", log = "y")
        lines(w2, y2, type = "b", cex = 0.3, col = "blue", lty = 3)
        lines(w2, z2, type = "b", cex = 0.5, col = "red")
      } else {
        plot(w2, exp(x2), type = "p", cex = 0.7, xlim = c(0, maxx), ylim = c(exp(miny2), exp(maxy1)),
             xlab = "Time (Circle=DV, Dotted Line=PRED, Line=IPRE)", ylab = "DV,PRED,IPRE (original scale)")
        lines(w2, exp(y2), type = "b", cex = 0.3, col = "blue", lty = 3)
        lines(w2, exp(z2), type = "b", cex = 0.5, col = "red")
      }

      if (n.CMT == 1) {
        mtext(outer = TRUE, side = 3, paste("Individiual Fitting Curve, ID=", cID, sep = ""))
      } else {
        mtext(outer = TRUE, side = 3, paste("Individiual Fitting Curve, ID=", cID, ", CMT=", CMTs[m], sep = ""))
      }
    }
  }

  ClosePDF()
  message("S3-Predictions.PDF generated.")
}
