#' Generate Prediction Diagnostics Report (S3-Predictions.PDF)
#'
#' Generates a PDF report with spaghetti plots, DV vs PRED/IPRE,
#' and individual fitting curves.
#'
#' @param run_dir character, path to the NONMEM run directory (default: current directory)
#' @param file character, output PDF filename
#' @param model character, run/model name; NULL auto-detects via GetCurModelName()
#' @return Invisibly, a list with the output \code{file}, number of \code{pages},
#'   compartment count, and subject count.
#' @export
nmw_report_pred <- function(run_dir = getwd(), file = "S3-Predictions.PDF", model = NULL) {
  owd <- setwd(run_dir); on.exit(setwd(owd))

  defpar <- par(no.readonly = TRUE); defpar$new <- NULL
  on.exit(par(defpar), add = TRUE)

  ## ---- parsing + statistics: reused verbatim from nmw ----------------------
  CtlName   <- if (is.null(model)) GetCurModelName() else model
  RunNumber <- CtlName

  fLogDV <- tryCatch(LogDV2(), error = function(e) FALSE)

  sdtab <- read.table("sdtab", header = TRUE, skip = 1, na.strings = "***********")
  nRec <- length(sdtab[, 1])
  FDATA <- ReadFDATA()
  FDATA[is.na(FDATA)] <- 0

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

  ## ---- spaghetti / goodness-of-fit pages (DxPlotPost) via sp_figure_page ---
  ## Each of the three original mfrow blocks becomes one full figure page.
  .plot_spaghetti <- function(doc, FITDATA2, RunNumber, fLogDV, cmt_label = "") {
    # Page 1: DV vs PRED / DV vs IPRE  (mfrow = c(2,1))
    sp_figure_page(doc, {
      DxPlotPost(FITDATA2$PRED, FITDATA2$DV, mat = FITDATA2, xlbl = "PRED", ylbl = "DV", smooth = "F")
      DxPlotPost(FITDATA2$IPRE, FITDATA2$DV, mat = FITDATA2, xlbl = "IPRE", ylbl = "DV", smooth = "F")
      mtext(outer = TRUE, side = 3, paste("Spaghetti Plot of Model", toupper(RunNumber), cmt_label))
    }, mfrow = c(2, 1), oma = c(1, 1, 3, 1), mar = c(4, 4, 2, 1))

    # Page 2: DV / PRED / IPRE vs Time after Latest Dose  (mfrow = c(3,1))
    sp_figure_page(doc, {
      DxPlotPost(FITDATA2$TALD, FITDATA2$DV,   mat = FITDATA2, xlbl = "Time after Latest Dose", ylbl = "DV",   smooth = "F")
      DxPlotPost(FITDATA2$TALD, FITDATA2$PRED, mat = FITDATA2, xlbl = "Time after Latest Dose", ylbl = "PRED", smooth = "F")
      DxPlotPost(FITDATA2$TALD, FITDATA2$IPRE, mat = FITDATA2, xlbl = "Time after Latest Dose", ylbl = "IPRE", smooth = "F")
      mtext(outer = TRUE, side = 3, paste("Spaghetti Plot of Model", toupper(RunNumber), cmt_label))
    }, mfrow = c(3, 1), oma = c(1, 1, 3, 1), mar = c(4, 4, 2, 1))

    # Page 3: same but log / original scale on Y  (mfrow = c(3,1))
    sp_figure_page(doc, {
      if (fLogDV == FALSE) {
        DxPlotPost(FITDATA2$TALD, FITDATA2$DV,   mat = FITDATA2, xlbl = "Time after Latest Dose", ylbl = "DV",   smooth = "F", Log = "y")
        DxPlotPost(FITDATA2$TALD, FITDATA2$PRED, mat = FITDATA2, xlbl = "Time after Latest Dose", ylbl = "PRED", smooth = "F", Log = "y")
        DxPlotPost(FITDATA2$TALD, FITDATA2$IPRE, mat = FITDATA2, xlbl = "Time after Latest Dose", ylbl = "IPRE", smooth = "F", Log = "y")
        mtext(outer = TRUE, side = 3, paste("Spaghetti Plot of Model", toupper(RunNumber), "(log scale on Y)", cmt_label))
      } else {
        DxPlotPost(FITDATA2$TALD, exp(FITDATA2$DV),   mat = FITDATA2, xlbl = "Time after Latest Dose", ylbl = "DV",   smooth = "F")
        DxPlotPost(FITDATA2$TALD, exp(FITDATA2$PRED), mat = FITDATA2, xlbl = "Time after Latest Dose", ylbl = "PRED", smooth = "F")
        DxPlotPost(FITDATA2$TALD, exp(FITDATA2$IPRE), mat = FITDATA2, xlbl = "Time after Latest Dose", ylbl = "IPRE", smooth = "F")
        mtext(outer = TRUE, side = 3, paste("Spaghetti Plot of Model", toupper(RunNumber), "(original scale on Y)", cmt_label))
      }
    }, mfrow = c(3, 1), oma = c(1, 1, 3, 1), mar = c(4, 4, 2, 1))
  }

  ## ---- per-subject individual fitting curve page (2 panels) ----------------
  ## Values are passed as arguments and force()d so the deferred figure expr
  ## captures this call's data, not the loop's final iteration.
  .plot_indiv <- function(doc, fLogDV, maxx, maxy1, miny2,
                          w1, x1, y1, z1, w2, x2, y2, z2,
                          DoseLbl, titleTxt) {
    force(fLogDV); force(maxx); force(maxy1); force(miny2)
    force(w1); force(x1); force(y1); force(z1)
    force(w2); force(x2); force(y2); force(z2)
    force(DoseLbl); force(titleTxt)
    sp_figure_page(doc, {
      # ---- Panel 1 (linear / log-scale-labelled) ----
      if (fLogDV == FALSE) {
        plot(w1, x1, type = "p", cex = 0.7, xlim = c(0, maxx), ylim = c(0, maxy1),
             xlab = "Time (Circle=DV, Dotted Line=PRED, Line=IPRE)", ylab = "DV,PRED,IPRE")
      } else {
        plot(w1, x1, type = "p", cex = 0.7, xlim = c(0, maxx), ylim = c(miny2, maxy1),
             xlab = "Time (Circle=DV, Dotted Line=PRED, Line=IPRE)", ylab = "DV,PRED,IPRE (log scale)")
      }
      lines(w1, y1, type = "b", cex = 0.3, col = "blue", lty = 3)
      lines(w1, z1, type = "b", cex = 0.5, col = "red")
      if (!is.null(DoseLbl) && nrow(DoseLbl) > 0) {
        for (r in seq_len(nrow(DoseLbl))) {
          text(DoseLbl$t[r], 0, DoseLbl$a[r], adj = c(0, 1), cex = 0.7)
        }
      }
      mtext(outer = FALSE, side = 1, line = 4, cex = 0.8, "*Numbers within graph are dosing amounts.")

      # ---- Panel 2 (log-y interval / original scale) ----
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

      mtext(outer = TRUE, side = 3, titleTxt)
    }, mfrow = c(2, 1), oma = c(1, 1, 3, 1), mar = c(4, 4, 2, 1))
  }

  ## ---- build the document --------------------------------------------------
  doc <- sp_new(file, paper = "letter", family = "Courier", size = 10)
  frame_set(doc, top = doc$H - 45, bottom = 45, left = 45, right = doc$W - 45)

  # Title / intro flow page (carries the running CONFIDENTIAL footer)
  flow_run(doc,
    list(
      block_para("Summary 3 - Predictions", size = 16, font = 2),
      block_rule(),
      block_para(sprintf("Model: %s", toupper(RunNumber)), size = 11),
      block_para(paste0("Spaghetti / goodness-of-fit plots (DV vs PRED, DV vs IPRE, ",
                        "concentration-time) and per-subject individual fitting curves."),
                 size = 10),
      block_para(sprintf("Number of compartments plotted: %d", n.CMT), size = 10)),
    footer = block_para("CONFIDENTIAL", size = 8, align = "center"))

  ## ---- spaghetti pages (guard empty CMT subsets) ---------------------------
  if (n.CMT == 1) {
    if (nrow(FITDATA) > 0) {
      .plot_spaghetti(doc, FITDATA, RunNumber, fLogDV)
    }
  } else {
    for (i in 1:n.CMT) {
      FITDATA2 <- FITDATA[FITDATA$CMT == CMTs[i], ]
      if (nrow(FITDATA2) == 0) next          # EDGE-CASE FIX: skip empty CMT
      .plot_spaghetti(doc, FITDATA2, RunNumber, fLogDV,
                      cmt_label = paste(", CMT=", CMTs[i], sep = ""))
    }
  }

  ## ---- individual fitting curves (verbatim extraction + empty guards) ------
  sdtab[, "ID"] <- FDATA[, "ID"]
  IDs <- unique(sdtab[sdtab$DV > 0, "ID"])
  nID <- length(IDs)

  if (!exists("CMTs")) CMTs <- 1

  for (i in seq_len(nID)) {
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

      # EDGE-CASE FIX: subject/CMT with NO observation rows -> skip before any
      # plotting/labelling (prevents "zero-length labels" in plot()/text()).
      if (length(x) == 0) next

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

      # Precompute dosing-amount labels (both ADDL-expanded and simple cases),
      # each guarded against empty dosing histories.
      DoseLbl <- NULL
      if (length(intersect(names(FDATA), "ADDL")) == 1) {
        DosingHist <- FDATA[FDATA[, "ID"] == cID & FDATA[, "AMT"] > 0 & FDATA[, "MDV"] == 1,
                            c("AMT", "TIME", "II", "ADDL")]
        nDoseRec <- nrow(DosingHist)
        if (nDoseRec > 0) {
          for (j in 1:nDoseRec) {
            cAMT <- DosingHist[j, "AMT"]; cTIME <- DosingHist[j, "TIME"]
            cII <- DosingHist[j, "II"];  cADDL <- DosingHist[j, "ADDL"]
            for (k in 0:cADDL) {
              DoseLbl <- rbind(DoseLbl, data.frame(t = cTIME + k * cII, a = cAMT))
            }
          }
        }
      } else if ("AMT" %in% colnames(FDATA)) {
        a <- FDATA[FDATA[, "ID"] == cID & FDATA[, "MDV"] == 1 & FDATA[, "AMT"] > 0, "AMT"]
        d <- FDATA[FDATA[, "ID"] == cID & FDATA[, "MDV"] == 1 & FDATA[, "AMT"] > 0, "TIME"]
        if (length(d) > 0) DoseLbl <- data.frame(t = d, a = a)
      }

      titleTxt <- if (n.CMT == 1) {
        paste("Individual Fitting Curve, ID=", cID, sep = "")
      } else {
        paste("Individual Fitting Curve, ID=", cID, ", CMT=", CMTs[m], sep = "")
      }

      .plot_indiv(doc, fLogDV, maxx, maxy1, miny2,
                  w1, x1, y1, z1, w2, x2, y2, z2, DoseLbl, titleTxt)
    }
  }

  sp_close(doc)
  message(paste0(file, " generated."))
  invisible(list(file = file, pages = doc$page_no, n.CMT = n.CMT, nID = nID))
}
