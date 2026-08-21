#' Generate Residual Diagnostics Report (S4-Residuals.PDF)
#'
#' Generates a PDF report with WRES, CWRES, IWRES diagnostics including
#' scatter plots, histograms, normality tests, binomial count tests,
#' and individual residual curves with run tests.
#'
#' @param run_dir character, path to the NONMEM run directory (default: current directory)
#' @param file character, output PDF filename
#' @param model character, run/model name; NULL auto-detects via GetCurModelName()
#' @return Invisibly, a list with the output \code{file}, number of \code{pages},
#'   subject and bad-run counts, and usable/omitted residual names.
#' @export
nmw_report_resid <- function(run_dir = getwd(), file = "S4-Residuals.PDF", model = NULL) {
  owd <- setwd(run_dir); on.exit(setwd(owd))

  ## ---- parsing + statistics: reused verbatim from nmw report_resid.R --------
  RunNumber <- if (is.null(model)) GetCurModelName() else model

  PHI <- ReadLastTable(paste0(RunNumber, ".phi"))   # chained-$EST safe (last table)
  XML <- readLines(paste0(RunNumber, ".xml"))

  nRec <- GetNRecFromXML(XML)

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

  # Some runs omit one or more residual columns. Preserve a fixed data schema,
  # but record which residuals contain usable values so every downstream
  # statistic and plot can omit unavailable diagnostics explicitly.
  requiredCols <- c("MDV", "TIME", "DV", "PRED", "IPRE")
  missingRequired <- setdiff(requiredCols, colnames(sdtab))
  if (length(missingRequired) > 0) {
    stop("sdtab is missing required column(s): ",
         paste(missingRequired, collapse = ", "))
  }

  residualSource <- c(WRES = "WRES", CWRE = "CWRES", IWRE = "IWRE")
  obsRows <- !is.na(sdtab$MDV) & sdtab$MDV == 0
  availableResiduals <- names(residualSource)[vapply(
    residualSource,
    function(cn) {
      cn %in% colnames(sdtab) && any(is.finite(sdtab[obsRows, cn]))
    },
    logical(1)
  )]
  omittedResiduals <- setdiff(names(residualSource), availableResiduals)

  for (cn in setdiff(unname(residualSource), colnames(sdtab))) {
    sdtab[[cn]] <- NA_real_
  }

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

  ## ---- empty-data guard (mirrors report_pred's defensive style) -------------
  doc <- simPDF::sp_new(file, paper = "letter", family = "Courier", size = 10)
  simPDF::frame_set(doc, top = doc$H - 45, bottom = 45, left = 45, right = doc$W - 45)
  footer <- simPDF::block_para("CONFIDENTIAL", size = 8, align = "center")

  if (nrow(CWRES) == 0 || nrow(FITDATA3) == 0) {
    simPDF::flow_run(doc, list(
      simPDF::block_para("Summary 4 - Residuals", size = 16, font = 2),
      simPDF::block_para(paste("Model:", toupper(RunNumber)), size = 10),
      simPDF::block_rule(),
      simPDF::block_para("No observation records (MDV == 0) found; nothing to plot.",
                         size = 11)), footer = footer)
    np <- doc$page_no
    simPDF::sp_close(doc)
    message(paste(file, "generated (no observations)."))
    return(invisible(list(file = file, pages = np, nID = 0L,
                          nBadRun = NA_integer_, nBadRun2 = NA_integer_,
                          residuals = character(),
                          omittedResiduals = names(residualSource))))
  }

  if (length(availableResiduals) == 0) {
    simPDF::flow_run(doc, list(
      simPDF::block_para("Summary 4 - Residuals", size = 16, font = 2),
      simPDF::block_para(paste("Model:", toupper(RunNumber)), size = 10),
      simPDF::block_rule(),
      simPDF::block_para(
        "No usable WRES, CWRES, or IWRES values were found; nothing to plot.",
        size = 11)), footer = footer)
    np <- doc$page_no
    simPDF::sp_close(doc)
    message(paste(file, "generated (no usable residuals)."))
    return(invisible(list(file = file, pages = np,
                          nID = length(unique(CWRES$ID)),
                          nBadRun = NA_integer_, nBadRun2 = NA_integer_,
                          residuals = character(),
                          omittedResiduals = names(residualSource))))
  }

  # Run the count test separately for each residual so occasional missing
  # values do not distort another residual's expected counts.
  BTests <- setNames(lapply(availableResiduals, function(vname) {
    keep <- is.finite(CWRES[, vname])
    ResTest(CWRES[keep, , drop = FALSE], TestCols = vname)
  }), availableResiduals)

  # Per-subject run-test statistics (lifted verbatim out of the old plot loop
  # so the bad-run summary tables can be emitted with the other text)
  ID2s <- unique(CWRES[, "ID"])
  nID2 <- length(ID2s)
  resRun <- rep(NA_real_, nID2)
  resRun2 <- rep(NA_real_, nID2)
  if ("CWRE" %in% availableResiduals) {
    for (i in seq_len(nID2)) {
      ID <- ID2s[i]
      cCWRE <- CWRES[CWRES[, "ID"] == ID, "CWRE"]
      cCWRE <- cCWRE[is.finite(cCWRE)]
      nCWRE <- length(cCWRE)
      if (nCWRE > 0) {
        resRun[i] <- run.test.nm(cCWRE)
        if (nCWRE > 1) {
          cCWRE2 <- as.numeric(cCWRE[2:nCWRE] >= cCWRE[1:(nCWRE - 1)]) - 0.5
        } else {
          cCWRE2 <- 0.5
        }
        resRun2[i] <- run.test.nm(cCWRE2)
      }
    }
  }

  resRunM <- cbind(ID2s, resRun)
  colnames(resRunM) <- c("ID", "p.value")
  nBadRun <- if ("CWRE" %in% availableResiduals) {
    sum(resRunM[, "p.value"] < 0.05, na.rm = TRUE)
  } else {
    NA_integer_
  }

  resRun2M <- cbind(ID2s, resRun2)
  colnames(resRun2M) <- c("ID", "p.value")
  nBadRun2 <- if ("CWRE" %in% availableResiduals) {
    sum(resRun2M[, "p.value"] < 0.05, na.rm = TRUE)
  } else {
    NA_integer_
  }

  ## ---- drawing helpers (base-R plot expressions, run at draw time) ----------
  residualLabels <- c(WRES = "WRES", CWRE = "CWRES", IWRE = "IWRES")
  residualX <- c(WRES = "PRED", CWRE = "IPRE", IWRE = "IPRE")

  plotUnavailable <- function(label) {
    plot.new()
    text(0.5, 0.5, paste(label, "not available"), cex = 0.9)
  }

  # Count-ratio graph for one residual type.
  drawCountPlot <- function(BTest, vname) {
    usable <- is.finite(BTest[, "Expected"]) & BTest[, "Expected"] > 0
    tab <- BTest[usable, , drop = FALSE]
    x <- tab[, "Z-value"]
    yvals <- tab[, paste(vname, "Cnt")] / tab[, "Expected"]
    lower <- tab[, "LB"] / tab[, "Expected"]
    upper <- tab[, "UB"] / tab[, "Expected"]
    ylim <- range(c(0, yvals, lower, upper), finite = TRUE)
    if (diff(ylim) == 0) ylim <- ylim + c(-0.5, 0.5)
    ch_out <- substr(residualLabels[[vname]], 1, 1)
    ch_in <- tolower(ch_out)

    plot(0, 0, type = "n", xlim = c(-3, 3), ylim = ylim,
         xlab = "Z-value", ylab = "Observed Count/Expected Count",
         main = paste("Ratio of", residualLabels[[vname]], "Count to Expected Count"))
    for (i in seq_along(yvals)) {
      if (yvals[i] < lower[i] || yvals[i] > upper[i]) {
        text(x[i], yvals[i], ch_out, col = "red")
      }
      else text(x[i], yvals[i], ch_in)
    }
    lines(x, yvals, lty = 3)
    lines(x, lower, lty = 1)
    lines(x, upper, lty = 1)
    abline(h = 1, lty = 2)
  }

  makeCountPlotBlock <- function(BTest, vname) {
    force(BTest)
    force(vname)
    simPDF::block_plot(drawCountPlot(BTest, vname), height = 230)
  }

  drawResidualPair <- function(data, timeData, vname) {
    xname <- residualX[[vname]]
    label <- residualLabels[[vname]]
    ok <- is.finite(data[, xname]) & is.finite(data[, vname])
    if (any(ok)) {
      x <- data[ok, xname]
      y <- data[ok, vname]
      ord <- order(x)
      plot(x[ord], y[ord], type = "b", cex = 0.5,
           xlab = xname, ylab = label)
      abline(h = 0, lty = 3)
    } else {
      plotUnavailable(label)
    }

    okTime <- is.finite(timeData[, "TALD"]) & is.finite(timeData[, vname])
    if (any(okTime)) {
      DxPlotPost(timeData[okTime, "TALD"], timeData[okTime, vname],
                 mat = timeData[okTime, , drop = FALSE],
                 xlbl = "Time after Latest Dose", ylbl = label,
                 smooth = if (sum(okTime) > 1) "T" else "F")
    } else {
      plotUnavailable(label)
    }
  }

  # One full per-subject residual page. Missing diagnostics use explanatory
  # placeholders so the stable 3 x 2 layout is retained.
  drawSubject <- function(i) {
    ID <- ID2s[i]
    data <- CWRES[CWRES[, "ID"] == ID, , drop = FALSE]
    FITi <- FITDATA3[floor(FITDATA3[, "ID"]) == ID, ]
    for (vname in names(residualSource)) {
      if (vname %in% availableResiduals) {
        drawResidualPair(data, FITi, vname)
      } else {
        plotUnavailable(residualLabels[[vname]])
        plotUnavailable(residualLabels[[vname]])
      }
    }
    mtext(outer = FALSE, side = 1, line = 4, cex = 0.7, "*Fractional number means dosing occasion.")

    mtext(outer = TRUE, side = 3, paste("Individual Residual Curve, ID =", ID2s[i]))

    Col <- "black"
    if (is.finite(resRun[i]) && resRun[i] < 0.05) Col <- "red"
    runText <- if (is.finite(resRun[i])) {
      paste("Run Test(CWRES) p =", format(resRun[i], digits = 3))
    } else {
      "Run Test(CWRES): not available"
    }
    mtext(outer = TRUE, side = 3, line = -1, col = Col, cex = 0.8,
          paste("IOFV =", format(IDStat5[IDStat5[, "ID"] == ID, "iOFV"], digits = 6),
                ", OFV per DV =", format(IDStat5[IDStat5[, "ID"] == ID, "OFVpDV"], digits = 5),
                ",", runText))
  }

  # Round only the DISPLAY copy of an extreme-value subset (never the values fed
  # to any computation -- CWRES itself is untouched).
  fmtSub <- function(sub) {
    for (cn in c("PRED", "IPRE", "WRES", "CWRE", "IWRE")) {
      if (cn %in% colnames(sub)) sub[, cn] <- round(sub[, cn], 4)
    }
    sub
  }
  ## ---- text / tables via one flow_run() ------------------------------------
  blocks <- list(
    simPDF::block_para("Summary 4 - Residuals", size = 16, font = 2),
    simPDF::block_para(paste("Model:", toupper(RunNumber)), size = 10),
    simPDF::block_rule(),
    simPDF::block_para(
      paste("Usable residuals:",
            paste(residualLabels[availableResiduals], collapse = ", ")),
      size = 9),
    if (length(omittedResiduals) > 0) simPDF::block_para(
      paste("Omitted residuals (missing or without finite values):",
            paste(residualLabels[omittedResiduals], collapse = ", ")),
      size = 9, font = 3),
    simPDF::block_spacer(6),
    simPDF::block_para("Tests of residual counts using binomial distributions",
                       size = 13, font = 2))

  for (vname in availableResiduals) {
    blocks <- c(blocks, list(
      simPDF::block_keep(list(
        simPDF::block_para(residualLabels[[vname]], size = 11, font = 2),
        simPDF::block_matrix(round(BTests[[vname]], 3), size = 8))),
      simPDF::block_spacer(4),
      simPDF::block_para("Ratio of Observed to Expected Residual Counts",
                         size = 11, font = 2),
      simPDF::block_para("(letter outside the confidence band is drawn red)", size = 8),
      makeCountPlotBlock(BTests[[vname]], vname),
      simPDF::block_spacer(6)))
  }

  # Extreme residual value tables
  ExtCols <- c("ID", "TIME", "DV", "PRED", "IPRE", availableResiduals)
  thresholds <- list(WRES = c(3, 2), CWRE = 2, IWRE = 2)
  extremeDefs <- list()
  for (vname in availableResiduals) {
    for (threshold in thresholds[[vname]]) {
      keep <- is.finite(CWRES[, vname]) & abs(CWRES[, vname]) > threshold
      extremeDefs[[length(extremeDefs) + 1L]] <- list(
        sub = CWRES[keep, ExtCols, drop = FALSE],
        hdr = paste("Extreme", residualLabels[[vname]],
                    "values larger than", threshold))
    }
  }

  blocks <- c(blocks, list(simPDF::block_spacer(6)))
  for (d in extremeDefs) {
    blocks <- c(blocks, list(simPDF::block_para(d$hdr, size = 12, font = 2)))
    if (nrow(d$sub) > 0) {
      blocks <- c(blocks, list(simPDF::block_table(fmtSub(d$sub), size = 8, family = "mono")))
    } else {
      blocks <- c(blocks, list(simPDF::block_para("None", size = 10)))
    }
    blocks <- c(blocks, list(simPDF::block_spacer(4)))
  }

  # Bad run test summaries require CWRES.
  if ("CWRE" %in% availableResiduals) {
    blocks <- c(blocks, list(
      simPDF::block_spacer(6),
      simPDF::block_para("Subjects with Bad Run Test Result on CWRES", size = 13, font = 2),
      if (nBadRun > 0)
        simPDF::block_table(round(resRunM[resRunM[, "p.value"] < 0.05, , drop = FALSE], 5),
                            size = 9, family = "mono")
      else simPDF::block_para("None", size = 10),
      simPDF::block_spacer(6),
      simPDF::block_para("Subjects with Bad Run Test 2 Result on CWRES", size = 13, font = 2),
      if (nBadRun2 > 0)
        simPDF::block_table(round(resRun2M[resRun2M[, "p.value"] < 0.05, , drop = FALSE], 5),
                            size = 9, family = "mono")
      else simPDF::block_para("None", size = 10)))
  } else {
    blocks <- c(blocks, list(
      simPDF::block_spacer(6),
      simPDF::block_para("CWRES is unavailable; subject run tests were omitted.",
                         size = 10, font = 3)))
  }

  simPDF::flow_run(doc, Filter(Negate(is.null), blocks), footer = footer)

  ## ---- multi-panel diagnostic plot pages via sp_figure_page ----------------
  drawOverviewPair <- function(vname, absolute = FALSE) {
    label <- residualLabels[[vname]]
    if (!(vname %in% availableResiduals)) {
      plotUnavailable(label)
      plotUnavailable(label)
      return(invisible(NULL))
    }

    xname <- residualX[[vname]]
    y <- FITDATA3[, vname]
    if (absolute) {
      y <- abs(y)
      label <- paste0("|", label, "|")
    }

    ok <- is.finite(FITDATA3[, xname]) & is.finite(y)
    if (any(ok)) {
      yLabel <- if (absolute) {
        label
      } else {
        paste(label, "SD =", format(sd(y[ok]), digits = 4))
      }
      mat <- FITDATA3[ok, , drop = FALSE]
      DxPlotPost(FITDATA3[ok, xname], y[ok], mat = mat,
                 xlbl = xname, ylbl = yLabel,
                 smooth = if (sum(ok) > 1) "T" else "F")
    } else {
      plotUnavailable(label)
    }

    okTime <- is.finite(FITDATA3[, "TALD"]) & is.finite(y)
    if (any(okTime)) {
      mat <- FITDATA3[okTime, , drop = FALSE]
      DxPlotPost(FITDATA3[okTime, "TALD"], y[okTime], mat = mat,
                 xlbl = "Time after latest dose", ylbl = label,
                 smooth = if (sum(okTime) > 1) "T" else "F")
    } else {
      plotUnavailable(label)
    }
    invisible(NULL)
  }

  drawDistribution <- function(vname) {
    label <- residualLabels[[vname]]
    if (!(vname %in% availableResiduals)) {
      plotUnavailable(label)
      plotUnavailable(label)
      return(invisible(NULL))
    }

    var.data <- FITDATA3[, vname]
    var.data <- var.data[is.finite(var.data)]
    if (length(var.data) == 0) {
      plotUnavailable(label)
      plotUnavailable(label)
      return(invisible(NULL))
    }

    h.res <- hist(var.data, plot = FALSE)
    canDensity <- length(var.data) >= 2 && sd(var.data) > 0
    if (canDensity) {
      d.res <- density(var.data)
      h.rat <- max(h.res$counts) / max(h.res$density)
      xrange <- range(d.res$x)
      yrange <- c(0, max(h.res$counts, max(d.res$y) * h.rat))
      plot(h.res, xlim = xrange, ylim = yrange,
           xlab = paste(label, "SD =", format(sd(var.data), digits = 4)),
           main = "")
      lines(d.res$x, h.rat * d.res$y)
    } else {
      plot(h.res, xlab = paste(label, "SD = NA"), main = "")
    }

    if (length(var.data) >= 3 && length(var.data) <= 5000 && sd(var.data) > 0) {
      p.value <- format(shapiro.test(var.data)$p.value, digits = 4)
      mtext(paste("Shapiro-Wilk test p-value =", p.value),
            cex = 0.7, line = 0, adj = 0)
    }
    qqnorm(var.data, main = "", ylab = label, cex = 0.4)
    invisible(NULL)
  }

  # Residual scatter plots (3 x 2 grid)
  simPDF::sp_figure_page(doc, {
    for (vname in names(residualSource)) drawOverviewPair(vname)
    mtext(outer = TRUE, side = 3, paste("Residuals of Model", toupper(RunNumber)))
  }, mfrow = c(3, 2), oma = c(1, 1, 3, 1))

  # Absolute residuals (3 x 2 grid)
  simPDF::sp_figure_page(doc, {
    for (vname in names(residualSource)) drawOverviewPair(vname, absolute = TRUE)
    mtext(outer = TRUE, side = 3, paste("Absolute Values of Residuals of Model", toupper(RunNumber)))
  }, mfrow = c(3, 2), oma = c(1, 1, 3, 1))

  # Histograms + QQ plots (3 x 2 grid: one row per residual type)
  simPDF::sp_figure_page(doc, {
    for (vname in names(residualSource)) drawDistribution(vname)
    mtext(outer = TRUE, side = 3, paste("Distribution of Residuals of Model", toupper(RunNumber)))
  }, mfrow = c(3, 2), oma = c(1, 1, 3, 1))

  # Individual residual curves (one 3 x 2 page per subject)
  for (i in seq_len(nID2)) {
    simPDF::sp_figure_page(doc, drawSubject(i), mfrow = c(3, 2), oma = c(1, 1, 3, 1))
  }

  np <- doc$page_no
  simPDF::sp_close(doc)
  message(paste(file, "generated."))
  invisible(list(file = file, pages = np, nID = nID2,
                 nBadRun = nBadRun, nBadRun2 = nBadRun2,
                 residuals = availableResiduals,
                 omittedResiduals = omittedResiduals))
}
