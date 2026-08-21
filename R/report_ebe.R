#' Generate EBE Diagnostics Report (S5-EBE.PDF)
#'
#' Generates a PDF report with empirical Bayes estimate (EBE) analysis
#' including ETA distributions, shrinkage, covariate relationships,
#' correlations, and multiple linear regression with influence diagnostics.
#'
#' @param run_dir character, path to the NONMEM run directory (default: current directory)
#' @param file character, output PDF filename
#' @param model character, run/model name; NULL auto-detects via GetCurModelName()
#' @return Invisibly, a list with the output \code{file}, number of \code{pages},
#'   eta count, and subject count.
#' @export
nmw_report_ebe <- function(run_dir = getwd(), file = "S5-EBE.PDF", model = NULL) {
  owd <- setwd(run_dir); on.exit(setwd(owd))
  oldwarn  <- options(warn = -1); on.exit(options(oldwarn),  add = TRUE)
  oldwidth <- options(width = 95); on.exit(options(oldwidth), add = TRUE)

  CtlName   <- if (is.null(model)) GetCurModelName() else model
  RunNumber <- CtlName

  ## ==== parsing + statistics: reused VERBATIM from report_ebe.R ==============
  FDATA <- ReadFDATA()
  FDATA[is.na(FDATA)] <- 0
  VarStat <- NMVarStat(FDATA)

  XML <- readLines(paste0(RunNumber, ".xml"))
  PHI <- ReadLastTable(paste0(RunNumber, ".phi"))   # chained-$EST safe (last table)
  EXT <- ReadLastTable(paste0(RunNumber, ".ext"))
  EXT <- EXT[EXT[, "ITERATION"] >= 0, ]

  params  <- CountEXTParams(EXT)
  nEtaAll <- params$nEta
  nEpsAll <- params$nEps

  OMEGA   <- BtwTagMat("omega",   XML, nEtaAll)
  OMEGAse <- BtwTagMat("omegase", XML, nEtaAll)

  OMa  <- rbind(OMEGA, OMEGAse)
  nEta <- length(OMa[1, ])
  OM   <- OMa[1:nEta, , drop = FALSE]
  SeOM <- OMa[(nEta + 1):(2 * nEta), , drop = FALSE]

  EtaNames    <- character()
  EtaNamesRaw <- character()
  EtaErrorRaw <- character()
  for (i in 1:nEta) {
    EtaNames    <- c(EtaNames,    paste("Eta", i))
    EtaNamesRaw <- c(EtaNamesRaw, paste0("ETA.", i, "."))
    EtaErrorRaw <- c(EtaErrorRaw, paste0("ETC.", i, ".", i, "."))
  }
  rownames(OM)   <- EtaNames; colnames(OM)   <- EtaNames
  rownames(SeOM) <- EtaNames; colnames(SeOM) <- EtaNames

  sdtab <- read.table("sdtab", header = TRUE, skip = 1, na.strings = "***********")
  if (nrow(sdtab) == nrow(FDATA) & length(unique(sdtab$ID)) != length(unique(FDATA$ID))) {
    sdtab[, "ID"] <- FDATA[, "ID"]
    warning("\nCheck sdtab's ID column! Consider using FORMAT=s1PE16.8\n")
  }
  IDs <- unique(sdtab[, "ID"])
  nID <- length(IDs)

  EtaNames2  <- "ID"
  seEtaNames <- character()
  for (i in 1:nEta) {
    EtaNames2  <- c(EtaNames2,  paste0("ETA", i))
    seEtaNames <- c(seEtaNames, paste0("seETA", i))
  }

  tabEta <- PHI[, c("ID", EtaNamesRaw)]
  colnames(tabEta) <- EtaNames2

  if (var(tabEta[, "ETA1"]) == 0) {
    patab  <- read.table("patab", header = TRUE, skip = 1, na.strings = "***********")
    tabEta <- unique(patab[, EtaNames2])
  }

  seEtaRaw <- PHI[, c("ID", EtaErrorRaw)]
  colnames(seEtaRaw) <- c("ID", seEtaNames)
  seEta <- PHI[, EtaErrorRaw, drop = FALSE]
  colnames(seEta) <- seEtaNames

  EtaRep <- cbind(tabEta, PHI[, EtaErrorRaw, drop = FALSE])

  ## ---- covariate classification + augmentation (verbatim) ------------------
  cc       <- ClassifyCovariates(VarStat)
  ResWords <- cc$ResWords
  RmvList  <- cc$RmvList
  CatList  <- cc$CatList
  ContList <- cc$ContList

  CovariateList <- union(CatList, ContList)
  tabCovariate  <- matrix(nrow = nID, ncol = length(CovariateList))
  colnames(tabCovariate) <- CovariateList
  nRec <- dim(FDATA)[1]

  if (length(CovariateList) > 0) {
    PrevID <- 0; j <- 1
    for (i in 1:nRec) {
      if (PrevID != FDATA[i, "ID"]) {
        for (k in 1:length(CovariateList)) {
          tabCovariate[j, CovariateList[k]] <- FDATA[i, CovariateList[k]]
        }
        j <- j + 1
      }
      PrevID <- FDATA[i, "ID"]
    }
  }

  tabEta <- cbind(tabEta, tabCovariate)
  VarEta <- names(tabEta)
  Cat1   <- intersect(VarEta, CatList)
  Cont1  <- intersect(VarEta, ContList)
  nCat   <- length(Cat1)
  nCont  <- length(Cont1)

  Cat2Rmv <- NULL
  if (nCat > 0) {
    for (i in 1:nEta) {
      for (j in 1:nCat) {
        if (length(unique(tabEta[, Cat1[j]])) == 1) Cat2Rmv <- union(Cat2Rmv, Cat1[j])
      }
    }
  }
  Cat1 <- setdiff(Cat1, Cat2Rmv)
  nCat <- length(Cat1)

  ## ---- correlation / covariance matrices (verbatim) ------------------------
  corEtaA <- cor(tabEta[, c(Cont1, EtaNames2[-1])])
  covEta2 <- cov(tabEta[, EtaNames2[-1], drop = FALSE])
  corEta2 <- cov2cor(covEta2)
  corOM   <- cov2cor(OM)

  oAll <- list(corEtaA, covEta2, OM, covEta2 / OM, corEta2, corOM, corEta2 / corOM)
  names(oAll) <- c("Correlation of Covariates and EBE", "Covariance of EBE", "Omega Matrix",
                   "Ratios of Cov(EBE)/OM", "Correlation of EBE", "Correlation from Omega Matrix",
                   "Ratios of Cor(EBE)/(Cor from OM)")

  ## ---- multiple linear regression + influence diagnostics (verbatim) -------
  tabEta2 <- RemoveNA(tabEta)

  EtaZero <- rep(FALSE, nEta)
  for (j in 1:nEta) {
    EtaZero[j] <- (all.equal(tabEta[, paste0("ETA", j)], rep(0, nrow(tabEta))) == TRUE)
  }

  CovarList <- c(Cat1, Cont1)

  ESQ_SSE_THRESHOLD   <- 0.15
  DFFITS_THRESHOLD    <- 1
  COVRATIO_MULTIPLIER <- 3

  regr_ok <- (length(CovarList) > 0 & !all(EtaZero))
  result  <- NULL
  result2 <- NULL
  if (regr_ok) {
    modelstr <- character()
    for (i in 1:length(CovarList)) {
      modelstr <- paste(modelstr, CovarList[i])
      if (i < length(CovarList)) modelstr <- paste(modelstr, "+")
    }

    result  <- list()
    result2 <- list()
    for (j in 1:nEta) {
      if (!EtaZero[j]) {
        formstr    <- paste0("ETA", j, "~", modelstr)
        result[j]  <- list(summary(lm(formula(formstr), data = tabEta)))
        result2[j] <- list(mlr2(tabEta2[, paste0("ETA", j)], tabEta2[, CovarList]))
      } else {
        result[j]  <- NA
        result2[j] <- NA
      }
    }
    names(result)  <- EtaNames
    names(result2) <- EtaNames
  }

  ## ---- individual ETA confidence intervals / shrinkage (verbatim) ----------
  LL   <- matrix(nrow = nID, ncol = nEta)
  UL   <- matrix(nrow = nID, ncol = nEta)
  ZERO <- matrix(nrow = nID, ncol = nEta)
  RSE  <- matrix(nrow = nID, ncol = nEta)
  SHR  <- matrix(nrow = nID, ncol = nEta)

  LLnames <- character(); ULnames <- character()
  ZEROnames <- character(); RSEnames <- character(); SHRnames <- character()
  for (i in 1:nEta) {
    LLnames   <- c(LLnames,   paste0("LL", i))
    ULnames   <- c(ULnames,   paste0("UL", i))
    ZEROnames <- c(ZEROnames, paste0("ZERO", i))
    RSEnames  <- c(RSEnames,  paste0("RSE", i))
    SHRnames  <- c(SHRnames,  paste0("SHR", i))
  }
  colnames(LL) <- LLnames; colnames(UL) <- ULnames
  colnames(ZERO) <- ZEROnames; colnames(RSE) <- RSEnames; colnames(SHR) <- SHRnames

  for (i in 1:nEta) {
    LL[, i]   <- EtaRep[, i + 1] - 2 * EtaRep[, nEta + i + 1]
    UL[, i]   <- EtaRep[, i + 1] + 2 * EtaRep[, nEta + i + 1]
    ZERO[, i] <- (LL[, i] * UL[, i] < 0)
    RSE[, i]  <- EtaRep[, nEta + i + 1] / abs(EtaRep[, i + 1])
    SHR[, i]  <- EtaRep[, nEta + i + 1] / sqrt(OM[i, i])
  }

  write.csv(cbind(EtaRep, LL, UL, ZERO, RSE, SHR), file = "EBEALL.CSV", row.names = FALSE)

  ## ==== TEXT / TABLE pages: simPDF blocks (replaces PrinMTxt/AddPage) ========
  blocks <- list(
    block_para("Summary 5 - Empirical Bayes Estimates (EBE)", size = 16, font = 2),
    block_rule(),
    block_para(sprintf("Model: %s    Subjects: %d    Etas: %d", RunNumber, nID, nEta), size = 10)
  )

  ## Section: Estimation vs EBE (correlation / covariance / omega matrices)
  blocks <- c(blocks, list(
    block_spacer(4),
    block_para("Estimation vs EBE", size = 13, font = 2)))
  for (nm in names(oAll)) {
    blocks <- c(blocks, list(
      block_keep(list(
        block_para(nm, size = 11, font = 2),
        block_matrix(round(oAll[[nm]], 5), size = 9))),
      block_spacer(4)))
  }

  ## Section: Multiple linear regression + influence (text dumps)
  if (regr_ok) {
    blocks <- c(blocks, list(
      block_spacer(4),
      block_para("Multiple Linear Regression", size = 13, font = 2)))
    for (i in 1:nEta) {
      if (!EtaZero[i]) {
        sOut  <- capture.output(result[[i]])
        sOut2 <- capture.output(result2[[i]])
        blocks <- c(blocks, list(
          block_keep(list(
            block_para(paste("Multiple Linear Regression : ETA", i), size = 11, font = 2),
            block_pre(sOut[4:length(sOut)], size = 8))),
          block_spacer(2),
          block_para(paste("Multiple Linear Regression - Influence : ETA", i), size = 11, font = 2),
          block_pre(sOut2, size = 8),
          block_spacer(6)))
      }
    }
  }

  ## Section: Individual ETA confidence intervals (per-eta tables)
  blocks <- c(blocks, list(
    block_spacer(4),
    block_para("Individual ETA Confidence Intervals", size = 13, font = 2),
    block_para("LL/UL: Estimate +/- 2*SE   ZERO: CI covers 0   RSE: SE/|ETA|   SHR: SE/sqrt(OMEGA)",
               size = 8)))
  for (i in 1:nEta) {
    tabHeader <- c("ID", paste0("ETA", i), paste0("seETA", i), paste0("LL", i),
                   paste0("UL", i), paste0("ZERO", i), paste0("RSE", i), paste0("SHR", i))
    tabEtaOut <- cbind(EtaRep[, c(1, i + 1, i + nEta + 1)], LL[, i], UL[, i], ZERO[, i], RSE[, i], SHR[, i])
    colnames(tabEtaOut) <- tabHeader
    rownames(tabEtaOut) <- NULL
    blocks <- c(blocks, list(
      block_keep(list(block_para(paste("ETA", i), size = 11, font = 2))),
      block_table(as.data.frame(tabEtaOut), size = 8, family = "mono"),
      block_spacer(6)))
  }

  ## ==== emit document ========================================================
  doc <- sp_new(file, paper = "letter", family = "Courier", size = 10)
  frame_set(doc, top = doc$H - 45, bottom = 45, left = 45, right = doc$W - 45)
  flow_run(doc, Filter(Negate(is.null), blocks),
           footer = block_para("CONFIDENTIAL", size = 8, align = "center"))

  ## ==== FIGURE pages: sp_figure_page grids (replace per-section par(mfrow)) ==

  ## (A) ETA distributions: PlotDistribution() = 3 panels/eta, 4 etas/page.
  RowPerPage <- 4
  eta_chunks <- split(seq_len(nEta), ceiling(seq_len(nEta) / RowPerPage))
  for (ch in eta_chunks) {
    sp_figure_page(doc, {
      for (ii in ch) {
        var.data <- EtaRep[, paste0("ETA", ii)]
        var.data <- var.data[var.data != 0]
        PlotDistribution(var.data, paste("Eta", ii),
                         show_shrinkage = TRUE, omega_var = OM[ii, ii],
                         show_ttest = TRUE)
      }
      npad <- (RowPerPage - length(ch)) * 3
      if (npad > 0) for (q in seq_len(npad)) plot.new()
      mtext(outer = TRUE, side = 3, "Normality and Population Shrinkage of Etas")
    }, mfrow = c(RowPerPage, RowPerPage - 1), oma = c(1, 1, 3, 1), mar = c(4, 4, 2, 1))
  }

  ## (B) ETA vs categorical covariates: nEta x nCatChunk boxplot grids.
  if (nCat > 0) {
    MAX_CELLS_PER_PAGE <- 35
    max_cats_per_page  <- max(1, floor(MAX_CELLS_PER_PAGE / nEta))
    nPages_cat <- ceiling(nCat / max_cats_per_page)
    chunks <- split(seq_len(nCat), ceiling(seq_len(nCat) / max_cats_per_page))
    for (pg in seq_along(chunks)) {
      chunk     <- chunks[[pg]]
      nCatChunk <- length(chunk)
      title_str <- "ETA vs Categorical Variables"
      if (nPages_cat > 1) title_str <- paste0(title_str, " (page ", pg, "/", nPages_cat, ")")
      sp_figure_page(doc, {
        for (i in 1:nEta) {
          for (j_idx in seq_along(chunk)) {
            j    <- chunk[j_idx]
            form <- paste0("ETA", i, "~", Cat1[j])
            boxplot(formula = formula(form), data = tabEta,
                    xlab = Cat1[j], ylab = paste("ETA", i))
            Anova.pvalue <- anova(lm(formula(form), data = tabEta))[[5]][1]
            if (is.finite(Anova.pvalue) & Anova.pvalue < 0.05) { txtcol <- "red" }
            else { txtcol <- "black" }
            mtext(side = 3, cex = 0.6, col = txtcol,
                  paste("ANOVA p =", format(Anova.pvalue, digits = 4)))
          }
        }
        mtext(outer = TRUE, side = 3, title_str)
      }, mfrow = c(nEta, nCatChunk), oma = c(1, 1, 3, 1), mar = c(4, 4, 2, 2))
    }
  }

  ## (C) Covariate vs ETA pair plot (single full-page pairs()).
  sp_figure_page(doc, {
    EBEpair(tabEta[, c("ID", Cont1, EtaNames2[-1])], RunNumber)
  }, mfrow = c(1, 1), oma = c(0, 0, 2, 0), mar = c(4, 4, 2, 1))

  ## (D) ETA vs continuous covariates grid (verbatim inner drawing).
  IDdep         <- rownames(VarStat[VarStat[, "DepID"] == 1, ])
  IDTimeMDVdep  <- rownames(VarStat[VarStat[, "DepIDTimeMDV"] == 1, ])
  IDTimeEVIDdep <- rownames(VarStat[VarStat[, "DepIDTimeEVID"] == 1, ])
  OnlyIDTime0   <- union(IDTimeMDVdep, IDTimeEVIDdep)
  OnlyIDTime1   <- setdiff(OnlyIDTime0, IDdep)
  OnlyIDTime2   <- setdiff(OnlyIDTime1, ResWords)
  OUTLIER_PCT   <- 0.025

  sp_figure_page(doc, {
    for (i in 1:nEta) {
      for (j in 0:nCont) {
        if (j == 0) {
          plot(0, 0, type = "n", ylim = c(0, 1), xlim = c(0, 1), xaxt = "n", yaxt = "n",
               ylab = "", xlab = "", bty = "n")
          text(0, 0.5, paste0("ETA", i), adj = 0)
        } else {
          Xaxt <- ifelse(i < nEta, "n", "s")
          Xlab <- ifelse(i < nEta, "", Cont1[j])
          Yaxt <- ifelse(j == 1, "s", "n")
          Ylab <- ifelse(j == 1, paste("ETA", i), "")

          x <- tabEta[, Cont1[j]]
          minx <- min(x, na.rm = TRUE); maxx <- max(x, na.rm = TRUE)
          y <- tabEta[, paste0("ETA", i)]

          IDList <- tabEta[, "ID"]
          n <- length(x)
          ordx <- order(x); ordy <- order(y)

          OutSide <- union(union(ordx[1:max(1, round(OUTLIER_PCT * n))], ordx[round((1 - OUTLIER_PCT) * n):n]),
                           union(ordy[1:max(1, round(OUTLIER_PCT * n))], ordy[round((1 - OUTLIER_PCT) * n):n]))
          InSide <- intersect(ordx[max(2, round(OUTLIER_PCT * n)):min(round((1 - OUTLIER_PCT) * n), n - 1)],
                              ordy[max(2, round(OUTLIER_PCT * n)):min(round((1 - OUTLIER_PCT) * n), n - 1)])
          plot(x, y, type = "n", xlim = c(minx - 0.1 * (maxx - minx), maxx + 0.1 * (maxx - minx)),
               xaxt = Xaxt, xlab = Xlab, yaxt = Yaxt, ylab = Ylab)
          text(x[OutSide], y[OutSide], IDList[OutSide], cex = 0.8)
          points(x[InSide], y[InSide], cex = 0.5)

          if (length(intersect(Cont1[j], OnlyIDTime2)) == 1) {
            for (k in 1:nID) {
              cID <- IDs[k]
              x2 <- unique(FDATA[FDATA[, "ID"] == cID, Cont1[j]])
              nx2 <- length(x2)
              y2 <- rep(tabEta[tabEta[, "ID"] == cID, paste0("ETA", i)], nx2)
              lines(x2, y2, lty = 1)
              points(x2[-1], y2[-1], cex = 0.5, pch = 2)
            }
          }
          abline(h = 0, lty = 2)
          lines(lowess(x, y), col = "red", lty = 1)
        }
      }
    }

    plot(0, 0, type = "n", ylim = c(0, 1), xlim = c(0, 1), xaxt = "n", yaxt = "n",
         ylab = "", xlab = "", bty = "n")
    if (nCont > 0) {
      for (j in 1:nCont) {
        plot(0, 0, type = "n", ylim = c(0, 1), xlim = c(0, 1), xaxt = "n", yaxt = "n",
             ylab = "", xlab = "", bty = "n")
        text(0.5, 0.8, Cont1[j])
      }
    }
  }, mfrow = c(nEta + 1, nCont + 1), oma = c(1, 1, 3, 1), mar = c(0, 0, 0, 0))

  ## (E) Influence diagnostics per eta (2x2 per page, verbatim inner drawing).
  if (regr_ok) {
    for (i in 1:nEta) {
      if (!EtaZero[i]) {
        sp_figure_page(doc, {
          D        <- result2[[i]][[2]][, "Cook's D"]
          y.hat    <- result2[[i]][[2]][, "Yhat"]
          sdr      <- result2[[i]][[2]][, "R-Student"]
          h        <- result2[[i]][[2]][, "hat"]
          e        <- result2[[i]][[2]][, "Residual"]
          SSE      <- result2[[i]][[6]]
          COVRATIO <- result2[[i]][[2]][, "COV-Ratio"]
          DFFITS   <- result2[[i]][[2]][, "DFFITS"]
          n        <- result2[[i]][[3]]
          p        <- result2[[i]][[4]]

          plot(D, type = "n", xlab = "Index", ylab = "Cook's Distance")
          for (jj in 1:n) {
            if (D[jj] == max(D)) text(jj, D[jj], jj)
            else points(jj, D[jj])
          }

          plot(y.hat, sdr, type = "n", xlab = "Predicted Value", ylab = "Studentized deleted residuals")
          for (jj in 1:n) {
            if (abs(sdr[jj]) > 2) text(y.hat[jj], sdr[jj], jj)
            else points(y.hat[jj], sdr[jj])
          }

          plot(h, e^2 / SSE, type = "n", xlab = "hat", ylab = "e^2/SSE")
          for (jj in 1:n) {
            if (e[jj]^2 / SSE > ESQ_SSE_THRESHOLD) text(h[jj], e[jj]^2 / SSE, jj)
            else points(h[jj], e[jj]^2 / SSE)
          }

          if (all(is.finite(COVRATIO))) {
            plot(COVRATIO, DFFITS, type = "n", xlab = "COVRATIO", ylab = "DFFITS")
            for (jj in 1:n) {
              if (abs(DFFITS[jj]) > DFFITS_THRESHOLD | abs(COVRATIO[jj] - 1) > COVRATIO_MULTIPLIER * p / n) {
                text(COVRATIO[jj], DFFITS[jj], jj)
              } else {
                points(COVRATIO[jj], DFFITS[jj])
              }
            }
          } else {
            plot.new()
          }

          mtext(paste("Influence Diagnostics on Eta", i), outer = TRUE, side = 3)
        }, mfrow = c(2, 2), oma = c(1, 1, 3, 1), mar = c(4, 4, 2, 1))
      }
    }
  }

  np <- doc$page_no
  sp_close(doc)
  message(paste(file, "generated."))
  invisible(list(file = file, pages = np, nEta = nEta, nID = nID))
}
