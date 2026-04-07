# Internal helper functions shared across report_*.R files
# These functions reduce code duplication in report generation.


#' Count THETA, ETA, EPS Parameters from EXT Data Frame
#'
#' @param EXT data.frame from .ext file (already read)
#' @return list with nTheta, nEta, nEps
#' @keywords internal
CountEXTParams <- function(EXT) {
  cols <- colnames(EXT)
  nTheta <- sum(startsWith(cols, "THETA"))
  nOmega <- sum(startsWith(cols, "OMEGA"))
  nSigma <- sum(startsWith(cols, "SIGMA"))
  nEta <- as.integer((sqrt(1 + 8 * nOmega) - 1) / 2)
  nEps <- as.integer((sqrt(1 + 8 * nSigma) - 1) / 2)
  list(nTheta = nTheta, nEta = nEta, nEps = nEps)
}


#' Merge IDStat with Individual OFV from PHI
#'
#' Common pattern used in report_ofv.R and report_resid.R.
#'
#' @param FDATA data.frame of NONMEM data
#' @param PHI data.frame from .phi file
#' @return list with IDStat3 (merged), VarStat, IDStat
#' @keywords internal
MergeIDStatOFV <- function(FDATA, PHI) {
  VarStat <- NMVarStat(FDATA)
  IDStat <- NMIDStat(FDATA)
  tOFV <- PHI[, c("ID", "OBJ")]
  colnames(tOFV) <- c("ID", "iOFV")
  IDStat2 <- cbind(rownames(IDStat[[1]]), IDStat[[1]])
  IDStat2Name <- names(IDStat2)
  IDStat2Name[1] <- "ID"
  names(IDStat2) <- IDStat2Name
  IDStat3 <- merge(tOFV, IDStat2)
  list(IDStat3 = IDStat3, VarStat = VarStat, IDStat = IDStat)
}


#' Calculate TaLD (Time after Latest Dose) for Report Generation
#'
#' Shared logic between report_pred.R and report_resid.R.
#'
#' @param FDATA data.frame of NONMEM data
#' @param nRec integer, number of records
#' @return list with ToLD, TaLD, DOCC, NewID vectors
#' @keywords internal
CalcTaLDForReport <- function(FDATA, nRec) {
  LastID <- 0
  ToLD <- vector(length = nRec)
  TaLD <- vector(length = nRec)
  DOCCvec <- vector(length = nRec)
  NewID <- vector(length = nRec)

  pToLD <- 0
  DOCC <- -1
  DosingHist2 <- NULL

  for (i in 1:nRec) {
    CurID <- FDATA[i, "ID"]
    if (LastID != CurID) {
      pToLD <- 0
      DOCC <- -1
      if (length(intersect(names(FDATA), "ADDL")) == 1) {
        DosingHist <- FDATA[FDATA[, "ID"] == CurID & FDATA[, "AMT"] > 0 & FDATA[, "MDV"] == 1,
                            c("AMT", "TIME", "II", "ADDL")]
        nDoseRec <- length(DosingHist[, "AMT"])
        for (j in 1:nDoseRec) {
          cADDL <- DosingHist[j, "ADDL"]
          if (cADDL > 0) {
            cAMT <- DosingHist[j, "AMT"]
            cTIME <- DosingHist[j, "TIME"]
            cII <- DosingHist[j, "II"]
            for (k in 1:cADDL) {
              DosingHist <- rbind(DosingHist, c(cAMT, cTIME + k * cII, NA, NA))
            }
          }
        }
        DosingHist2 <- DosingHist[order(DosingHist[, "TIME"]), c("AMT", "TIME")]
      }
    }

    if (FDATA[i, "AMT"] > 0 & FDATA[i, "MDV"] == 1) {
      cToLD <- FDATA[i, "TIME"]
      DOCC <- DOCC + 1
    } else {
      if (length(intersect(names(FDATA), "ADDL")) == 1 && !is.null(DosingHist2)) {
        CurTime <- FDATA[i, "TIME"]
        matching <- DosingHist2[DosingHist2[, "TIME"] < CurTime, "TIME"]
        cToLD <- if (length(matching) > 0) max(matching) else pToLD
      } else {
        cToLD <- pToLD
      }
    }
    ToLD[i] <- cToLD
    TaLD[i] <- FDATA[i, "TIME"] - cToLD
    NewID[i] <- FDATA[i, "ID"] + max(0, DOCC) / 100
    DOCCvec[i] <- DOCC
    LastID <- CurID
    pToLD <- cToLD
  }

  list(ToLD = ToLD, TaLD = TaLD, DOCC = DOCCvec, NewID = NewID)
}


#' Plot Distribution Panel (Summary + Histogram + QQ)
#'
#' Common 3-panel distribution plot used in report_ebe.R,
#' report_indipk.R, and report_resid.R.
#'
#' @param var.data numeric vector of data
#' @param var.name character, variable name for labeling
#' @param show_shrinkage logical, whether to show shrinkage
#' @param omega_var numeric, omega variance for shrinkage calc
#' @param show_ttest logical, whether to show t-test p-value
#' @param show_cv logical, whether to show CV percent
#' @keywords internal
PlotDistribution <- function(var.data, var.name, show_shrinkage = FALSE,
                             omega_var = NULL, show_ttest = FALSE,
                             show_cv = FALSE) {
  sdvar <- sd(var.data)

  # Panel 1: Summary statistics
  plot(0, 0, type = "n", ylim = c(1, 8), xlim = c(0, 10), xaxt = "n", yaxt = "n",
       ylab = "", xlab = "", bty = "n")
  sData <- summary(var.data)
  text(0, 8, paste("Minimum :", sData[1]), adj = 0)
  text(0, 7, paste("1st Qu. :", sData[2]), adj = 0)
  text(0, 6, paste("Median  :", sData[3]), adj = 0)
  text(0, 5, paste("Mean    :", sData[4]), adj = 0)
  text(0, 4, paste("3rd Qu. :", sData[5]), adj = 0)
  text(0, 3, paste("Maximum :", sData[6]), adj = 0)
  text(0, 2, paste("Std Dev :", format(sdvar, digits = 4)), adj = 0)
  if (show_ttest) {
    ttest.pvalue <- format(t.test(var.data)$p.value, digits = 4)
    text(0, 1, paste("t-test p=", ttest.pvalue), adj = 0)
  }
  if (show_cv) {
    text(0, 1, paste("CV (%)  :", format(sdvar / sData[4] * 100, digits = 4)), adj = 0)
  }
  mtext(side = 3, line = 1, var.name, adj = 0)

  # Panel 2: Histogram with density overlay
  h.res <- hist(var.data, plot = FALSE)
  d.res <- density(var.data, na.rm = TRUE)
  h.rat <- max(h.res$counts) / max(h.res$density)
  xrange <- range(d.res$x)
  yrange <- c(0, max(h.res$counts, max(d.res$y) * h.rat))
  plot(h.res, xlim = xrange, ylim = yrange, xlab = var.name, main = "")
  lines(d.res$x, h.rat * d.res$y)

  if (sdvar > 0 && length(var.data) < 5000) {
    p.value <- format(shapiro.test(var.data)$p.value, digits = 4)
    mtext(paste("Shapiro-Wilk test p-value =", p.value), cex = 0.7, line = 0, adj = 0)
  }

  if (show_shrinkage && !is.null(omega_var) && omega_var > 0) {
    shrink <- sdvar / sqrt(omega_var)
    mtext(paste("Shrinkage =", format(shrink, digits = 4)), cex = 0.7, line = 1, adj = 0)
  }

  # Panel 3: QQ plot
  qqnorm(var.data, main = "", ylab = var.name, cex = 0.4)
}


#' OFV Report Split Screen Layout
#'
#' Layout matrix for the S1-OFV report's split screen display.
#'
#' @keywords internal
OFV_SCREEN_LAYOUT <- matrix(
  c(0, 1, 0, 1,
    0.55, 0.97, 0.6, 0.95,
    0.55, 0.97, 0.45, 0.8,
    0.55, 0.97, 0.1, 0.45),
  nrow = 4, ncol = 4, byrow = TRUE
)


#' Read EXT File from Current Directory
#'
#' @return data.frame or NULL if file not found
#' @keywords internal
ReadEXTFile <- function() {
  EXTName <- paste0(GetCurModelName(), ".ext")
  if (!file.exists(EXTName)) return(NULL)
  read.table(EXTName, skip = 1, header = TRUE)
}
