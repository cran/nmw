#' Post-Hoc OFV and CWRES Calculation
#'
#' Calculates individual OFV and CWRES using METHOD=COND or METHOD=ZERO.
#'
#' @param ctlname character, control file name (without extension)
#' @param sdtab data.frame, standard output table
#' @return list with OFV, OFVi, CWRES
#' @keywords internal
OBJfpost <- function(ctlname, sdtab) {
  FCON <- readLines("FCON")
  iESTM <- 0
  for (i in 1:length(FCON)) {
    if (substring(FCON[i], 1, 4) == "ESTM") iESTM <- i
  }
  ESTM <- read.fwf("FCON", rep(4, 20), skip = iESTM - 1, nrow = 1)
  fMETH <- ESTM[1, 9] + 1

  TDAT <- sdtab
  TDAT <- TDAT[TDAT[, "MDV"] == 0, ]
  IDs <- unique(TDAT[, "ID"])

  OM <- as.matrix(read.table("fort.53"))
  NETA <- length(OM[1, ])
  OM <- OM[1:NETA, ]

  if (file.exists("fort.54")) {
    SG <- as.matrix(read.table("fort.54"))
    NEPS <- length(SG[1, ])
    SG <- SG[1:NEPS, ]
  } else {
    NEPS <- 1
    SG <- matrix(1)
  }

  etanames <- "ID"
  for (i in 1:NETA) {
    etanames <- c(etanames, paste("ETA", i, sep = ""))
  }

  if (fMETH == 0) {
    ETAr <- matrix(rep(0, NETA * length(IDs)), ncol = NETA)
    ETA <- cbind(IDs, ETAr)
    colnames(ETA) <- c(etanames)
  } else {
    ETAs <- read.table(paste(ctlname, ".ETA", sep = ""), header = TRUE, skip = 1)
    ETA <- as.matrix(ETAs[, etanames])
  }

  gnames <- vector()
  for (i in 1:NETA) gnames <- append(gnames, paste("G", i, "1", sep = ""))
  hnames <- vector()
  for (i in 1:NEPS) hnames <- append(hnames, paste("H", i, "1", sep = ""))

  CWRE <- vector()
  OFVi <- cbind(IDs, rep(NA, length(IDs)))

  DATcolnames <- colnames(TDAT)
  nDATcol <- length(DATcolnames)

  for (i in 1:length(IDs)) {
    DATi <- matrix(as.matrix(TDAT[TDAT[, "ID"] == IDs[i], ]), ncol = nDATcol)
    colnames(DATi) <- DATcolnames
    Yi <- DATi[, "DV"]
    F0i <- DATi[, "PRED"]
    F1i <- DATi[, "IPRE"]
    Gi <- matrix(as.matrix(DATi[, gnames]), ncol = NETA)
    Hi <- matrix(as.matrix(DATi[, hnames]), ncol = NEPS)
    RES0i <- Yi - F0i
    RES1i <- Yi - F1i
    COVi <- Gi %*% OM %*% t(Gi) + diag(diag(Hi %*% SG %*% t(Hi)), nrow = length(Yi))

    if (fMETH == 0) {
      WRESi <- SqrtInvCov(COVi) %*% RES0i
    } else {
      ETAi <- ETA[ETA[, "ID"] == IDs[i], 2:(NETA + 1)]
      WRESi <- SqrtInvCov(COVi) %*% (RES1i + Gi %*% ETAi)
    }
    OFVi[i, 2] <- determinant(COVi, logarithm = TRUE)$modulus[1] + t(WRESi) %*% WRESi
    CWRE <- append(CWRE, WRESi)
  }
  output <- list(sum(OFVi[, 2]), OFVi, cbind(TDAT, CWRE))
  names(output) <- list("OFV", "OFVi", "CWRES")
  return(output)
}


#' Check if Model Uses Log-Transformed DV (from Control File)
#'
#' @param CtlFileName character, path to NONMEM control file
#' @return logical, TRUE if LOG(F) is found after $PRED or $ERROR
#' @export
LogDV <- function(CtlFileName) {
  ModFile <- toupper(readLines(CtlFileName))

  iPRED <- grep("$PRED", ModFile, fixed = TRUE)
  iERROR <- grep("$ERROR", ModFile, fixed = TRUE)
  iLOGF <- grep("LOG(F)", ModFile, fixed = TRUE)

  if (length(iPRED) == 1) {
    if (iLOGF > iPRED) return(TRUE)
  }

  if (length(iERROR) == 1) {
    if (iLOGF > iERROR) return(TRUE)
  }
  return(FALSE)
}


#' Check if Model Uses Log-Transformed DV (from FSUBS)
#'
#' @return logical, TRUE if LOG(F) is found after SUBROUTINE PRED or ERROR in FSUBS
#' @export
LogDV2 <- function() {
  if (!file.exists("FSUBS")) return(FALSE)
  ModFile <- readLines("FSUBS")

  iPRED <- grep("SUBROUTINE PRED (", ModFile, fixed = TRUE)
  iERROR <- grep("SUBROUTINE ERROR (", ModFile, fixed = TRUE)
  iLOGF <- grep("LOG(F)", ModFile, fixed = TRUE)

  if (length(iLOGF) == 1) {
    if (length(iPRED) == 1) {
      if (iLOGF > iPRED) return(TRUE)
    }

    if (length(iERROR) == 1) {
      if (iLOGF > iERROR) return(TRUE)
    }
  }

  return(FALSE)
}
