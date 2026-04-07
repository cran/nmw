#' Generate Parameter Estimates Report (S2-Parameters.PDF)
#'
#' Generates a PDF report with theta, omega, and sigma parameter estimates,
#' standard errors, confidence intervals, and significance flags.
#'
#' @param run_dir character, path to the NONMEM run directory (default: current directory)
#' @export
nmw_report_param <- function(run_dir = getwd()) {
  owd <- setwd(run_dir)
  on.exit(setwd(owd))

  defpar <- par(no.readonly = TRUE)
  defpar$new <- NULL
  on.exit(par(defpar), add = TRUE)

  CtlName <- GetCurModelName()

  XML <- readLines(paste0(CtlName, ".xml"))
  EXT <- read.table(paste0(CtlName, ".ext"), skip = 1, header = TRUE)
  EXT <- EXT[EXT[, "ITERATION"] >= 0, ]

  params <- CountEXTParams(EXT)
  nThetaAll <- params$nTheta
  nEtaAll <- params$nEta
  nEpsAll <- params$nEps

  THETA <- as.double(BtwTagVals("nm:theta", XML))
  THETASE <- as.double(BtwTagVals("nm:thetase", XML))
  OMEGA <- BtwTagMat("omega", XML, nEtaAll)
  OMEGAse <- BtwTagMat("omegase", XML, nEtaAll)
  SIGMA <- BtwTagMat("sigma", XML, nEpsAll)
  SIGMAse <- BtwTagMat("sigmase", XML, nEpsAll)

  Thetas <- cbind(THETA, THETASE)
  OMa <- rbind(OMEGA, OMEGAse)
  SGa <- rbind(SIGMA, SIGMAse)

  nThAll <- length(Thetas[, 1])
  Fixed <- vector()
  Unfixed <- vector()

  for (i in 1:nThAll) {
    if (Thetas[i, 2] == 1e+10) {
      Fixed <- c(Fixed, i)
    } else {
      Unfixed <- c(Unfixed, i)
    }
  }
  nFixedTh <- length(Fixed)
  nUnfixedTh <- length(Unfixed)
  ThRowName <- character()
  for (i in 1:nThAll) {
    ThRowName <- c(ThRowName, paste("Theta", i))
  }
  rownames(Thetas) <- ThRowName
  colnames(Thetas) <- c("Point Estitmate", "Standard Error")
  LL <- Thetas[, 1] - 2 * Thetas[, 2]
  UL <- Thetas[, 1] + 2 * Thetas[, 2]
  ZERO <- Thetas[, 2] / abs(Thetas[, 1]) > 0.5
  ONE <- (Thetas[, 1] - 2 * Thetas[, 2] - 1) * (Thetas[, 1] + 2 * Thetas[, 2] - 1) < 0 |
    (Thetas[, 1] - 2 * Thetas[, 2] + 1) * (Thetas[, 1] + 2 * Thetas[, 2] + 1) < 0

  Thetas <- cbind(Thetas, LL, UL, ZERO, ONE)
  UnfixedThetas <- Thetas[Unfixed, ]

  nEta <- length(OMa[1, ])
  OM <- OMa[1:nEta, , drop = FALSE]
  SeOM <- OMa[(nEta + 1):(2 * nEta), , drop = FALSE]
  EtaNames <- character()

  for (i in 1:nEta) {
    EtaNames <- c(EtaNames, paste("Eta", i))
  }
  rownames(OM) <- EtaNames
  colnames(OM) <- EtaNames
  rownames(SeOM) <- EtaNames
  colnames(SeOM) <- EtaNames

  RSEOM <- SeOM / abs(OM) * 100

  for (i in 1:nEta) {
    for (j in i:nEta) {
      if (j > i) OM[i, j] <- OM[i, j] / sqrt(OM[i, i] * OM[j, j])
    }
  }

  nEps <- length(SGa[1, ])
  SG <- SGa[1:nEps, ]
  SeSG <- SGa[(nEps + 1):(2 * nEps), ]

  # --- PDF Generation ---
  PrepPDF("S2-Parameters.PDF")

  AddPage()
  PrinTxt(1, 1, "Summary 2 - Parameters", Cex = 1.2)
  PrinTxt(3, 1, "Thetas", Cex = 1.0)
  PrinTxt(5, 3, paste("Number of All Thetas     :", nThAll))
  PrinTxt(6, 3, paste("Number of Fixed Thetas   :", nFixedTh))
  PrinTxt(7, 3, paste("Number of Unfixed Thetas :", nUnfixedTh))

  if (nFixedTh > 0) {
    PrinTxt(9, 2, "Fixed Theta Values", Cex = 0.9)
    for (i in 1:nFixedTh) {
      PrinTxt(9 + i, 5, paste("Theta", Fixed[i], ":", Thetas[Fixed[i], 1]))
    }
  }

  PrinTxt(9 + nFixedTh + 2, 2, "Estimated Thetas", Cex = 0.9)
  sUnfixed <- capture.output(UnfixedThetas)
  for (i in 1:length(sUnfixed)) {
    PrinTxt(9 + nFixedTh + 2 + i, 5, sUnfixed[i])
  }

  PrinTxt(nThAll + 13.5, 6, "*LL  : Lower Limit", Cex = 0.7)
  PrinTxt(nThAll + 14, 6, " UL  : Upper Limit", Cex = 0.7)
  PrinTxt(nThAll + 14.5, 6, " ZERO: Is this maybe zero? 0:No, 1:Yes", Cex = 0.7)
  PrinTxt(nThAll + 15, 6, " ONE : Is this maybe one?  0:No, 1:Yes", Cex = 0.7)

  AddPage()
  PrinTxt(3, 1, "Omegas", Cex = 1.0)
  PrinTxt(5, 3, paste("Number of Etas           :", nEta))

  PrinTxt(7, 2, "Omega Matrix", Cex = 0.9)
  sOM <- capture.output(OM)
  for (i in 1:length(sOM)) {
    PrinTxt(8 + i, 5, sOM[i])
  }
  PrinTxt(nEta + 10.5, 6, "*Lower triangle is covariance matrix.", Cex = 0.7)
  PrinTxt(nEta + 11, 6, " Upper triangle is correlation matrix.", Cex = 0.7)
  PrinTxt(nEta + 11.5, 6, " Diagonal elements are variances.", Cex = 0.7)

  PrinTxt(nEta + 13, 3, "Interindividual Variability (CV) in case of exp(eta) model (x100)")
  for (i in 1:nEta) {
    PrinTxt(nEta + 14, i * 8, paste("Eta", i))
    PrinTxt(nEta + 15, i * 8, format(sqrt(exp(OM[i, i]) - 1) * 100, digits = 4))
  }

  PrinTxt(nEta + 18, 2, "Standard Error of Omega Matrix", Cex = 0.9)
  sSeOM <- capture.output(SeOM)
  for (i in 1:length(sSeOM)) {
    PrinTxt(nEta + 19 + i, 5, sSeOM[i])
  }

  PrinTxt(2 * nEta + 23, 2, "Relative Standard Error(%) of Omega Matrix", Cex = 0.9)
  sRSEOM <- capture.output(RSEOM)
  for (i in 1:length(sRSEOM)) {
    PrinTxt(2 * nEta + 24 + i, 5, sRSEOM[i])
  }

  PrinTxt(3 * nEta + 27, 1, "Sigmas", Cex = 1)
  PrinTxt(3 * nEta + 29, 3, paste("Number of Epsilons       :", nEps))
  if (nEps == 1 && SG == 1 & SeSG == 1e+10) {
    PrinTxt(3 * nEta + 31, 3, "Fixed as 1")
  } else {
    sSG <- capture.output(SG)
    for (i in 1:length(sSG)) {
      PrinTxt(3 * nEta + 30 + i, 5, sSG[i])
    }
  }

  ClosePDF()
  message("S2-Parameters.PDF generated.")
}
