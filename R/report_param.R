#' Generate Parameter Estimates Report (S2-Parameters.PDF)
#'
#' Generates a PDF report with theta, omega, and sigma parameter estimates,
#' standard errors, confidence intervals, and significance flags.
#'
#' @param run_dir character, path to the NONMEM run directory (default: current directory)
#' @param file character, output PDF filename
#' @param model character, run/model name; NULL auto-detects via GetCurModelName()
#' @return Invisibly, a list with the output \code{file}, number of \code{pages},
#'   and counts of theta, eta, and epsilon parameters.
#' @export
nmw_report_param <- function(run_dir = getwd(), file = "S2-Parameters.PDF", model = NULL) {
  owd <- setwd(run_dir); on.exit(setwd(owd))

  ## ---- parsing + statistics: reused verbatim from nmw ----------------------
  ## (GetCurModelName, ReadLastTable, CountEXTParams, BtwTagVals, BtwTagMat are
  ##  nmw functions; in the real package these are simply the package's own.)
  CtlName <- if (is.null(model)) GetCurModelName() else model
  XML <- readLines(paste0(CtlName, ".xml"))
  EXT <- ReadLastTable(paste0(CtlName, ".ext"))
  EXT <- EXT[EXT[, "ITERATION"] >= 0, ]
  params <- CountEXTParams(EXT)
  nThetaAll <- params$nTheta; nEtaAll <- params$nEta; nEpsAll <- params$nEps

  THETA   <- as.double(BtwTagVals("nm:theta",   XML))
  THETASE <- as.double(BtwTagVals("nm:thetase", XML))
  OMEGA   <- BtwTagMat("omega",   XML, nEtaAll)
  OMEGAse <- BtwTagMat("omegase", XML, nEtaAll)
  SIGMA   <- BtwTagMat("sigma",   XML, nEpsAll)
  SIGMAse <- BtwTagMat("sigmase", XML, nEpsAll)

  Thetas <- cbind(THETA, THETASE); OMa <- rbind(OMEGA, OMEGAse); SGa <- rbind(SIGMA, SIGMAse)
  nThAll <- length(Thetas[, 1]); Fixed <- vector(); Unfixed <- vector()
  for (i in 1:nThAll) if (Thetas[i, 2] == 1e+10) Fixed <- c(Fixed, i) else Unfixed <- c(Unfixed, i)
  nFixedTh <- length(Fixed); nUnfixedTh <- length(Unfixed)
  ThRowName <- character(); for (i in 1:nThAll) ThRowName <- c(ThRowName, paste("Theta", i))
  rownames(Thetas) <- ThRowName; colnames(Thetas) <- c("Point Estimate", "Standard Error")
  LL <- Thetas[, 1] - 2 * Thetas[, 2]; UL <- Thetas[, 1] + 2 * Thetas[, 2]
  ZERO <- Thetas[, 2] / abs(Thetas[, 1]) > 0.5
  ONE <- (Thetas[, 1] - 2 * Thetas[, 2] - 1) * (Thetas[, 1] + 2 * Thetas[, 2] - 1) < 0 |
    (Thetas[, 1] - 2 * Thetas[, 2] + 1) * (Thetas[, 1] + 2 * Thetas[, 2] + 1) < 0
  Thetas <- cbind(Thetas, LL, UL, ZERO, ONE); UnfixedThetas <- Thetas[Unfixed, , drop = FALSE]
  nEta <- length(OMa[1, ]); OM <- OMa[1:nEta, , drop = FALSE]; SeOM <- OMa[(nEta + 1):(2 * nEta), , drop = FALSE]
  EtaNames <- character(); for (i in 1:nEta) EtaNames <- c(EtaNames, paste("Eta", i))
  rownames(OM) <- colnames(OM) <- rownames(SeOM) <- colnames(SeOM) <- EtaNames
  RSEOM <- SeOM / abs(OM) * 100
  for (i in 1:nEta) for (j in i:nEta) if (j > i) OM[i, j] <- OM[i, j] / sqrt(OM[i, i] * OM[j, j])
  nEps <- length(SGa[1, ]); SG <- SGa[1:nEps, ]; SeSG <- SGa[(nEps + 1):(2 * nEps), ]

  ## ---- PDF via simPDF (replaces 72 lines of dead-reckoned PrinTxt/AddPage) --
  EtaCV <- matrix(round(sqrt(exp(diag(OM)) - 1) * 100, 4), nrow = 1,
                  dimnames = list("CV(%)", EtaNames))

  blocks <- list(
    simPDF::block_para("Summary 2 - Parameters", size = 16, font = 2),
    simPDF::block_rule(),
    simPDF::block_keep(list(
      simPDF::block_para("Thetas", size = 13, font = 2),
      simPDF::block_para(sprintf("Number of All / Fixed / Unfixed Thetas : %d / %d / %d",
                                 nThAll, nFixedTh, nUnfixedTh), size = 10))),
    if (nFixedTh > 0) simPDF::block_keep(list(
      simPDF::block_para("Fixed Theta Values", size = 11, font = 2),
      simPDF::block_pre(sprintf("Theta %d : %s", Fixed, format(Thetas[Fixed, 1])), size = 9))),
    simPDF::block_para("Estimated Thetas", size = 11, font = 2),
    simPDF::block_table(round(UnfixedThetas, 5), size = 9, family = "mono"),
    simPDF::block_para("*LL: Lower Limit   UL: Upper Limit   (Point Estimate +/- 2*SE)", size = 8),
    simPDF::block_para(" ZERO: maybe zero? 0=No 1=Yes    ONE: maybe one? 0=No 1=Yes", size = 8),
    simPDF::block_spacer(6),
    simPDF::block_keep(list(
      simPDF::block_para("Omegas", size = 13, font = 2),
      simPDF::block_para(sprintf("Number of Etas : %d", nEta), size = 10))),
    simPDF::block_para("Omega Matrix (lower = covariance, upper = correlation, diag = variance)",
                       size = 11, font = 2),
    simPDF::block_matrix(round(OM, 5), size = 9),
    simPDF::block_spacer(4),
    simPDF::block_para("Interindividual Variability CV(%) for exp(eta) model", size = 11, font = 2),
    simPDF::block_matrix(EtaCV, size = 9),
    simPDF::block_spacer(6),
    simPDF::block_para("Standard Error of Omega Matrix", size = 11, font = 2),
    simPDF::block_matrix(round(SeOM, 5), size = 9),
    simPDF::block_spacer(4),
    simPDF::block_para("Relative Standard Error(%) of Omega Matrix", size = 11, font = 2),
    simPDF::block_matrix(round(RSEOM, 4), size = 9),
    simPDF::block_spacer(6),
    simPDF::block_keep(list(
      simPDF::block_para("Sigmas", size = 13, font = 2),
      simPDF::block_para(sprintf("Number of Epsilons : %d", nEps), size = 10))))

  if (nEps == 1 && all(SG == 1) && all(SeSG == 1e+10)) {
    blocks <- c(blocks, list(simPDF::block_para("Fixed as 1", size = 10)))
  } else {
    SGm <- if (is.matrix(SG)) SG else matrix(SG, nrow = 1, dimnames = list("Eps", NULL))
    blocks <- c(blocks, list(simPDF::block_matrix(round(SGm, 5), size = 9)))
  }

  doc <- simPDF::sp_new(file, paper = "letter", family = "Courier", size = 10)
  simPDF::frame_set(doc, top = doc$H - 45, bottom = 45, left = 45, right = doc$W - 45)
  simPDF::flow_run(doc, Filter(Negate(is.null), blocks),
                   footer = simPDF::block_para("CONFIDENTIAL", size = 8, align = "center"))
  np <- doc$page_no
  simPDF::sp_close(doc)
  invisible(list(file = file, pages = np,
                 nTheta = nThAll, nEta = nEta, nEps = nEps))
}
