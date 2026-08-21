#' Histogram Panel for Pairs Plot
#'
#' @param x numeric vector
#' @param ... additional arguments
#' @keywords internal
PanelHist <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))

  par(usr = c(usr[1:2], 0, 1.5))
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks
  nB <- length(breaks)
  y <- h$counts
  y <- y / max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}


#' Correlation Panel for Pairs Plot
#'
#' @param x numeric vector
#' @param y numeric vector
#' @param digits integer, number of digits for correlation
#' @param prefix character, prefix for text
#' @param cex.cor numeric, character expansion for correlation
#' @keywords internal
PanelCor <- function(x, y, digits = 2, prefix = "", cex.cor) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste(prefix, txt, sep = "")
  cex <- if (missing(cex.cor)) 0.8 / strwidth(txt) else cex.cor
  text(0.5, 0.5, txt, cex = min(1, cex * abs(r) + 0.5))
}


#' Character Panel for Pairs Plot
#'
#' @param x numeric vector
#' @param y numeric vector
#' @keywords internal
PanelChar <- function(x, y) {
  IDs <- get("ParDat", envir = nmw_env)[, "ID"]
  lines(lowess(x, y), col = "red")
  points(x, y, cex = 0.5)
}


#' EBE Pair Plot
#'
#' Creates pairs plot of ETAs vs covariates with histograms and correlations.
#'
#' @param tabEta data.frame with ID, covariates, and ETA columns
#' @param RunNumber character, model/run identifier
#' @param each logical, if TRUE creates separate plots per covariate
#' @export
EBEpair <- function(tabEta, RunNumber, each = FALSE) {
  assign("ParDat", tabEta, envir = nmw_env)
  NameList <- names(tabEta)
  nVar <- length(NameList)
  nEta <- sum(substr(names(tabEta), 1, 3) == "ETA")
  EtaNames <- vector()
  for (i in 1:nEta) EtaNames <- append(EtaNames, paste("ETA", i, sep = ""))

  if (each == FALSE) {
    pairs(tabEta[, -1], main = paste("Covariate vs ETA of", RunNumber), gap = 0,
          lower.panel = PanelChar,
          diag.panel = PanelHist,
          upper.panel = PanelCor)
  } else {
    for (i in 1:nVar) {
      if (substr(NameList[i], 1, 3) != "ETA" & NameList[i] != "ID") {
        pairs(tabEta[, c(NameList[i], EtaNames)],
              main = paste(NameList[i], "vs ETAs of Model", RunNumber), gap = 0,
              lower.panel = PanelChar,
              diag.panel = PanelHist,
              upper.panel = PanelCor)
      }
    }
  }
}

# Package-level environment for pair plot data sharing
nmw_env <- new.env(parent = emptyenv())
