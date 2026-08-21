#' Diagnostic Plot for Post-Processing
#'
#' Creates spaghetti-style diagnostic plots with individual ID labels.
#'
#' @param x numeric vector, x-axis values
#' @param y numeric vector, y-axis values
#' @param mat data.frame with an ID column
#' @param xlbl character, x-axis label
#' @param ylbl character, y-axis label
#' @param smooth character, "T" for lowess smoothing, "F" for identity line
#' @param xlm numeric vector of length 2, x-axis limits
#' @param ylm numeric vector of length 2, y-axis limits
#' @param Log character, log transformation for axes (e.g., "y")
#' @export
DxPlotPost <- function(x, y, mat, xlbl, ylbl, smooth, xlm = "", ylm = "", Log = "") {
  if (Log == "y") {
    x <- x[y > 0]
    y <- y[y > 0]
  }

  if (ylbl == "DV" & (xlbl == "PRED" | xlbl == "IPRE")) {
    xlm <- ylm <- c(min(x, y), max(x, y))
  } else {
    if (identical(xlm, "")) {
      xlm <- c(min(x, na.rm = TRUE), max(x, na.rm = TRUE))
    }
    if (identical(ylm, "")) {
      ylm <- c(min(y, na.rm = TRUE), max(y, na.rm = TRUE))
    }
  }

  plot(0.001, 0.001, type = "n", bty = "o", xlim = xlm, ylim = ylm,
       xlab = xlbl, ylab = ylbl, log = Log)
  if (smooth == "T") {
    lines(lowess(x, y), lty = 2)
    abline(h = 0, lty = 3)
  } else if (xlbl == "PRED" | xlbl == "IPRE") {
    abline(0, 1, lty = 3)
  }
  SUBID <- unique(mat$ID)
  for (i in SUBID) {
    SelID <- mat$ID == i
    x1 <- x[SelID]
    y1 <- y[SelID]
    xord <- sort.list(x1)
    x1 <- x1[xord]
    y1 <- y1[xord]
    lines(x1, y1)
    for (j in 1:length(x1)) text(x1[j], y1[j], i, cex = 0.75)
  }
}
