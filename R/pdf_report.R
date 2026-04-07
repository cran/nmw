#' Prepare PDF Output
#'
#' Opens a PDF device for NONMEM report generation.
#'
#' @param FileName character, output PDF filename
#' @param Paper character, paper size (default "letter")
#' @param FontFamily character, font family (default "Courier")
#' @export
PrepPDF <- function(FileName, Paper = "letter", FontFamily = "Courier") {
  pdf(FileName, paper = Paper, width = 8.5, height = 11, family = FontFamily,
      title = "Report of NONMEM Single-run")
}


#' Print Text at Position
#'
#' Prints text at a specific row and column position on the PDF page.
#'
#' @param Row numeric, row position
#' @param Col numeric, column position
#' @param Text character, text to print
#' @param Cex numeric, character expansion factor
#' @export
PrinTxt <- function(Row, Col, Text, Cex = 0.8) {
  text(Col - 1, 2 * Row - 1, Text, cex = Cex, offset = 0)
}


#' Add a New Page to PDF
#'
#' Creates a new page in the PDF output with optional headers and footers.
#'
#' @param Cex numeric, character expansion factor (0.8 or 0.6)
#' @param Header1 character, left header
#' @param Header2 character, center header
#' @param Header3 character, right header
#' @param Footer1 character, left footer
#' @param Footer2 character, center footer
#' @param Footer3 character, right footer
#' @param PrintRowNum logical, whether to print row numbers
#' @param StartRowNum integer, starting row number
#' @export
AddPage <- function(Cex = 0.8, Header1 = "", Header2 = "", Header3 = "",
                    Footer1 = "", Footer2 = "", Footer3 = "",
                    PrintRowNum = FALSE, StartRowNum = 1) {
  if (Cex == 0.8) {
    nRow <- 55
    nCol <- 90
    options(width = 100)
  }
  if (Cex == 0.6) {
    nRow <- 80
    nCol <- 128
    options(width = 150)
  }

  par(oma = c(0, 0, 0, 0), mfrow = c(1, 1), mar = c(0, 0, 0, 0), adj = 0, cex = Cex)
  plot(0, 0, type = "n", ylim = c(2 * nRow - 1, 0), xlim = c(0, nCol - 1),
       xaxt = "n", yaxt = "n", ylab = "", xlab = "", bty = "n")

  if (PrintRowNum == TRUE) {
    for (j in 1:nRow) {
      text(-2 - 1 + 3 - floor(log10(j + StartRowNum - 1)), 2 * j - 1,
           paste(j + StartRowNum - 1, ":", sep = ""), offset = 0)
    }
  }

  text(0, -3, Header1, pos = 4)
  text(nCol / 2, -3, Header2, pos = 3)
  text(nCol, -3, Header3, pos = 2)

  text(0, 2 * nRow + 1, Footer1, pos = 4)
  text(nCol / 2, 2 * nRow + 1, Footer2, pos = 1)
  text(nCol, 2 * nRow + 1, Footer3, pos = 2)
}


#' Print Multiple Lines of Text
#'
#' Prints a character vector of text lines across multiple pages as needed.
#'
#' @param MTxt character vector of text lines
#' @param Cex numeric, character expansion factor
#' @param Header1 character, left header
#' @param Header2 character, center header
#' @param Header3 character, right header
#' @param Footer1 character, left footer
#' @param Footer2 character, center footer
#' @param Footer3 character, right footer
#' @param PrintRowNum logical, whether to print row numbers
#' @export
PrinMTxt <- function(MTxt, Cex = 0.8, Header1 = "", Header2 = "", Header3 = "",
                     Footer1 = "", Footer2 = "", Footer3 = "",
                     PrintRowNum = FALSE) {
  if (Cex == 0.8) {
    nRow <- 55
  }
  if (Cex == 0.6) {
    nRow <- 80
  }
  for (i in 1:length(MTxt)) {
    if (i %% nRow == 1) {
      AddPage(Cex = Cex, Header1 = Header1, Header2 = Header2, Header3 = Header3,
              Footer1 = Footer1, Footer2 = Footer2, Footer3 = Footer3,
              PrintRowNum = PrintRowNum)
    }
    PrinTxt((i - 1) %% nRow + 1, 1, MTxt[i])
  }
}


#' Close PDF Output
#'
#' Closes the current PDF device.
#'
#' @export
ClosePDF <- function() dev.off()


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
