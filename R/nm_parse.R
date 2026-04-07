#' Parse Value After Tag from NONMEM Output
#'
#' Extracts the value following a tag string (after the colon) in NONMEM output lines.
#'
#' @param TagString character, the tag to search for
#' @param RawRead character vector of lines from NONMEM output
#' @return character, the value after the tag
#' @export
ParseOut <- function(TagString, RawRead) {
  Inter <- grep(TagString, RawRead, value = TRUE)
  if (length(Inter) > 0) {
    return(strsplit(Inter, ":")[[1]][2])
  } else {
    return(0)
  }
}


#' Parse Data Item Names from NONMEM Output
#'
#' @param RawRead character vector of lines from NONMEM output
#' @return character vector of item names
#' @keywords internal
ParseItemName <- function(RawRead) {
  Tag1 <- "0LABELS FOR DATA ITEMS:"
  Tag2 <- "LABELS FOR PRED-DEFINED ITEMS:"

  Line1 <- grep(Tag1, RawRead)
  Line2 <- grep(Tag2, RawRead)

  Ans <- vector()
  for (i in (Line1 + 1):(Line2 - 1)) {
    Ans <- c(Ans, strsplit(RawRead[i], " ")[[1]][-1])
  }
  return(Ans)
}


#' Parse Fixed-Width Format Item Widths from NONMEM Output
#'
#' @param RawRead character vector of lines from NONMEM output
#' @return list with LineCount and Widths
#' @keywords internal
ParseItemWidth <- function(RawRead) {
  Line1 <- grep("0FORMAT FOR DATA:", RawRead)
  FORM <- trimws(RawRead[Line1 + 1])
  nChar <- nchar(FORM)
  FORM <- substr(FORM, 2, nChar - 1)
  f2 <- strsplit(FORM, ",")[[1]]
  nF <- length(f2)
  newF <- vector()
  for (i in 1:nF) {
    Pos <- grep("(\\()", f2[i])
    if (length(Pos) > 0) {
      tStr <- f2[i]
      nChar <- nchar(tStr)
      for (j in 1:nChar) {
        if (substr(tStr, j, j) == "(") Pos1 <- j
        if (substr(tStr, j, j) == ")") Pos2 <- j
      }
      Mul <- as.numeric(substr(tStr, 1, Pos1 - 1))
      Term <- substr(tStr, Pos1 + 1, Pos2 - 1)
      for (j in 1:Mul) {
        newF <- c(newF, Term)
      }
    } else {
      newF <- c(newF, f2[i])
    }
  }

  LineCount <- 1
  for (i in 1:length(newF)) {
    if (substr(newF[i], nchar(newF[i]), nchar(newF[i])) == "/") {
      LineCount <- LineCount + 1
      newF[i] <- substr(newF[i], 1, nchar(newF[i]) - 1)
    }
  }

  Fmts <- chartr("F", "E", newF)
  n <- length(Fmts)
  Widths <- vector()
  for (i in 1:n) {
    x <- strsplit(Fmts[i], "E")[[1]]
    if (x[1] == "") {
      r <- 1
    } else {
      r <- as.integer(x[1])
    }
    for (j in 1:r) {
      Widths <- c(Widths, as.integer(x[2]))
    }
  }
  Result <- list(LineCount, Widths)
  names(Result) <- c("Lines", "Widths")
  return(Result)
}


#' Extract Value Between XML Tags
#'
#' @param Tag1 opening tag
#' @param Tag2 closing tag
#' @param RawRead character vector of lines
#' @return character, the value between tags
#' @keywords internal
BtwTagVal <- function(Tag1, Tag2, RawRead) {
  return(sub(Tag2, "", sub(Tag1, "", grep(paste(Tag1, "(.*?)", Tag2, sep = ""),
                                          RawRead, value = TRUE))))
}


#' Extract Lines Between XML Tags
#'
#' @param Tag1 opening tag
#' @param Tag2 closing tag
#' @param RawRead character vector of lines
#' @return character vector of lines between tags
#' @keywords internal
BtwTagLines <- function(Tag1, Tag2, RawRead) {
  Line1 <- grep(Tag1, RawRead)
  Line2a <- grep(Tag2, RawRead)

  if (length(Line1) > 0) {
    if (length(Line2a) > 1) {
      i <- 1
      while (Line2a[i] < Line1) i <- i + 1
      Line2 <- Line2a[i]
    } else {
      Line2 <- Line2a
    }
    return(RawRead[(Line1 + 1):(Line2 - 1)])
  } else {
    return(0)
  }
}


#' Extract Vector of Values Between XML Tags
#'
#' Parses NONMEM XML output to extract theta, omega, or sigma vectors.
#'
#' @param Tag character, the tag name (e.g., "nm:theta")
#' @param RawRead character vector of XML lines
#' @return numeric vector of values
#' @export
BtwTagVals <- function(Tag, RawRead) {
  Tag1 <- paste("<", Tag, ">", sep = "")
  Tag2 <- paste("</", Tag, ">", sep = "")
  Line1 <- grep(Tag1, RawRead)
  Line2 <- grep(Tag2, RawRead)
  if (length(Line1) != 0 && length(Line2) != 0) {
    TmpRaw <- RawRead[(Line1 + 1):(Line2 - 1)]
    nVal <- length(TmpRaw)
    Vals <- rep(0, nVal)
    for (i in 1:nVal) {
      Vals[i] <- as.numeric(BtwTagVal(paste("<nm:val nm:name='", i, "'>", sep = ""),
                                      "</nm:val>", TmpRaw[i]))
    }
    return(Vals)
  } else {
    return(0)
  }
}


#' Extract Matrix from NONMEM XML Output
#'
#' Parses NONMEM XML output to extract omega or sigma matrices.
#'
#' @param Tag character, the tag name (e.g., "omega", "sigma")
#' @param RawRead character vector of XML lines
#' @param nRow integer, dimension of the matrix
#' @return numeric matrix
#' @export
BtwTagMat <- function(Tag, RawRead, nRow) {
  Tag1 <- paste("<nm:", Tag, ">", sep = "")
  Tag2 <- paste("</nm:", Tag, ">", sep = "")
  MatLines <- BtwTagLines(Tag1, Tag2, RawRead)
  RetMat <- matrix(rep(1e+10, nRow * nRow), nrow = nRow, ncol = nRow)

  if (length(MatLines) > 0 & MatLines[1] != 0) {
    for (i in 1:nRow) {
      MatRow <- BtwTagLines(paste("<nm:row nm:rname='", i, "'>", sep = ""),
                            "</nm:row>", MatLines)
      for (j in 1:i) {
        RetMat[i, j] <- RetMat[j, i] <- as.double(
          BtwTagVal(paste("<nm:col nm:cname='", j, "'>", sep = ""),
                    "</nm:col>", MatRow))
      }
    }
  }
  return(RetMat)
}
