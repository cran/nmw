#' Get Current Model Name from Working Directory
#'
#' Detects model name by finding the .xml file in the current directory.
#' Falls back to extracting from folder name if no .xml file is found.
#'
#' @return character, the model name
#' @keywords internal
GetCurModelName <- function() {
  xmlFiles <- list.files(pattern = "\\.xml$", ignore.case = TRUE)
  if (length(xmlFiles) >= 1) {
    return(sub("\\.xml$", "", xmlFiles[1], ignore.case = TRUE))
  }
  vPath <- strsplit(getwd(), "/")[[1]]
  n <- length(vPath)
  return(strsplit(vPath[n], "(\\.)")[[1]][1])
}


#' Get Objective Function Value
#'
#' Retrieves the final OFV from the .ext file in the current directory.
#'
#' @return numeric, the objective function value
#' @export
GetOFV <- function() {
  EXTName <- paste0(GetCurModelName(), ".ext")
  if (length(intersect(toupper(list.files()), toupper(EXTName))) == 1) {
    EXT <- ReadLastTable(EXTName)              # chained-$EST safe (last table)
    if ("OBJ" %in% colnames(EXT)) {
      EXT0 <- EXT[EXT[, "ITERATION"] >= 0, ]
      nRow <- NROW(EXT0)
      if (nRow == 0) {
        return(EXT[1, "OBJ"])
      } else {
        return(tail(EXT0[, "OBJ"], 1))
      }
    }
  }
  return(NA)
}


#' Get Estimation Method
#'
#' Identifies the estimation method from the .ext file header.
#'
#' @return character, abbreviation of the estimation method
#' @export
GetEstMethod <- function() {
  EXTName <- paste0(GetCurModelName(), ".ext")
  if (length(intersect(toupper(list.files()), toupper(EXTName))) == 1) {
    TabLines <- grep("TABLE NO", readLines(EXTName), value = TRUE)   # chained-safe: last $EST method
    EstString <- strsplit(TabLines[length(TabLines)], ":")[[1]][2]
    EstMethod <- c("FO", "FOI", "FOCE", "FOCEI", "L", "LI",
                   "FO", "FOI", "FOCE", "FOCEI", "L", "LI")
    names(EstMethod) <- c(
      " First Order",
      " First Order with Interaction",
      " First Order Conditional Estimation",
      " First Order Conditional Estimation with Interaction",
      " Laplacian Conditional Estimation",
      " Laplacian Conditional Estimation with Interaction",
      " First Order (Evaluation)",
      " First Order with Interaction (Evaluation)",
      " First Order Conditional Estimation (Evaluation)",
      " First Order Conditional Estimation with Interaction (Evaluation)",
      " Laplacian Conditional Estimation (Evaluation)",
      " Laplacian Conditional Estimation with Interaction (Evaluation)")
    if (!is.na(EstString)) return(EstMethod[[EstString]])
  }
  return(NA)
}


#' Get PROBLEM Statement Value
#'
#' @param Tag character, the tag to search for in FCON
#' @return character, the extracted value
#' @keywords internal
GetProbVal <- function(Tag) {
  if (length(intersect(toupper(list.files()), "FCON")) == 1) {
    FCON <- readLines("FCON")
    PROB <- strsplit(trimws(substr(FCON[grep("PROB", FCON)], 9, 80)), " ")[[1]]
    Loc <- grep(Tag, PROB)
    # FIX (v0.3.2): grep returns integer(0) when Tag is absent (e.g. CTL has no
    # P:/F: tag, or NM-TRAN truncated $PROB at 80 chars in FCON dropping the tag).
    # Length-zero comparison `if (Loc > 0)` previously raised
    # "argument is of length zero".
    if (length(Loc) > 0 && Loc[1] > 0) {
      return(substr(PROB[Loc[1]], 3, nchar(PROB[Loc[1]])))
    } else {
      return("")
    }
  } else {
    return(NA)
  }
}


#' Check for Tag in File
#'
#' @param FileName character, filename
#' @param Tag character, text to search for
#' @return logical or NA
#' @keywords internal
FileTag <- function(FileName, Tag) {
  if (length(intersect(toupper(list.files()), toupper(FileName))) == 1) {
    if (!is.na(pmatch(Tag, readLines(FileName)))) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    return(NA)
  }
}


#' Check Minimization Success
#'
#' @return logical, TRUE if minimization was successful
#' @export
MinSuccess <- function() {
  return(FileTag("PRINT.OUT", "MINIMIZATION SUCCESSFUL"))
}


#' Check Standard Error Success
#'
#' @return logical, TRUE if standard errors were computed successfully
#' @export
SESuccess <- function() {
  return(FileTag("PRINT.OUT", "********************                            STANDARD ERROR OF ESTIMATE"))
}


#' Get Count of Parameters
#'
#' @return integer, total number of estimated parameters, or NA if file not found
#' @export
GetCountPara <- function() {
  GRDName <- paste0(GetCurModelName(), ".grd")
  if (!file.exists(GRDName)) return(NA)
  GRD <- ReadLastTable(GRDName)              # chained-$EST safe (last table)
  ncol(GRD) - 1
}


#' Get Count of All Thetas
#'
#' @return integer, total number of thetas, or NA if file not found
#' @keywords internal
GetCountAllTheta <- function() {
  EXT <- ReadEXTFile()
  if (is.null(EXT)) return(NA)
  length(grep("THETA", colnames(EXT)))
}


#' Get Count of Unfixed Thetas
#'
#' @return integer
#' @keywords internal
GetCountUnfixedTheta <- function() {
  return(GetCountPara() - GetCountOmega() - GetCountUnfixedEps())
}


#' Get Count of Fixed Thetas
#'
#' @return integer
#' @keywords internal
GetCountFixedTheta <- function() {
  return(GetCountAllTheta() - GetCountUnfixedTheta())
}


#' Get Count of Etas
#'
#' @return integer, number of eta parameters, or NA if file not found
#' @export
GetCountEta <- function() {
  PHIName <- paste0(GetCurModelName(), ".phi")
  if (!file.exists(PHIName)) return(NA)
  PHI <- ReadLastTable(PHIName)              # chained-$EST safe (last table)
  length(grep("ETA", colnames(PHI)))
}


#' Get Count of Unfixed Omega Elements
#'
#' @return integer, or NA if file not found
#' @keywords internal
GetCountOmega <- function() {
  EXT <- ReadEXTFile()
  if (is.null(EXT)) return(NA)
  tOM <- EXT[EXT$ITERATION == -1000000000, grep("OMEGA", colnames(EXT))]
  sum(tOM != 0)
}


#' Get Off-Diagonal Omega Count
#'
#' @return integer
#' @keywords internal
GetOffDiagOmega <- function() {
  Parameter <- GetCountPara()
  Theta <- GetCountAllTheta()
  FixedTheta <- GetCountFixedTheta()
  nEta <- GetCountEta()
  return(Parameter - (Theta - FixedTheta) - nEta)
}


#' Get Count of Epsilons
#'
#' @return integer, number of epsilon parameters, or NA if file not found
#' @export
GetCountEps <- function() {
  EXT <- ReadEXTFile()
  if (is.null(EXT)) return(NA)
  length(grep("SIGMA", colnames(EXT)))
}


#' Get Count of Unfixed Epsilons
#'
#' @return integer, or NA if file not found
#' @keywords internal
GetCountUnfixedEps <- function() {
  EXT <- ReadEXTFile()
  if (is.null(EXT)) return(NA)
  if (!("ITERATION" %in% colnames(EXT))) return(NA)
  EXT <- EXT[EXT[, "ITERATION"] >= 0, ]
  nRow <- NROW(EXT)
  if (nRow == 0) return(0)
  ColNo <- grep("SIGMA", colnames(EXT))
  nCol <- length(ColNo)

  nUnfixed <- 0
  for (i in 1:nCol) {
    for (j in 2:nRow) {
      if (EXT[j, ColNo[i]] != EXT[j - 1, ColNo[i]]) {
        nUnfixed <- nUnfixed + 1
        break
      }
    }
  }
  nUnfixed
}


#' Get Count of Observations
#'
#' @return integer, total number of observation records
#' @export
GetCountObs <- function() {
  if (file.exists("PRINT.OUT")) {
    PRINTOUT <- readLines("PRINT.OUT")
    Tag <- "TOT. NO. OF OBS RECS:"
    LineNo <- grep(Tag, PRINTOUT)
    if (length(LineNo) > 0) {
      return(as.integer(trimws(substr(PRINTOUT[LineNo[1]], nchar(Tag) + 1, 80))))
    } else {
      return(NA)
    }
  } else {
    return(NA)
  }
}


#' Get Corrected AIC (AICc)
#'
#' @return numeric, the corrected AIC value
#' @export
GetAICc <- function() {
  n <- GetCountObs()
  p <- GetCountPara()
  ofv <- GetOFV()
  return(ofv + 2 * p + 2 * p * (p + 1) / (n - p - 1))
}
