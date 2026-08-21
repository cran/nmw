#' Read FDATA.CSV with Auto-Detection of Header Format
#'
#' Reads FDATA.CSV, auto-detecting whether it has a NONMEM-style
#' 2-row header (skip=2) or a standard 1-row header (skip=0).
#'
#' @param file character, path to the CSV file (default: "FDATA.CSV")
#' @return data.frame
#' @keywords internal
ReadFDATA <- function(file = "FDATA.CSV") {
  first_line <- readLines(file, n = 1)
  if (grepl("^[A-Za-z]", first_line)) {
    return(read.csv(file))
  } else {
    return(read.csv(file, skip = 2))
  }
}


#' Add Dose Number, Time of Latest Dose, and Time after Latest Dose
#'
#' Adds DoNo (Dosing Occasion Number), ToLD (Time of Latest Dose), and
#' TaLD (Time after Latest Dose) columns to a NONMEM dataset.
#'
#' @param NMData data.frame of NONMEM dataset
#' @param ID character, column name for subject ID
#' @param TIME character, column name for time
#' @param AMT character, column name for dose amount
#' @param II character, column name for interdose interval
#' @param ADDL character, column name for additional doses
#' @param MDV character, column name for missing dependent variable flag
#' @return data.frame with DoNo, ToLD, TaLD columns added
#' @export
AddDoNoTaLD <- function(NMData, ID = "ID", TIME = "TIME", AMT = "AMT",
                        II = "II", ADDL = "ADDL", MDV = "MDV") {
  nRec <- NROW(NMData)
  DoNo <- vector(length = nRec)
  ToLD <- vector(length = nRec)
  TaLD <- vector(length = nRec)

  LastID <- -1
  for (i in seq_len(nRec)) {
    CurID <- NMData[i, ID]
    if (LastID != CurID) {
      pToLD <- 0
      cDoNo <- 0
    }

    if (length(intersect(names(NMData), ADDL)) == 1) {
      DosingHist <- NMData[NMData[, ID] == CurID & NMData[, AMT] > 0 & NMData[, MDV] == 1,
                           c(TIME, AMT, II, ADDL)]
      nDoseRec <- length(DosingHist[, AMT])
      for (j in seq_len(nDoseRec)) {
        cADDL <- DosingHist[j, ADDL]
        if (cADDL > 0) {
          cAMT <- DosingHist[j, AMT]
          cTIME <- DosingHist[j, TIME]
          cII <- DosingHist[j, II]
          for (k in 1:cADDL) DosingHist <- rbind(DosingHist, c(cTIME + k * cII, cAMT, NA, NA))
        }
      }
      DosingHist2 <- DosingHist[order(DosingHist[, TIME]), c(TIME, AMT)]
    }

    if (NMData[i, AMT] > 0 & NMData[i, MDV] == 1) {
      cToLD <- NMData[i, TIME]
      cDoNo <- cDoNo + 1
    } else {
      if (length(intersect(names(NMData), ADDL)) == 1) {
        cTime <- NMData[i, TIME]
        cToLD <- max(DosingHist2[DosingHist2[, TIME] < cTime, TIME])
        cDoNo <- length(DosingHist2[DosingHist2[, TIME] < cTime, TIME])
      } else {
        cToLD <- pToLD
      }
    }
    DoNo[i] <- cDoNo
    ToLD[i] <- cToLD
    LastID <- CurID
    pToLD <- cToLD
  }

  TaLD <- NMData[, TIME] - ToLD
  return(cbind(NMData, DoNo, ToLD, TaLD))
}


#' Expand Dose History with ADDL/II Records
#'
#' Expands compressed dosing records (ADDL/II) into individual dose records.
#'
#' @param DoseHistTab data.frame with columns TIME, AMT, II, ADDL
#' @return data.frame with expanded dose records
#' @export
ExpandDoseHist <- function(DoseHistTab) {
  RetTab <- DoseHistTab
  nDoseRec <- dim(RetTab)[1]
  TempRec <- matrix(nrow = 1, ncol = dim(RetTab)[2])
  colnames(TempRec) <- colnames(RetTab)

  for (j in seq_len(nDoseRec)) {
    cADDL <- RetTab[j, "ADDL"]
    if (cADDL > 0) {
      cTIME <- RetTab[j, "TIME"]
      cAMT <- RetTab[j, "AMT"]
      cII <- RetTab[j, "II"]
      for (k in 1:cADDL) {
        RetTab[j, "II"] <- 0
        RetTab[j, "ADDL"] <- 0
        TempRec <- RetTab[j, ]
        TempRec[1, "TIME"] <- cTIME + k * cII
        TempRec[1, "AMT"] <- cAMT
        RetTab <- rbind(RetTab, TempRec)
      }
    }
  }

  RetTab <- RetTab[order(RetTab[, "TIME"]), ]

  # Check too short interval
  IIs <- unique(DoseHistTab[DoseHistTab[, "II"] > 0, "II"])
  mII <- min(IIs, na.rm = TRUE)

  nRetRec <- dim(RetTab)[1]
  pTIME <- -1 * max(IIs, na.rm = TRUE)
  pAMT <- 0
  for (j in seq_len(nRetRec)) {
    cTIME <- RetTab[j, "TIME"]
    cAMT <- RetTab[j, "AMT"]
    if ((cTIME - pTIME) < 0.5 * mII & cAMT > 0 & pAMT > 0) {
      message("Warning: Maybe too short interval")
      message(paste(capture.output(RetTab[j, ]), collapse = "\n"))
    }
    pTIME <- cTIME
    pAMT <- cAMT
  }

  return(RetTab)
}


#' Remove Rows with NA Values
#'
#' @param RawData data.frame
#' @return data.frame with rows containing NA removed
#' @export
RemoveNA <- function(RawData) {
  nRow <- nrow(RawData)
  Index <- vector(length = nRow)
  for (i in seq_len(nRow)) {
    Index[i] <- !any(is.na(RawData[i, ]))
  }
  return(RawData[Index, ])
}


#' Remove Columns from Table
#'
#' @param Tab data.frame
#' @param OldCols character vector of column names to remove
#' @return data.frame with specified columns removed
#' @keywords internal
RmvCol <- function(Tab, OldCols) {
  return(Tab[, setdiff(colnames(Tab), OldCols)])
}


#' Rename Columns
#'
#' @param ColList character vector of column names
#' @param OldCol old column name
#' @param NewCol new column name
#' @return character vector with renamed columns
#' @keywords internal
RenCol <- function(ColList, OldCol, NewCol) {
  nCol <- length(ColList)
  for (i in seq_len(nCol)) if (ColList[i] == OldCol) ColList[i] <- NewCol
  return(ColList)
}


#' Remove Fixed (Non-Varying) Columns
#'
#' Removes columns where all values are the same across rows.
#'
#' @param Table data.frame or matrix
#' @return data.frame with fixed columns removed
#' @keywords internal
RmvFixed <- function(Table) {
  ColNames <- colnames(Table)
  NotFixed <- rep(FALSE, length(ColNames))
  for (i in seq_len(dim(Table)[1])) {
    for (j in seq_along(ColNames)) {
      if (Table[i, ColNames[j]] != Table[1, ColNames[j]]) NotFixed[j] <- TRUE
    }
  }
  return(Table[, ColNames[NotFixed]])
}


#' Remove Zero Rows/Cols from Symmetric Matrix
#'
#' @param SymmMat symmetric matrix
#' @return matrix with zero rows/columns removed
#' @keywords internal
RmvZero <- function(SymmMat) {
  nRow <- dim(SymmMat)[1]
  NonZero <- rep(FALSE, nRow)

  for (i in seq_len(nRow)) {
    for (j in 1:i) {
      if (SymmMat[i, j] != 0) NonZero[i] <- TRUE
    }
  }
  return(as.matrix(SymmMat[NonZero, NonZero]))
}
