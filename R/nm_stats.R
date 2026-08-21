#' Check if All Values are NA
#'
#' @param Column vector
#' @return logical
#' @keywords internal
AllNA <- function(Column) all(is.na(Column))


#' Check if All Values are the Same
#'
#' @param Column vector
#' @return logical
#' @keywords internal
AllSame <- function(Column) (length(unique(Column)) == 1)


#' Test Functional Dependency
#'
#' Tests whether DepNameList columns are functionally dependent on
#' DetNameList columns in a NONMEM table.
#'
#' @param NMTable data.frame
#' @param DetNameList character vector of determinant column names
#' @param DepNameList character vector of dependent column names
#' @return logical, TRUE if functionally dependent
#' @keywords internal
FuncDep <- function(NMTable, DetNameList, DepNameList) {
  LenDetName <- length(DetNameList)

  UniqDet <- as.data.frame(unique(NMTable[, DetNameList]))
  UniqTot <- as.data.frame(unique(NMTable[, c(DetNameList, DepNameList)]))
  LenUniqDet <- length(UniqDet[, 1])
  LenUniqTot <- length(UniqTot[, 1])

  if (LenUniqDet != LenUniqTot) return(FALSE)

  if (sum(UniqDet == UniqTot[, 1:LenDetName], na.rm = TRUE) == LenDetName * LenUniqDet) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


#' NONMEM Variable Statistics
#'
#' Computes statistics for each variable in a NONMEM dataset including
#' counts of NA, zero, positive values, unique values, and functional
#' dependency checks.
#'
#' @param NMTable data.frame of NONMEM data
#' @return matrix with variable statistics
#' @export
NMVarStat <- function(NMTable) {
  VarNames <- names(NMTable)
  StatNames <- c("nNA", "nZero", "nPos", "nUniq", "Min", "Max",
                 "ALLNA", "ALLSame", "ALLInt", "ALLReal",
                 "DepID", "DepIDTimeMDV", "DepIDTimeEVID")
  nVar <- length(VarNames)
  nStat <- length(StatNames)

  if (length(intersect("sDT", VarNames)) == 1) {
    TLabel <- "sDT"
  } else {
    TLabel <- "TIME"
  }

  VarStat <- matrix(nrow = nVar, ncol = nStat, dimnames = list(VarNames, StatNames))

  for (i in 1:nVar) {
    TotData <- NMTable[, paste(VarNames[i])]
    UniqData <- unique(TotData)
    nUniq <- length(UniqData)

    VarStat[i, "nNA"] <- sum(is.na(TotData))
    VarStat[i, "nZero"] <- sum(TotData == 0, na.rm = TRUE)
    VarStat[i, "nPos"] <- sum(TotData > 0, na.rm = TRUE)
    VarStat[i, "nUniq"] <- nUniq
    VarStat[i, "Min"] <- min(UniqData, na.rm = TRUE)
    VarStat[i, "Max"] <- max(UniqData, na.rm = TRUE)
    VarStat[i, "ALLNA"] <- AllNA(UniqData)
    VarStat[i, "ALLSame"] <- AllSame(UniqData)
    VarStat[i, "ALLInt"] <- is.integer(UniqData)
    VarStat[i, "ALLReal"] <- is.double(UniqData)
    VarStat[i, "DepID"] <- FuncDep(NMTable, "ID", paste(VarNames[i]))
    if (length(intersect("MDV", VarNames)) == 1) {
      VarStat[i, "DepIDTimeMDV"] <- FuncDep(NMTable, c("ID", TLabel, "MDV"), paste(VarNames[i]))
    }
    if (length(intersect("EVID", VarNames)) == 1) {
      VarStat[i, "DepIDTimeEVID"] <- FuncDep(NMTable, c("ID", TLabel, "EVID"), paste(VarNames[i]))
    }
  }

  return(VarStat)
}


#' Integer Variable Statistics
#'
#' Prints unique values for all integer variables.
#'
#' @param NMTable data.frame
#' @param VarStat matrix from NMVarStat
#' @keywords internal
IntStat <- function(NMTable, VarStat) {
  VarNames <- rownames(VarStat[VarStat[, "ALLInt"] == TRUE, ])
  nVar <- length(VarNames)
  for (i in 1:nVar) {
    message(VarNames[i])
    message(paste(sort(unique(NMTable[, paste(VarNames[i])])), collapse = " "))
  }
}


#' NONMEM Individual (ID) Statistics
#'
#' Computes per-subject statistics including record counts,
#' DV counts, dosing information, and sorting checks.
#'
#' @param NMTable data.frame of NONMEM data
#' @return list with IDStat data.frame, sortedID logical, sortedDT logical
#' @export
NMIDStat <- function(NMTable) {
  sortedID <- TRUE
  sortedDT <- TRUE

  VarNames <- names(NMTable)
  if (length(intersect("sDT", VarNames)) == 1) {
    TLabel <- "sDT"
  } else {
    TLabel <- "TIME"
  }

  IDList <- unique(NMTable[, "ID"])
  StatNames <- c("nRec", "nDV", "nMDV", "nAMT", "nEVID1", "nEVID2",
                 "nEVID3", "nEVID4", "FRec", "FDRec", "FRDT", "FDDT")
  nID <- length(IDList)
  nStat <- length(StatNames)

  IDStat <- as.data.frame(matrix(data = 0, nrow = nID, ncol = nStat,
                                 dimnames = list(IDList, StatNames)))

  nTotRec <- length(NMTable[, 1])
  pID <- ""
  pDT <- ""

  for (i in 1:nTotRec) {
    cID <- NMTable[i, "ID"]
    cDT <- NMTable[i, TLabel]

    IDStat[paste(cID), "nRec"] <- IDStat[paste(cID), "nRec"] + 1
    if (length(intersect("MDV", VarNames)) == 1) {
      if (NMTable[i, "MDV"] == 0) IDStat[paste(cID), "nDV"] <- IDStat[paste(cID), "nDV"] + 1
      if (NMTable[i, "MDV"] == 1) {
        IDStat[paste(cID), "nMDV"] <- IDStat[paste(cID), "nMDV"] + 1
        if (length(intersect("EVID", VarNames)) == 1) {
          if (NMTable[i, "EVID"] == 1) {
            IDStat[paste(cID), "nEVID1"] <- IDStat[paste(cID), "nEVID1"] + 1
            if (NMTable[i, "AMT"] > 0) IDStat[paste(cID), "nAMT"] <- IDStat[paste(cID), "nAMT"] + 1
          }
          if (NMTable[i, "EVID"] == 2) {
            IDStat[paste(cID), "nEVID2"] <- IDStat[paste(cID), "nEVID2"] + 1
          }
          if (NMTable[i, "EVID"] == 3) {
            IDStat[paste(cID), "nEVID3"] <- IDStat[paste(cID), "nEVID3"] + 1
          }
          if (NMTable[i, "EVID"] == 4) {
            IDStat[paste(cID), "nEVID4"] <- IDStat[paste(cID), "nEVID4"] + 1
            if (NMTable[i, "AMT"] > 0) IDStat[paste(cID), "nAMT"] <- IDStat[paste(cID), "nAMT"] + 1
          }
        } else {
          if (NMTable[i, "AMT"] > 0) IDStat[paste(cID), "nAMT"] <- IDStat[paste(cID), "nAMT"] + 1
        }
      }
    } else {
      IDStat[paste(cID), "nDV"] <- IDStat[paste(cID), "nDV"] + 1
    }

    if (IDStat[paste(cID), "nRec"] == 1) {
      IDStat[paste(cID), "FRec"] <- i
      IDStat[paste(cID), "FRDT"] <- cDT
    }
    if (length(intersect("AMT", VarNames)) == 1) {
      if (IDStat[paste(cID), "FDRec"] == 0 & NMTable[i, "AMT"] > 0) {
        IDStat[paste(cID), "FDRec"] <- i
        IDStat[paste(cID), "FDDT"] <- cDT
      }
    }

    if (pID > cID) sortedID <- FALSE
    if (pID == cID & pDT > cDT) sortedDT <- FALSE

    pID <- cID
    pDT <- cDT
  }
  Result <- list(IDStat, sortedID, sortedDT)
  names(Result) <- c("Individuals", "ID Sorted", "Time Sorted")
  return(Result)
}
