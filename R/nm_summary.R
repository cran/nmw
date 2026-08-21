#' Get Model Names from Folder Names
#'
#' Extracts model names by splitting folder names on period.
#'
#' @param vName character vector of folder names
#' @return character vector of model names
#' @keywords internal
GetModelNames <- function(vName) {
  ModelNames <- character()
  n <- length(vName)
  for (i in 1:n) {
    vStr <- strsplit(vName[i], "(\\.)")[[1]]
    if (length(vStr) != 2) stop("Folder Name should not have more than one period!")
    ModelNames <- c(ModelNames, vStr[1])
  }
  return(ModelNames)
}


#' Get Termination Reason from PRINT.OUT
#'
#' @param folderName character, optional folder to set as working directory
#' @return character, termination reason code
#' @keywords internal
GetReason <- function(folderName) {
  if (!missing(folderName)) {
    oWorkDir <- setwd(folderName)
    on.exit(setwd(oWorkDir))
  }
  if (!file.exists("PRINT.OUT")) return("UNKNOWN")
  OUT <- readLines("PRINT.OUT")
  ln0 <- grep("#TERM", OUT)
  if (length(ln0) < 1) return("")

  res0 <- OUT[ln0[1] + 3]

  if (substr(res0, 1, 34) == "NO. OF FUNCTION EVALUATIONS USED: ") {
    res1 <- "SUCCESS"
  } else if (res0 == "HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.") {
    res1 <- "HOWEVER"
  } else if (substr(res0, 1, 32) == "DUE TO ROUNDING ERRORS (ERROR=13") {
    res1 <- "ROUNDING"
  } else if (res0 == "DUE TO MAX. NO. OF FUNCTION EVALUATIONS EXCEEDED") {
    res1 <- "MAXEVAL"
  } else if (res0 %in% c("DUE TO PROXIMITY OF NEXT ITERATION EST. TO A VALUE",
                          "DUE TO PROXIMITY OF LAST ITERATION EST. TO A VALUE")) {
    res1 <- "INFINITY"
  } else if (substr(res0, 1, 43) == "N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION: ") {
    res0 <- OUT[ln0[1] - 2]
    if (res0 == "HESSIAN OF POSTERIOR DENSITY IS NON-POSITIVE-DEFINITE DURING SEARCH") {
      res1 <- "POSTERIOR"
    } else if (res0 == "AT INITIAL OBJ. FUNCTION EVALUATION") {
      res0 <- OUT[ln0[1] - 7]
      if (substr(res0, 1, 32) == "OCCURS DURING SEARCH FOR ETA AT ") {
        res1 <- "EBE"
      } else {
        res1 <- res0
      }
    } else if (substr(res0, 1, 34) == "Second Error signal from LINV1PD: ") {
      res1 <- "LINV1PD"
    } else {
      res1 <- res0
    }
  } else if (is.numeric(as.numeric(substr(trimws(res0), 1, 10)))) {
    res1 <- "UNKNOWN"
  } else {
    res1 <- res0
  }
  return(res1)
}


#' Summarize Single NONMEM Run
#'
#' Collects key metrics from a single NONMEM run folder.
#'
#' @param ModelName character, model name
#' @return data.frame with one row of summary metrics
#' @keywords internal
SumOut1 <- function(ModelName) {
  OutName <- ModelName
  Parent <- GetProbVal("P:")
  EstMethod <- GetEstMethod()
  Formula <- GetProbVal("F:")
  Minimize <- MinSuccess()
  Reason <- GetReason()
  OFV <- GetOFV()
  SE <- SESuccess()
  Parameter <- GetCountPara()
  Theta <- GetCountAllTheta()
  FixedTheta <- GetCountFixedTheta()
  Eta <- GetCountEta()
  Eps <- GetCountEps()
  OffDiagOmega <- paste(GetOffDiagOmega(), Eta * (Eta - 1) / 2, sep = "/")
  nDV <- GetCountObs()
  AICc <- GetAICc()
  return(data.frame(OutName, Parent, EstMethod, Formula, Minimize, Reason,
                    OFV, SE, Parameter, Theta, FixedTheta, Eta, Eps,
                    OffDiagOmega, nDV, AICc, stringsAsFactors = FALSE))
}


#' Summarize Multiple NONMEM Runs
#'
#' Scans subdirectories matching RunExt and collects summary metrics
#' for model development flow analysis.
#'
#' @param FileExt character, control file extension (default ".CTL")
#' @param RunExt character, run folder extension matching the NONMEM
#'   `nmfe<ver>.bat` output directory.  Default ".R76" for NONMEM 7.6.
#'   Use ".R75" or ".R74" for older versions.
#' @param OutExt character, output file extension (default ".OUT")
#' @return data.frame with one row per run; columns include OutName,
#'   Parent, EstMethod, Formula, Minimize, Reason, OFV, SE, Parameter,
#'   Theta, FixedTheta, Eta, Eps, OffDiagOmega, nDV, AICc, plus tree
#'   layout columns (NextSib, FirstKid, LastKid, LOrder, Depth, Pos).
#'   Returns invisible NULL with a warning if no run folders are found.
#' @export
SumOut <- function(FileExt = ".CTL", RunExt = ".R76", OutExt = ".OUT") {
  WorkDir <- getwd()
  on.exit(setwd(WorkDir))

  Folders <- list.dirs(path = WorkDir, full.names = FALSE, recursive = FALSE)
  # FIX (v0.3.2): MatchEnd() returns a LOGICAL mask; use it to subset Folders
  RunFolders <- Folders[MatchEnd(Folders, RunExt)]
  nRun <- length(RunFolders)
  if (nRun == 0) {
    warning(sprintf("SumOut: no folders ending in '%s' under %s", RunExt, WorkDir))
    return(invisible(NULL))
  }
  ModelNames <- GetModelNames(RunFolders)

  MDL <- data.frame(OutName = character(), Parent = character(), EstMethod = character(),
                    Formula = character(), Minimize = logical(), Reason = character(),
                    OFV = double(), SE = logical(), Parameter = integer(),
                    Theta = integer(), FixedTheta = integer(), Eta = integer(),
                    Eps = integer(), OffDiagOmega = character(), nDV = integer(),
                    AICc = double())

  for (i in 1:nRun) {
    setwd(paste(WorkDir, RunFolders[i], sep = "/"))
    MDL <- rbind(MDL, SumOut1(ModelNames[i]))
  }

  nCtl <- nrow(MDL)

  MDL <- cbind(MDL, data.frame(NextSib = character(length = nCtl),
                               FirstKid = character(length = nCtl),
                               LastKid = character(length = nCtl),
                               LOrder = integer(length = nCtl),
                               Depth = integer(length = nCtl),
                               Pos = double(length = nCtl), stringsAsFactors = FALSE))

  for (i in 1:nCtl) {
    if (is.na(MDL[i, "Parent"])) MDL[i, "Parent"] <- ""
    if (toupper(MDL[i, "Parent"]) == "ROOT") MDL[i, "Parent"] <- ""
    if (length(intersect(MDL[i, "Parent"], MDL[, "OutName"])) == 0) MDL[i, "Parent"] <- ""
  }

  for (i in 1:nCtl) {
    cCtl <- MDL[i, "OutName"]
    Kids <- MDL[MDL[, "Parent"] == cCtl, ]
    if (nrow(Kids) > 0) {
      MDL[i, "FirstKid"] <- Kids[1, "OutName"]
      MDL[i, "LastKid"] <- tail(Kids[, "OutName"], 1)
    }
    cPar <- MDL[i, "Parent"]
    Sibs <- MDL[MDL[, "Parent"] == cPar, ]
    if (nrow(Sibs) > 1) {
      Loc <- which(Sibs[, "OutName"] == MDL[i, "OutName"]) + 1
      if (Loc <= nrow(Sibs)) MDL[i, "NextSib"] <- Sibs[Loc, "OutName"]
    }
  }

  # Robust forest layout: assign Depth, Pos, LOrder by DFS over ALL roots.
  # Handles a forest (multiple roots, e.g. parents living in another folder),
  # disconnected components, and parent-cycles. The previous single-traversal
  # left unvisited nodes at the initial (Depth=0, Pos=0) -> nodes overlapped.
  for (i in 1:nCtl) {                      # self-parent guard (would be unreachable)
    if (!is.na(MDL[i, "Parent"]) && MDL[i, "Parent"] == MDL[i, "OutName"]) {
      MDL[i, "Parent"] <- ""
    }
  }
  Visited <- logical(nCtl)
  cLOrder <- 0L
  cLeaf <- 0                               # global leaf counter -> unique column per leaf
  LayoutNode <- function(idx, depth) {
    if (Visited[idx]) return(invisible(NULL))   # cycle guard
    Visited[idx] <<- TRUE
    MDL[idx, "Depth"] <<- depth
    MDL[idx, "LOrder"] <<- cLOrder
    cLOrder <<- cLOrder + 1L
    KidIdx <- which(MDL[, "Parent"] == MDL[idx, "OutName"] & !Visited)
    if (length(KidIdx) == 0) {             # leaf: take next free column
      MDL[idx, "Pos"] <<- cLeaf
      cLeaf <<- cLeaf + 1
    } else {                               # internal: centre over children
      for (k in KidIdx) LayoutNode(k, depth + 1)
      KidPos <- MDL[KidIdx, "Pos"]
      MDL[idx, "Pos"] <<- (min(KidPos) + max(KidPos)) / 2
    }
  }
  for (r in which(MDL[, "Parent"] == "")) LayoutNode(r, 0)  # every root + its subtree
  for (i in which(!Visited))              LayoutNode(i, 0)  # any node left in a cycle

  setwd(WorkDir)
  write.csv(MDL, paste0(GetCurModelName(), ".csv"), row.names = FALSE, na = "")
  return(MDL)
}


#' Connect Two Points with Lines
#'
#' Draws connecting lines between parent and child nodes in flow diagram.
#'
#' @param pt1 numeric vector of length 2, first point (x, y)
#' @param pt2 numeric vector of length 2, second point (x, y)
#' @keywords internal
ConnPoint <- function(pt1, pt2) {
  x1 <- pt1[1]; y1 <- pt1[2]
  x2 <- pt2[1]; y2 <- pt2[2]
  if (x1 == x2) {
    lines(c(x1, x2), c(y1, y2))
  } else {
    yi <- (y1 + y2) / 2
    lines(c(x1, x1, x2, x2), c(y1, yi, yi, y2))
  }
}


#' Draw Model Development Flow Diagram
#'
#' Creates a graphical model development flow diagram from SumOut results.
#'
#' @param MDL data.frame from SumOut
#' @param MDLName character, model development name for title
#' @param out.start character, filter start
#' @param out.end character, filter end
#' @param Target character, "PDF" or other
#' @export
Outline <- function(MDL, MDLName, out.start = "0", out.end = "zzzzzzzz", Target = "PDF") {
  ttab <- MDL
  ttab <- ttab[ttab[, "OutName"] >= out.start & ttab[, "OutName"] <= out.end, ]
  MinAICc <- min(ttab$AICc, na.rm = TRUE)

  start.x <- min(ttab$Pos)
  start.y <- min(ttab$Depth)
  maxdepth <- max(ttab$Depth)
  maxwidth <- max(ttab$Pos)
  ttab$Pos <- ttab$Pos * 0.9

  if (maxdepth == 0) maxdepth <- 1

  par(oma = c(1, 0.5, 1, 0.5), mar = c(0.5, 0.2, 1.5, 0.2))

  plot(ttab$Pos, ttab$Depth,
       xlim = c(start.x - maxwidth * 0.05, maxwidth * 1.1),
       ylim = c(maxdepth * 1.1, start.y - maxdepth * 0.1),
       type = "n", xlab = "", ylab = "", axes = FALSE,
       main = paste("Model Development Flow of", toupper(MDLName)))
  mtext("*: Minimization Successful,  #: Std Error Available,  Theta - Eta - Off Diag - Total Parameters",
        side = 1, cex = 0.8)

  for (i in (1:length(ttab$OutName))) {
    sym1 <- ifelse(ttab$Minimize[i] == TRUE, "*", " ")
    col1 <- ifelse(ttab$Minimize[i] == TRUE, "red", "black")
    sym1 <- paste(sym1, substr(ttab$Reason[i], 1, 1))

    sym2 <- ifelse(ttab$SE[i] == TRUE, "#", " ")
    col2 <- ifelse(ttab$SE[i] == TRUE, "blue", "black")

    AICc <- ttab$AICc[i]
    sym3 <- ifelse(is.finite(MinAICc) & AICc == MinAICc, "v", " ")
    col3 <- ifelse(is.finite(MinAICc) & AICc == MinAICc, "brown", "black")

    text(ttab$Pos[i], ttab$Depth[i] - 0.02 * maxdepth,
         paste(ttab$OutName[i], sym1), cex = 0.7, col = col1)
    text(ttab$Pos[i], ttab$Depth[i],
         paste(format(ttab$OFV[i], digits = 6), sym2), cex = 0.65, col = col2)
    text(ttab$Pos[i], ttab$Depth[i] + 0.02 * maxdepth,
         paste(ttab$Theta[i], ttab$Eta[i], ttab$OffDiagOmega[i], ttab$Parameter[i], sep = "-"),
         cex = 0.7)
    text(ttab$Pos[i], ttab$Depth[i] + 0.04 * maxdepth,
         paste(format(AICc, digits = 6), sym3), cex = 0.7, col = col3)

    fs <- ttab[i, "Formula"]
    if (fs != "") {
      fs <- strsplit(fs, ",")[[1]]
      for (j in 1:length(fs)) {
        text(ttab$Pos[i], ttab$Depth[i] + ((j - 1) * 0.016 - 0.02) * maxdepth,
             paste("      ", fs[j]), cex = 0.65, adj = 0)
      }
    }

    pid <- ttab$Parent[i]
    if ((!is.na(pid)) & pid != "" & sum(ttab[, "OutName"] == pid, na.rm = TRUE) == 1) {
      ppt <- c(ttab[ttab$OutName == paste(pid), ]$Pos,
               ttab[ttab$OutName == paste(pid), ]$Depth + max(ttab$Depth) * 0.05)
      ConnPoint(c(ttab$Pos[i], ttab$Depth[i] - max(ttab$Depth) * 0.03), ppt)
    }
  }
}
