#' Match End of String
#' @param vName0 character vector of names
#' @param End0 character string to match at end
#' @param IgnoreCase logical, whether to ignore case
#' @return logical vector
#' @keywords internal
MatchEnd <- function(vName0, End0, IgnoreCase = TRUE) {
  if (IgnoreCase == TRUE) {
    vName <- toupper(vName0)
    End <- toupper(End0)
  } else {
    vName <- vName0
    End <- End0
  }

  n <- length(vName)
  vRet <- rep(FALSE, n)
  nc1 <- nchar(End)

  for (i in 1:n) {
    nc2 <- nchar(vName[i])
    if (nc2 >= nc1) {
      if (substring(vName[i], nc2 - nc1 + 1, nc2) == End) {
        vRet[i] <- TRUE
      }
    }
  }
  return(vRet)
}
