TrimOut = function(inFile, outFile="PRINT.OUT")
{
  if (missing(inFile)) {
    vPath = strsplit(getwd(), "/")[[1]]
    n = length(vPath)
    inFile = paste0(strsplit(vPath[n],"(\\.)")[[1]][1], ".OUT")
  }
  if (!file.exists(inFile)) stop("File does not exist!")

  O1 = readLines(inFile)
  x1 = which(startsWith(O1, "$PROB"))      # fo find local times
  x2 = which(startsWith(O1, "Stop Time:")) # to find local times  
  xLines = c(x1 - 1, x2 + 2)
  nLine = length(O1)
  MaxW = 1023

  x3 = which(startsWith(O1, "1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) "))
  if (is.na(x3[1])) {
    Line1 = nLine
  } else {
    Line1 = x3[1]
  }

  Res = O1[1]
  cLine = 2
  for (i in 2:nLine) {
    if (i %in% xLines | i < Line1 ) {
      Res[cLine] = O1[i]
      cLine = cLine + 1
    } else {
      Ch1 = substr(O1[i], 1, 1)
      if (Ch1 == "1" | Ch1 == "0") {
        Res[cLine] = ""
        Res[cLine + 1] = substr(O1[i], 2, MaxW)
        cLine = cLine + 2
      } else if (Ch1 == "+") {
        Res[cLine - 1] = paste0(Res[cLine - 1], substr(O1[i], nchar(O1[i - 1]) + 1, MaxW))
      } else if (Ch1 == " ") {
        Res[cLine] = substr(O1[i], 2, MaxW)
        cLine = cLine + 1
      } else {
        Res[cLine] = O1[i]
        cLine = cLine + 1      
      }
    }
  }

  if (outFile == "") { 
    return(Res)
  } else { 
    writeLines(Res, outFile)
  }
}
