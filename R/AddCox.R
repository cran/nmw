AddCox = function(nmData, coxData, coxCol, dateCol = "DATE", idCol = "ID")
{
  coxData = coxData[order(coxData[, idCol], coxData[, dateCol]), ]
  nmData[, coxCol] = NA
  nrNm = NROW(nmData)
  for (i in 1:nrNm) {
    cID = nmData[i, idCol]
    cDATE = nmData[i, dateCol] 
    cDAT = coxData[coxData[, idCol] == cID & coxData[, dateCol] <= cDATE, , drop=F]
    nrcDAT = NROW(cDAT)
    if (nrcDAT > 0) {
      nmData[i, coxCol] = cDAT[nrcDAT, coxCol]
    } else {
      cDAT = coxData[coxData[, idCol] == cID, , drop=F]
      if (NROW(cDAT) == 0) { nmData[i, coxCol] = NA
      } else { nmData[i, coxCol] = cDAT[1, coxCol] }
    }
  }
  
  if (sum(is.na(nmData[, coxCol])) > 0) warning("There is NA value in the result!")
  return(nmData)
}
