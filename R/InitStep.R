InitStep = function(DataAll, THETAinit, OMinit, SGinit, nTheta, LB=rep(0, nTheta), UB=rep(0, nTheta), Pred, METHOD=METHOD)
{
  e$METHOD = METHOD
  if (METHOD == "ZERO") {
    e$INTER = FALSE
  } else {
    e$INTER = TRUE
  }
  e$PRED = Pred
  e$IDs = unique(DataAll[,"ID"])
  e$nID = length(e$IDs)
  e$Oi = matrix(nrow=e$nID, ncol=3)
  colnames(e$Oi) = c("ID", "OFVi", "nRec")

  e$DataRef = list()
  for (i in 1:e$nID) {
    DATAi = DataAll[DataAll[,"ID"]==e$IDs[i],]
    e$Oi[i, "ID"] = e$IDs[i]
    e$Oi[i, "nRec"] = dim(DATAi)[1]
    e$DataRef[[i]] = DATAi
  }

  e$nTheta = nTheta
  e$nEta   = dim(OMinit)[1]
  e$nEps   = dim(SGinit)[1]
  e$nPara  = e$nTheta + e$nEta*(e$nEta + 1)/2 + e$nEps  # Assume Full Block Omega

  e$GNames = vector()
  for (i in 1:e$nEta) e$GNames[i] = paste0("G",i)

  e$HNames = vector()
  for (i in 1:e$nEps) e$HNames[i] = paste0("H",i)

  e$DNames = vector()
  for (i in 1:e$nEta) for (j in 1:i) e$DNames = append(e$DNames, paste0("D",i,j))

  e$EtaNames = vector(length=e$nEta)
  for (i in 1:e$nEta) e$EtaNames[i] = paste0("ETA",i)

  e$EBE = cbind(e$IDs, matrix(rep(0, e$nID*e$nEta), nrow=e$nID, ncol=e$nEta))
  colnames(e$EBE) = c("ID", e$EtaNames)

  IE = THETAinit
  e$alpha  = 0.1 - log((IE - LB)/(UB - LB)/(1 - (IE - LB)/(UB - LB)))
  e$OMscl  = ScaleVar(OMinit, e$nEta)
  e$SGscl  = ScaleVar(SGinit, e$nEps)

  e$METHOD = METHOD
  e$THETAinit = THETAinit
  e$OMinit = OMinit
  e$SGinit = SGinit
  e$LB     = LB
  e$UB     = UB
  e$OMindex = (e$nTheta + 1):(e$nTheta + e$nEta*(e$nEta + 1)/2)
  e$SGindex = (e$nPara - e$nEps + 1):e$nPara
}
