# Objective Functions
ObjEta = function(ETAi)
{
# External Variable : e$INTER, e$DATAi, e$THETA, e$invOM, e$SG, e$nEta, e$HNames
# External Function : e$PRED
  FGHDi = e$PRED(e$THETA, ETAi, e$DATAi)
  Ri    = e$DATAi[,"DV"] - FGHDi[,"F"]
  
  if (e$INTER == TRUE) {
    Hi = FGHDi[, e$HNames, drop=FALSE]
  } else {
    FGHD0 = e$PRED(e$THETA, rep(0, e$nEta), e$DATAi)
    Hi = FGHD0[, e$HNames, drop=FALSE]
  }
  
## Slower version
#  Vi = diag(diag(Hi %*% SG %*% t(Hi)))
#  iSum = log(det(Vi)) + t(Yi - Fi)%*%solve(Vi)%*%(Yi - Fi) + t(ETAi) %*% invOM %*% ETAi

## Faster version
  Vi    = diag(Hi %*% e$SG %*% t(Hi))
  return(sum(log(Vi) + Ri*Ri/Vi) + t(ETAi) %*% e$invOM %*% ETAi)
}


Oi = function(i)  # Combined version. You need not use this.
{
#Oi : Calculate individual Objective Function Value with Current Externals
# External Variables: Oi, THETA, ETAi, OM, invOM, SG, DataRef, GNmaes, HNames, DNames, Term123
# External Functions: PRED, ltv2mat
# METHOD = "ZERO" | "COND" | LAPL"
# INTER = TRUE | FALSE  : Here I do not consider. It assumes always with INTER option

# Oi requires basically
#  PRED, THETA, ETAi, DataRef, GNames, HNames, SG, nReci

# Oi requires additionally with METHOD="ZERO"
#  OM

# Oi requires  additionally with METHOD=="COND"
#  invOM, Term123

# Oi requires  additionally with METHOD=="LAPL"
#  invOM, Term123, DNames

  FGHDi = e$PRED(e$THETA, e$ETAi, e$DataRef[[i]]) # ETAi  = rep(0, nEta) # Fixed externally if METHOD="ZERO"
  Yi    = e$DataRef[[i]][, "DV"]
  Fi    = FGHDi[, "F"]
  Gi    = FGHDi[, e$GNames, drop=FALSE]
  Hi    = FGHDi[, e$HNames, drop=FALSE]
  Vi    = diag(Hi %*% e$SG %*% t(Hi))

  if (e$METHOD == "ZERO") {
    Ri     = Yi - Fi
    Ci     = Gi %*% e$OM %*% t(Gi) + diag(Vi, nrow=e$nReci)
    Result = determinant(Ci, logarithm=TRUE)$modulus[[1]] + t(Ri) %*% solve(Ci) %*% Ri
  } else if (e$METHOD == "COND") {
    Hsum = e$invOM + t(Gi) %*% solve(diag(Vi)) %*% Gi
    Result = e$Term123 + determinant(Hsum, logarithm=TRUE)$modulus[[1]]
  } else if (e$METHOD == "LAPL") {
    Hsum = e$invOM
    for (j in 1:e$nReci) {
      Dij = ltv2mat(FGHDi[j, e$DNames])
      Hsum = Hsum + (Gi[j,] %*% t(Gi[j,]) - (Yi[j] - Fi[j])*Dij)/Vi[j]
    }
    Result = e$Term123 + determinant(Hsum, logarithm=TRUE)$modulus[[1]]
  } else {
    print("Undefined Estimation Method")
    Result = NULL
  }
  return(Result)
}


Oi0 = function(i) # Calculate individual OFV with METHOD=="ZERO"
{
# External Variable: THETA, ETAi, DataRef, Gnames, HNames, OM, SG, nReci
# External Function: PRED
  FGHDi = e$PRED(e$THETA, e$ETAi, e$DataRef[[i]]) # ETAi = rep(0, nEta): fixed externally if METHOD=="ZERO"
  Ri    = e$DataRef[[i]][,"DV"] - FGHDi[,"F"]
  Gi    = FGHDi[, e$GNames, drop=FALSE]
  Hi    = FGHDi[, e$HNames, drop=FALSE]
  Ci    = Gi %*% e$OM %*% t(Gi) + diag(diag(Hi %*% e$SG %*% t(Hi)), nrow=e$nReci)
  return(determinant(Ci, logarithm=TRUE)$modulus[[1]] + t(Ri) %*% solve(Ci) %*% Ri)
}


Oi1 = function(i) # Calculate individual OFV with METHOD=="COND" or "LAPL"
{
# External Variable: THETA, ETAi, DataRef, Gnames, HNames, SG, invOM, Oi, Term123, METHOD, nReci
# External Function: PRED
  FGHDi = e$PRED(e$THETA, e$ETAi, e$DataRef[[i]]) # ETAi : predetermined by ObjEta if METHOD=="COND" or "LAPL"
  Yi    = e$DataRef[[i]][,"DV"]
  Fi    = FGHDi[,"F"]
  Gi    = FGHDi[, e$GNames, drop=FALSE]
  Hi    = FGHDi[, e$HNames, drop=FALSE]
  Vi    = diag(Hi %*% e$SG %*% t(Hi))
  if (e$METHOD == "COND") {
    Hsum = e$invOM + t(Gi) %*% solve(diag(Vi)) %*% Gi
  } else if (e$METHOD == "LAPL") {
    Hsum  = e$invOM
    for (j in 1:e$Oi[i,"nRec"]) {
      Dij = ltv2mat(FGHDi[j, e$DNames])
      Hsum = Hsum + (Gi[j,] %*% t(Gi[j,]) - (Yi[j] - Fi[j])*Dij)/Vi[j]
    }
  }
  return(e$Term123 + determinant(Hsum, logarithm=TRUE)$modulus[[1]])
}


Obj = function(vPara) # Calculate total OFV with vectorized parameter (THETA, OM, SG)
{
# External Variable: HNames for ObjEta, ETAi when METHOD=="ZERO", STEP=="EST" or "COV", Oi, EBE, OMindex, SGindex
# External Function: DECN

# Oi requires basically
#  PRED, THETA, ETAi, DataRef, GNames, HNames

# Oi requires  additionally with METHOD=="COND"
#  invOM, Term123

# Oi requires  additionally with METHOD=="LAPL"
#  invOM, Term123, DNames
  if (e$STEP=="EST") {  # "EST" STEP requires DECN, vPara is UCP
    vPara = DECN(vPara)
  }

  e$THETA = vPara[1:e$nTheta]                                  # PRED requires in ObjEta and Oi
  e$OM    = ltv2mat(vPara[e$OMindex])
  e$SG    = diag(vPara[e$SGindex], nrow=e$nEps) # ObjEta rquires

  if (e$METHOD == "ZERO") {  # If METHOD=="ZERO", ETAi should be fixed as rep(0, nEta)
    e$ETAi = rep(0, e$nEta)
    for (i in 1:e$nID)  {
      e$nReci <- e$Oi[i, "nRec"]
      e$Oi[i, "OFVi"] = Oi0(i)  # Oi requires additionally with METHOD="ZERO":  OM, SG
    }
  } else {  # METHOD=="COND" or "LAPL"
#    tryCatch(code = { e$invOM = solve(e$OM)}, error=function(e) { e$invOM = MASS::ginv(e$OM) })
#    tryCatch(code = { e$invOM = solve(e$OM)}, error=function(e) { e$invOM = diag(rep(1, as.integer(e$nEta))) })
    try(e$invOM <- solve(e$OM), silent=TRUE)                                    
    Term2 = determinant(e$OM, logarithm=TRUE)$modulus[[1]]     # If METHOD=="COND" or "LAPL", Oi requires
    for (i in 1:e$nID)  {
      e$DATAi = e$DataRef[[i]]
      e$nReci =- e$Oi[i, "nRec"]
      Res     = optim(e$EBE[i, e$EtaNames], ObjEta, method="BFGS") # ObjEta requires DATAi, THETA, invOM, SG, HNames # Regardless of "COV' or "EST" Step, use ObjEta
      e$ETAi  = as.matrix(Res$par, nrow=e$nEta)
      e$EBE[i, e$EtaNames] = Res$par
      e$Term123 = Term2 + Res$value                            # Oi requires Term123 for METHOD=="COND" or "LAPL"
      e$Oi[i, "OFVi"] = Oi1(i)
    }
  }
  return(sum(e$Oi[,"OFVi"]))
}


OiS0 = function(vPara)
{
# External Variable: DATAi, nTheta, nEta, nEps, nPara, GNames, HNames, ETAi
# External Function: PRED, ltv2mat
  THETA = vPara[1:e$nTheta]
  OM    = ltv2mat(vPara[e$OMindex])
  SG    = diag(vPara[e$SGindex], nrow=e$nEps)
  FGHDi = e$PRED(THETA, e$ETAi, e$DATAi)
  Ri    = e$DATAi[, "DV"] - FGHDi[, "F"]
  Gi    = FGHDi[, e$GNames, drop=FALSE]
  Hi    = FGHDi[, e$HNames, drop=FALSE]
  Ci    = Gi %*% OM %*% t(Gi) + diag(diag(Hi %*% SG %*% t(Hi)))
  return(determinant(Ci, logarithm=TRUE)$modulus[[1]] + t(Ri) %*% solve(Ci) %*% Ri)
}


OiS1 = function(vPara)
{
# External Variable: DATAi, nReci, nTheta, nEta, nEps, nPara, GNames, HNames, DNames
# External Function: PRED, ltv2mat
  THETA = vPara[1:e$nTheta]
  OM    = ltv2mat(vPara[e$OMindex])
  SG    = diag(vPara[e$SGindex], nrow=e$nEps)
  invOM = solve(OM)
  FGHDi = e$PRED(THETA, e$ETAi, e$DATAi)
  Yi    = e$DATAi[,"DV"]
  Fi    = FGHDi[,"F"]
  Gi    = FGHDi[, e$GNames, drop=FALSE]
  Hi    = FGHDi[, e$HNames, drop=FALSE]
  Vi    = diag(Hi %*% SG %*% t(Hi))
  Ri    = Yi - Fi
  Term1 = sum(log(Vi) + Ri*Ri/Vi) # = log(det(Vi)) + t(Yi - Fi)%*%solve(Vi)%*%(Yi - Fi)
  Term2 = determinant(OM, logarithm=TRUE)$modulus[[1]]
  Term3 = t(e$ETAi) %*% invOM %*% e$ETAi
  if (e$METHOD == "COND") {
    Hsum = e$invOM + t(Gi) %*% solve(diag(Vi)) %*% Gi
  } else if (e$METHOD == "LAPL") {
    Hsum  = e$invOM
    for (j in 1:e$nReci) {
      Dij = ltv2mat(FGHDi[j, e$DNames])
      Hsum = Hsum + (Gi[j,] %*% t(Gi[j,]) - (Yi[j] - Fi[j])*Dij)/Vi[j]
    }
  }
  Term4 = determinant(Hsum, logarithm=TRUE)$modulus[[1]]
  return(Term1 + Term2 + Term3 + Term4)
}


CalcSmat = function() # Calculate Smat with METHOD=="ZERO"
{
# External Variables: METHOD, nPara, nTheta, nEta, nEps, nPara, GNames, HNames, DNames
# External Functions: PRED
#  require(numDeriv)
  Smat = matrix(rep(0, e$nPara*e$nPara), nrow=e$nPara, ncol=e$nPara)
  for (i in 1:e$nID) {
    e$DATAi = e$DataRef[[i]]
    e$ETAi  = e$EBE[i, e$EtaNames]
    e$nReci = e$Oi[i,"nRec"]
    if (e$METHOD=="ZERO") {
      gr = grad(OiS0, e$FinalPara)
    } else {
      gr = grad(OiS1, e$FinalPara)
    }
    Smat = Smat + gr %*% t(gr)
  }
  return(Smat/4)
}


PostHocEta = function()
{
  e$THETA = e$FinalPara[1:e$nTheta]
  e$invOM = solve(ltv2mat(e$FinalPara[e$OMindex]))
  e$SG    = diag(e$FinalPara[e$SGindex], nrow=e$nEps)
  for (i in 1:e$nID) {
    e$DATAi = e$DataRef[[i]]
    e$EBE[i, e$EtaNames] = optim(rep(0, e$nEta), ObjEta, method="BFGS")$par
  }
  return(e$EBE)
}

