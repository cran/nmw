TabStep = function()
{
# Uses results in environment
  THETA = e$FinalPara[1:e$nTheta]
  OM = ltv2mat(e$FinalPara[(e$nTheta+1):(e$nTheta + e$nEta*(e$nEta + 1)/2)])
  SG = diag(e$FinalPara[(e$nTheta + e$nEta*(e$nEta + 1)/2 + 1):(e$nTheta + e$nEta*(e$nEta + 1)/2 + e$nEps)])

  tSD = vector()
  for (i in 1:e$nID) {
    cID  = e$IDs[i]
    DATAi = e$DataRef[[i]]
    FGH0 = e$PRED(THETA, rep(0, e$nEta), DATAi )
    G0i  = FGH0[, e$GNames, drop=FALSE]
    H0i  = FGH0[, e$HNames, drop=FALSE]
    R0i  = DATAi[, "DV"] - FGH0[,"F"]
    C0i  = G0i %*% OM %*% t(G0i) + diag(diag(H0i %*% SG %*% t(H0i)))
    WRES = SqrtInvCov(C0i) %*% R0i

    if (e$METHOD != "ZERO") {
      cEBE = e$EBE[e$EBE[,"ID"]==cID, 2:(e$nEta + 1)]
      FGH1 = e$PRED(THETA, cEBE, DATAi)
      G1i  = FGH1[, e$GNames, drop=FALSE]
      H1i  = FGH1[, e$HNames, drop=FALSE]
      R1i  = DATAi[, "DV"] - FGH1[,"F"]
      if (e$INTER == FALSE) {
        C1i  = G1i %*% OM %*% t(G1i) + diag(diag(H0i %*% SG %*% t(H0i)))        
      } else {
        C1i  = G1i %*% OM %*% t(G1i) + diag(diag(H1i %*% SG %*% t(H1i)))
      }
      CWRES = SqrtInvCov(C1i) %*% (R1i + G1i %*% cEBE)
    }

    if (e$METHOD == "ZERO") {
      tSD = rbind(tSD, cbind(DATAi[, c("ID", "TIME", "DV")], PRED=FGH0[,"F"], RES=R0i, WRES, G0i, H0i))
    } else {
      tSD = rbind(tSD, cbind(DATAi[, c("ID", "TIME", "DV")], PRED=FGH0[,"F"], RES=R0i, WRES, CIPREDI=FGH1[,"F"], CIRESI=R1i, CWRES, G1i, H1i))
    }
  }

#  Result = list(THETA, OM, SG, tSD)
#  names(Result) = c("THETA", "OMEGA", "SIGMA", "sdtab")
  return(tSD)
}
