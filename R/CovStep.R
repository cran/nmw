CovStep = function()
{
  e$STEP <- "COV"
  StartTime = Sys.time()
  Rmat = Hessian(Obj, e$FinalPara)/2  # FinalPara from EstStep()
  Smat = CalcSmat()
  invR = solve(Rmat)
  Cov = invR %*% Smat %*% invR
  SE = sqrt(diag(Cov))
  Correl = cov2cor(Cov)
#  InvCov = solve(Cov)
  InvCov = Rmat %*% solve(Smat) %*% Rmat
  EigenVal = sort(eigen(Correl)$values) # Sorted Eigenvalues

  RunTime = difftime(Sys.time(), StartTime)
  Result = list(RunTime, SE, Cov, Correl, InvCov, EigenVal, Rmat, Smat)
  names(Result) = list("Time", "Standard Error", "Covariance Matrix of Estimates", "Correlation Matrix of Estimates", "Inverse Covariance Matrix of Estimates", "Eigen Values", "R Matrix", "S Matrix")
  return(Result)
}
