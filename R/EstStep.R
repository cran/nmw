EstStep = function()
{
  e$STEP = "EST"
  OFVinit = Obj(rep(0.1, e$nPara))
  StartTime = Sys.time()
  Res = optim(rep(0.1, e$nPara), Obj, method="L-BFGS-B")
  RunTime = difftime(Sys.time(), StartTime)
  e$FinalPara = DECN(Res$par)
  Result = list(OFVinit, RunTime, Res, e$FinalPara)
  names(Result) = list("Initial OFV", "Time", "Optim", "Final Estimates")
  return(Result)
}
