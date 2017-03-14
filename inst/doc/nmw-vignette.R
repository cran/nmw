library(compiler)
enableJIT(3)
require(nmw)

DataAll = Theoph
colnames(DataAll) = c("ID", "BWT", "DOSE", "TIME", "DV")
DataAll[,"ID"] = as.numeric(as.character(DataAll[,"ID"]))

require(lattice)
xyplot(DV ~ TIME | as.factor(ID), data=DataAll, type="b")

nTheta = 3
nEta = 3
nEps = 2

THETAinit = c(2, 50, 0.1) # Initial estimate
# Omega matrix should be full block.
OMinit = matrix(c(0.2, 0.1, 0.1, 0.1, 0.2, 0.1, 0.1, 0.1, 0.2), nrow=nEta, ncol=nEta)
OMinit
# Sigma matrix should be diagonal.
SGinit = matrix(c(0.1, 0, 0, 0.1), nrow=nEps, ncol=nEps)
SGinit

LB = rep(0, nTheta) # Lower bound
UB = rep(1000000, nTheta) # Upper bound

PRED = function(THETA, ETA, DATAi) # Prediction function
{
  DOSE = 320
  TIME = DATAi[,"TIME"]

  KA = THETA[1]*exp(ETA[1])
  V  = THETA[2]*exp(ETA[2])
  K  = THETA[3]*exp(ETA[3])

  TERM1 = DOSE/V * KA/(KA - K)
  TERM2 = exp(-K*TIME)
  TERM3 = exp(-KA*TIME)

  F  = TERM1 * (TERM2 - TERM3)
  G1 = -F*K/(KA - K) + KA*TIME*TERM1*TERM3
  G2 = -F
  G3 = (F/(KA - K) - TIME*TERM1*TERM2) * K
  H1 = F
  H2 = 1

  if (METHOD=="LAPL") {
    D11 = DOSE*(KA*V**-1.0*(-1.0*KA*(-2.0*KA*(-1.0*K+KA)**-3.0*(-1.0*TERM3+TERM2)+ 
          KA*TIME*TERM3*(-1.0*K+KA)**-2.0)+
          -1.0*KA*(-1.0*K+KA)**-2.0*(-1.0*TERM3+TERM2)+
          KA*TIME*(-1.0*KA*TIME*TERM3*(-1.0*K+KA)**-1.0+
          -1.0*KA*TERM3*(-1.0*K+KA)**-2.0)+
          KA*TIME*TERM3*(-1.0*K+KA)**-1.0)+
          KA*V**-1.0*(-1.0*K+KA)**-1.0*(-1.0*TERM3+TERM2)+
          2.0*KA*V**-1.0*(-1.0*KA*(-1.0*K+KA)**-2.0*(-1.0*TERM3+TERM2)+
          KA*TIME*TERM3*(-1.0*K+KA)**-1.0))
    D21 = -G1
    D22 = F
    D31 = DOSE*(KA*V**-1.0*(KA*K*TIME*TERM2*(-1.0*K+KA)**-2.0+
          K*(-2.0*KA*(-1.0*K+KA)**-3.0*(-1.0*TERM3+TERM2)+
          KA*TIME*TERM3*(-1.0*K+KA)**-2.0))+
          KA*V**-1.0*(-1.0*K*TIME*TERM2*(-1.0*K+KA)**-1.0+
          K*(-1.0*K+KA)**-2.0*(-1.0*TERM3+TERM2)))
    D32 = -G3
    D33 = DOSE*KA*V**-1.0*(-1.0*K*TIME*(-1.0*K*TIME*TERM2*(-1.0*K+KA)**-1.0+
          K*TERM2*(-1.0*K+KA)**-2.0)+-1.0*K*TIME*TERM2*(-1.0*K+KA)**-1.0+
          K*(-1.0*K*TIME*TERM2*(-1.0*K+KA)**-2.0+
          2.0*K*(-1.0*K+KA)**-3.0*(-1.0*TERM3+TERM2))+
          K*(-1.0*K+KA)**-2.0*(-1.0*TERM3+TERM2))
  } else {
    D11 = 0
    D21 = 0
    D22 = 0
    D31 = 0
    D32 = 0
    D33 = 0
  }

  return(cbind(F, G1, G2, G3, H1, H2, D11, D21, D22, D31, D32, D33))
}

#########
METHOD = "ZERO" # PRED function refers this.
InitStep(DataAll, THETAinit=THETAinit, OMinit=OMinit, SGinit=SGinit, nTheta=nTheta, 
         LB=LB, UB=UB, METHOD=METHOD, Pred=PRED)
(EstRes = EstStep())            # About 3 secs
(CovRes = CovStep())            # About 1 sec
PostHocEta() # FinalPara from EstStep()

#########
METHOD = "COND" # PRED function refers this.
InitStep(DataAll, THETAinit=THETAinit, OMinit=OMinit, SGinit=SGinit, nTheta=nTheta, 
         LB=LB, UB=UB, METHOD=METHOD, Pred=PRED)
(EstRes = EstStep())            # About 4 mins
(CovRes = CovStep())            # About 40 secs
get("EBE", envir=e)

######### "LAPL" usually fails due to numerical difficulties.
METHOD = "LAPL" # PRED function refers this.
THETAinit = c(4, 50, 0.2) # It is changed for better convergence.
InitStep(DataAll, THETAinit=THETAinit, OMinit=OMinit, SGinit=SGinit, nTheta=nTheta, 
         LB=LB, UB=UB, METHOD=METHOD, Pred=PRED)
(EstRes = EstStep())            # About 3 mins, Succeeded with R-3.3.3 x64
(CovRes = CovStep())            # About 1 min
get("EBE", envir=e)
