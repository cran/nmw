\name{InitStep}
\alias{InitStep}
\title{Initialization Step}
\description{
It receives parameters for the estimation and stores them into \code{e} environment.
}
\usage{
InitStep(DataAll, THETAinit, OMinit, SGinit, LB, UB, Pred, METHOD)
}
\arguments{
  \item{DataAll}{Data for all subjects. It should contain columns which \code{Pred} function uses.}
  \item{THETAinit}{Theta initial values}
  \item{OMinit}{Omega matrix initial values}
  \item{SGinit}{Sigma matrix initial values}
  \item{LB}{Lower bounds for theta vector}
  \item{UB}{Upper bounds for theta vector}
  \item{Pred}{Prediction function name}
  \item{METHOD}{one of the estimation methods \code{"ZERO"}, \code{"COND"}, or \code{"LAPL"}}
}
\details{
Prediction function should return not only prediction values(F or IPRED) but also G (first derivative with respect to etas) and H (first derivative of Y with respect to epsilon).
For the \code{"LAPL"}, prediction function should return second derivative with respect to eta also.
"INTERACTION" is \code{TRUE} for \code{"COND"} and \code{"LAPL"} option, and \code{FALSE} for \code{"ZERO"}.
Omega matrix should be full block one.
Sigma matrix should be diagonal one.
}
\value{
This does not return values, but stores necessary values into the environment \code{e}.
}
\references{
NONMEM Users Guide
}
\author{Kyun-Seop Bae <k@acr.kr>}
\examples{
DataAll = Theoph
colnames(DataAll) = c("ID", "BWT", "DOSE", "TIME", "DV")
DataAll[,"ID"] = as.numeric(as.character(DataAll[,"ID"]))

nTheta = 3
nEta = 3
nEps = 2

THETAinit = c(2, 50, 0.1) # Initial estimate
OMinit = matrix(c(0.2, 0.1, 0.1, 0.1, 0.2, 0.1, 0.1, 0.1, 0.2), nrow=nEta, ncol=nEta)
OMinit
SGinit = diag(c(0.1, 0.1))
SGinit

LB = rep(0, nTheta) # Lower bound
UB = rep(1000000, nTheta) # Upper bound

FGD = deriv(~DOSE/(TH2*exp(ETA2))*TH1*exp(ETA1)/(TH1*exp(ETA1) - TH3*exp(ETA3))*
             (exp(-TH3*exp(ETA3)*TIME)-exp(-TH1*exp(ETA1)*TIME)),
            c("ETA1","ETA2","ETA3"),
            function.arg=c("TH1", "TH2", "TH3", "ETA1", "ETA2", "ETA3", "DOSE", "TIME"),
            func=TRUE, hessian=TRUE)
H = deriv(~F + F*EPS1 + EPS2, c("EPS1", "EPS2"), function.arg=c("F", "EPS1", "EPS2"), func=TRUE)

PRED = function(THETA, ETA, DATAi)
{
  FGDres = FGD(THETA[1], THETA[2], THETA[3], ETA[1], ETA[2], ETA[3], DOSE=320, DATAi[,"TIME"])
  Gres = attr(FGDres, "gradient")
  Hres = attr(H(FGDres, 0, 0), "gradient")

  if (e$METHOD == "LAPL") {
    Dres = attr(FGDres, "hessian")
    Res = cbind(FGDres, Gres, Hres, Dres[,1,1], Dres[,2,1], Dres[,2,2], Dres[,3,])
    colnames(Res) = c("F", "G1", "G2", "G3", "H1", "H2", "D11", "D21", "D22", "D31", "D32", "D33")
  } else {
    Res = cbind(FGDres, Gres, Hres)
    colnames(Res) = c("F", "G1", "G2", "G3", "H1", "H2")
  }
  return(Res)
}

######### First Order Approximation Method
InitStep(DataAll, THETAinit=THETAinit, OMinit=OMinit, SGinit=SGinit, LB=LB, UB=UB,
         Pred=PRED, METHOD="ZERO")

######### First Order Conditional Estimation with Interaction Method
InitStep(DataAll, THETAinit=THETAinit, OMinit=OMinit, SGinit=SGinit, LB=LB, UB=UB,
         Pred=PRED, METHOD="COND")

######### Laplacian Approximation with Interacton Method
InitStep(DataAll, THETAinit=THETAinit, OMinit=OMinit, SGinit=SGinit, LB=LB, UB=UB,
         Pred=PRED, METHOD="LAPL")
}
\keyword{Initialization Step}
