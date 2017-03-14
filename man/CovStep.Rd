\name{CovStep}
\alias{CovStep}
\title{Covariance Step}
\description{
It calculates standard erros and various variance matrices with the \code{e$FinalPara} after estimation step.
}
\usage{CovStep()}
\details{
Because \code{EstStep} uses nonlinear optimization, covariance step is separated from estimation step.
It caculcates variance-covariance matrix of estimates on the original scale.
}
\value{
  \item{Time}{consumed time}
  \item{Standard Error}{standard error of the estimates in the order of theta, omega, and sigma}
  \item{Covariance Matrix of Estimates}{covariance matrix of estimates in the order of theta, omega, and sigma. This is inverse(R) x S x inverse(R) by default.}
  \item{Correlation Matrix of Estimates}{correlation matrix of estimates in the order of theta, omega, and sigma}
  \item{Inverse Covariance Matrix of Estimates}{inverse covariance matrix of estimates in the order of theta, omega, and sigma}
  \item{Eigen Values}{eigen values of covariance matrix}
  \item{R Matrix}{R matrix of NONMEM, second derivative of log likelihood function with respect to esimation parameters}
  \item{S Matrix}{S matrix of NONMEM, sum of individual cross-product of first derivative of log likelihood function with respect to esimation parameters}
}
\references{NONMEM Users Guide}
\author{Kyun-Seop Bae <k@acr.kr>}
\seealso{
\code{\link{EstStep}}, \code{\link{InitStep}} 
}
\examples{
# Only after InitStep and EstStep
#CovStep()
}
\keyword{Covariance Step}
