\name{EstStep}
\alias{EstStep}
\title{Estimation Step}
\description{This estimates upon the conditions with \code{InitStep}.}
\usage{EstStep()}
\details{
It does not have arguments.
All necessary arguments are stored in the \code{e} environment.
It assumes "INTERACTION" between eta and epsilon for \code{"COND"} and \code{"LAPL"} options.
The output is basically same with NONMEM output.
}
\value{
  \item{Initial OFV}{initial value of objective function}
  \item{Time}{time consumed for this step}
  \item{Optim}{the raw output from \code{optim} function}
  \item{Final Estimates}{final estimates in the original scale}
}
\references{NONMEM Users Guide}
\seealso{\code{\link{InitStep}}}
\author{Kyun-Seop Bae <k@acr.kr>}
\examples{
# Only After InitStep
#EstStep()
}
\keyword{Estimation Step}
