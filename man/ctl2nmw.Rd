\name{ctl2nmw}
\alias{ctl2nmw}
\title{Generate nmw PRED and Initial Estimates from a NONMEM Control Stream}
\description{
Reads a closed-form \code{$PRED} NONMEM control stream and auto-generates the
\code{PRED} function (with its \code{deriv()} companions) and the initial
estimates (\code{THETAinit}, \code{OMinit}, \code{SGinit}, \code{LB}, \code{UB},
\code{METHOD}) required by \code{\link{InitStep}}.  This removes the need to
hand-write the \code{PRED} function shown in the package examples.
}
\usage{
ctl2nmw(ctl, eval = FALSE, dose = NULL, file = NULL)
}
\arguments{
  \item{ctl}{Path to a NONMEM control file, or a character vector/string of its lines.}
  \item{eval}{If \code{TRUE}, also build the live \code{PRED} closure (runs \code{deriv()}); if \code{FALSE} (default) only the code text and numeric initial estimates are returned.}
  \item{dose}{Optional numeric to override a constant \code{DOSE} assignment found in \code{$PRED}.}
  \item{file}{Optional path; if given, the generated R code is written there.}
}
\details{
The structural prediction \code{F} is flattened to a single expression by
inlining intermediate assignments (e.g. \code{KA}, \code{V}, \code{K}), because
\code{\link[stats]{deriv}} does not accept multi-statement blocks.  \code{FGD}
differentiates \code{F} with respect to every \code{ETA}; \code{H} differentiates
the error model \code{Y} with respect to every \code{EPS} while holding \code{F}
symbolic.  Free symbols are classified over both \code{F} and the error model, so
\code{THETA}s or data items that appear only in the residual error (for example a
combined additive-plus-proportional \code{W} term) are wired into \code{H}
correctly.

The generated \code{PRED} returns a matrix with columns \code{F}, \code{G1..Gk}
(with respect to \code{ETA}), \code{H1..Hm} (with respect to \code{EPS}) and, for
\code{METHOD = "LAPL"}, the second-derivative columns \code{D11, D21, D22, ...} in
the order \code{\link{InitStep}} expects.

nmw addresses data columns by name, so the supplied \code{DataAll} must contain
columns named \code{ID}, \code{TIME}, \code{DV} and any covariate that \code{$PRED}
references, independent of the \code{$INPUT} order.

\bold{Scope (v1).}  Only a single \code{$PROBLEM} with a closed-form \code{$PRED}
block is supported.  ODE models (\code{$PK}/\code{$DES}/\code{$SUBROUTINE}/
\code{ADVAN}), \code{IF}/\code{THEN} branching, MU-referencing, mixtures and
\code{FIX}ed parameters raise an informative error rather than generating
incorrect code.  Because definitions are inlined, assigning the same symbol more
than once in \code{$PRED} (sequential reassignment) or referencing a symbol
before it is defined also raises an informative error.  A diagonal or
block-diagonal \code{$OMEGA} is promoted to a full block (nmw estimates the full
lower triangle) with a note, and a non-diagonal \code{$SIGMA} is reduced to its
diagonal with a note.

\bold{Theta bounds.}  \code{\link{InitStep}} scales thetas with a logit-type
transform, so every initial estimate must lie strictly inside its bounds.  An
unbounded \code{$THETA} entry receives the bounds \code{(0, 1e6)} when its
initial estimate is positive, and \code{(-1e6, 1e6)} otherwise (with a note); an
initial estimate on or outside explicit \code{$THETA} bounds raises an
informative error.
}
\value{
An object of class \code{nmwCode} (a list); \code{print} shows the generated code.
  \item{code}{Character scalar: the generated R script (initial estimates + \code{PRED}).}
  \item{PRED}{The prediction closure when \code{eval = TRUE}, otherwise \code{NULL}.}
  \item{THETAinit, LB, UB}{Numeric initial estimate and bounds for theta.}
  \item{OMinit, SGinit}{Omega and (diagonal) sigma initial matrices.}
  \item{METHOD}{One of \code{"ZERO"}, \code{"COND"}, \code{"LAPL"}.}
  \item{nTheta, nEta, nEps}{Parameter counts.}
  \item{input}{Parsed \code{$INPUT} column information.}
  \item{notes}{Character vector of fidelity warnings (may be empty).}
}
\references{
NONMEM Users Guide
}
\author{Kyun-Seop Bae <k@acr.kr>}
\seealso{\code{\link{InitStep}}, \code{\link{EstStep}}, \code{\link{CovStep}}}
\examples{
ctl <- system.file("doc", "THEO-FOCEI.CTL", package = "nmw")
gen <- ctl2nmw(ctl)
gen                       # print the generated PRED + initial estimates
gen$THETAinit
gen$METHOD

\dontrun{
DataAll = Theoph
colnames(DataAll) = c("ID", "BWT", "DOSE", "TIME", "DV")
DataAll[,"ID"] = as.numeric(as.character(DataAll[,"ID"]))

gen = ctl2nmw(ctl, eval = TRUE)
InitStep(DataAll, THETAinit = gen$THETAinit, OMinit = gen$OMinit, SGinit = gen$SGinit,
         LB = gen$LB, UB = gen$UB, Pred = gen$PRED, METHOD = gen$METHOD)
EstStep()
CovStep()
}
}
\keyword{Population Pharmacokinetics}
