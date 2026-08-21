\name{TrimOut}
\alias{TrimOut}
\title{Trimming and beautifying NONMEM original OUTPUT file}
\description{
TrimOut removes unnecessary parts from NONMEM original OUTPUT file.
}
\usage{
  TrimOut(inFile, outFile="PRINT.OUT")
}
\arguments{
  \item{inFile}{NONMEM original untidy OUTPUT file name}
  \item{outFile}{Output file name to be written}
}

\details{
NONMEM original OUTPUT file contains unnecessary parts such as CONTROL file content, Start/End Time, License Info, Print control characters such as "+", "0", "1". This function trims those.
}
\value{
  If \code{outFile} is \code{""}, the trimmed lines are returned as a character
  vector.  Otherwise the trimmed output is written to \code{outFile} in the
  current working (or designated) folder and \code{NULL} is returned invisibly.
}
\author{Kyun-Seop Bae <k@acr.kr>}
\keyword{NONMEM OUTPUT}

