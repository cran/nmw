\name{TrimOut}
\alias{TrimOut}
\title{Trimming and beutifying NONMEM original OUTPUT file}
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
  outFile will be written in the current working folder or designated folder.
  Ths returns TRUE if the process was smooth.
}
\author{Kyun-Seop Bae <k@acr.kr>}
\keyword{NONMEM OUTPUT}

