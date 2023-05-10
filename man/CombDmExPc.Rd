\name{CombDmExPc}
\alias{CombDmExPc}
\title{Combine the demographics(DM), dosing(EX), and DV(PC) tables into a new NONMEM dataset}
\description{
A new NONMEM dataset can be created from the demographics, dosing, and DV tables.
}
\usage{
  CombDmExPc(dm, ex, pc)
}
\arguments{
  \item{dm}{A demographics table. It should contain a row per subject.}
  \item{ex}{An exposure table. Drug administration (dosing) history table.}
  \item{pc}{A DV(dependent variable) or PC(drug concentration) table}
}
\details{
  Combining a demographics, a dosing, and a concentration table can produce a new NONMEM dataset.
}
\value{
  A new NONMEM dataset
}
\author{Kyun-Seop Bae <k@acr.kr>}
\keyword{Data Preparation}
