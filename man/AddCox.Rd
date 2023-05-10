\name{AddCox}
\alias{AddCox}
\title{Add a Covariate Column to an Existing NONMEM dataset}
\description{
A new covariate column can be added to an existing NONMEM dataset.
}
\usage{
  AddCox(nmData, coxData, coxCol, dateCol = "DATE", idCol = "ID")
}
\arguments{
  \item{nmData}{an existing NONMEM dataset}
  \item{coxData}{a data table containing a covariate column}
  \item{coxCol}{the covariate column name in the coxData table}
  \item{dateCol}{date column name in the NONMEM dataset and the covariate data table}
  \item{idCol}{ID column name in the NONMEM dataset and the covariate data table}
}

\details{
  It first carry forward for the missing data. If NA is remained, it carry backward.
}
\value{
  A new NONMEM dataset containing the covariate column
}
\author{Kyun-Seop Bae <k@acr.kr>}
\keyword{Data Preparation}
