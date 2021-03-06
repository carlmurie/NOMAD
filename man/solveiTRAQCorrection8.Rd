% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NOMAD.R
\name{solveiTRAQCorrection8}
\alias{solveiTRAQCorrection8}
\title{Correct peptide abundances based on iTRAQ specification. The iTRAQTable
below characterizes the components of the abundance scores.}
\usage{
solveiTRAQCorrection8(y, iTRAQTable = NULL)
}
\arguments{
\item{y}{\code{logical} A vector of 8 peptide abundances ordered by iTRAQ (1 to 8)
for a single protein}

\item{iTRAQTable}{\code{logical}  A 8x8 matrix of isotope corrections. Probably shouldn't modify
unless you have a very good reason.}
}
\value{
vector of peptide abundances after correcting for isotope bleed over.
}
\description{
Correct peptide abundances based on iTRAQ specification. The iTRAQTable
below characterizes the components of the abundance scores.
}
