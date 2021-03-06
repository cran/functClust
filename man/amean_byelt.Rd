% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calibrating.R
\name{amean_byelt}
\alias{amean_byelt}
\title{Arithmetic mean
 by components occurring within an assembly motif}
\usage{
amean_byelt(fobs, mOccur)
}
\arguments{
\item{fobs}{a numeric vector. The vector \code{fobs} contains the
quantitative performances of assemblages.}

\item{mOccur}{a matrix of occurrence (occurrence of components).
Its first dimension equals to \code{length(fobs)}. Its second dimension
 equals to the number of components.}
}
\value{
Return a vector of \code{length(fobs)}.
Its values are computed as the average
of mean performances of assemblages that contain the same components
as the assemblage to calibrate \code{opt.model = "byelt"} .
}
\description{
Take a vector \code{fobs}
of performances of assemblages
that share a same assembly motif,
and return a vector of performances
predicted as the arithmetic mean
of performances of assemblages
that contain the same components as the assemblage to predict. \cr
}
\details{
Modelled performances are computed
using arithmetic mean (\code{opt.mean = "amean"}) of performances.
Assemblages share a same assembly motif.
Modelled performances are the average
of mean performances of assemblages that contain the same components
as the assemblage to calibrate (\code{opt.model = "byelt"}).
This procedure corresponds to a linear model with each assembly motif
based on the component occurrence in each assemblage.
}
\seealso{
\code{\link{amean_byelt}} using arithmetic mean.
 \code{\link{gmean_byelt}} using geometric mean.
}
\keyword{internal}
