% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validating_loo.R
\name{validate_gmean_byelt_LOO_xpr}
\alias{validate_gmean_byelt_LOO_xpr}
\title{Geometric mean
 by elements occurring within assembly motif
 using leave-one-out method
 over several experiments}
\usage{
validate_gmean_byelt_LOO_xpr(fobs, assMotif, mOccur, xpr)
}
\arguments{
\item{fobs}{a numeric vector. The vector \code{fobs} contains the
quantitative performances of assemblages.}

\item{assMotif}{a vector of labels of \code{length(fobs)}.
The vector \code{assMotif} contains the assembly motifs of assemblages.}

\item{mOccur}{a matrix of occurrence (occurrence of elements).
Its first dimension equals to \code{length(fobs)}. Its second dimension
 equals to the number of elements.}

\item{xpr}{a vector of numerics of \code{length(fobs)}.
The vector \code{xpr} contains the weight of each experiment,
and the labels (in \code{names(xpr)}) of different experiments.
The weigth of each experiment is used
in the computation of the Residual Sum of Squares
in the function \code{rss_clustering}.
The used formula is \code{rss}
if each experiment has the same weight.
The used formula is \code{wrss}
(barycenter of RSS for each experiment)
if each experiment has different weights.
All assemblages that belong to a given experiment
should then have a same weigth.
Each experiment is identified by its names (\code{names(xpr)})
and the RSS of each experiment is weighted by values of \code{xpr}.
The vector \code{} is generated
by the function \code{stats::setNames}.}
}
\value{
Return a vector of \code{length(fobs)}.
Its values are computed as the average
of mean performances of all assemblages,
experiment by experiment,
that contain the same elements
as the assemblage to predict (\code{opt.model = "byelt"}) .
}
\description{
Take a vector \code{fobs}
of performances of assemblages
over several experiments
that share a same assembly motif,
and return a vector of performances
predicted as the geometric mean
of performances of assemblages
that contain the same elements as the assemblage to predict.\cr

All assemblages share the same assembly motif.
}
\details{
Predicted performances are computed
using geometric mean (\code{opt.mean = "gmean"})
of performances.
Assemblages share a same assembly motif (\code{opt.model = "bymot"}).
Predicted performances are the average,
experiment by experiment,
of mean performances of assemblages
that contain the same elements
as the assemblage to predict,
except the assemblage to predict.
This procedure corresponds to a linear model with each assembly motif
based on the element occurrence in each assemblage.
}
\seealso{
\code{\link{validate_amean_bymot_LOO_xpr}} arithmetic mean
  byassembly motif \cr
 \code{\link{validate_amean_byelt_LOO_xpr}} arithmetic mean
  by elements occurring within assembly motif \cr
 \code{\link{validate_gmean_bymot_LOO_xpr}} geometric mean
  by assembly motif \cr
}
\keyword{internal}
