% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calibrating.R
\name{calibrate_byminrss}
\alias{calibrate_byminrss}
\title{Modelling of performances of assemblages}
\usage{
calibrate_byminrss(fobs, assMotif, mOccur, xpr,
opt.mean = "amean", opt.model = "bymot"  )
}
\arguments{
\item{fobs}{a numeric vector. The vector \code{fobs} contains the
quantitative performances of assemblages.}

\item{assMotif}{a vector of labels of \code{length(fobs)}.
The vector \code{assMotif} contains the assembly motifs of assemblages.}

\item{mOccur}{a matrix of occurrence (occurrence of components).
Its first dimension equals to \code{length(fobs)}. Its second dimension
 equals to the number of components.}

\item{xpr}{a vector of numerics of \code{length(fobs)}.
The vector \code{xpr} contains the weight of each experiment,
and the labels (in \code{names(xpr)}) of different experiments.
The weigth of each experiment is used
in the computation of the Residual Sum of Squares
in the function \code{rss_clustering}.
All assemblages that belong to a given experiment
should then have a same weigth.
Each experiment is identified by its names (\code{names(xpr)})
and the RSS of each experiment is weighted by values of \code{xpr}.
The vector \code{xpr} is generated
by the function \code{stats::setNames}.}

\item{opt.mean}{a character equals to \code{"amean"} or \code{"gmean"}.
Switchs to arithmetic formula if \code{opt.mean = "amean"}.
Switchs to geometric formula if \code{opt.mean = "gmean"}.}

\item{opt.model}{a character equals to \code{"bymot"} or \code{"byelt"}.
Switchs to simple mean by assembly motif if \code{opt.model = "bymot"}.
Switchs to linear model with assembly motif if \code{opt.model = "byelt"}.}
}
\value{
Return a vector of \code{length(fobs)}.
Its values are computed according to \code{opt.mean} and \code{opt.model}.
}
\description{
Take a vector \code{fobs} of assembly performances
over several experiments
and return a vector of performances
predicted as the mean performances of assemblages
that share the same assembly motif.\cr

Assembly motifs are labelled in the vector \code{assMotif}.
Experiments are labelled in the vector \code{xpr}.
Modelling options are indicated in \code{opt.mean} and \code{opt.model}.
Occurrence matrix \code{mOccur} is used if \code{opt.model =  "byelt"}.
}
\details{
Modelled performances are computed
using arithmetic mean (\code{opt.mean = "amean"})
or geometric mean (\code{opt.mean = "gmean"}).\cr

If \code{opt.model = "bymot"},
modelled performances are means
of performances of assemblages
that share a same assembly motif
by including all assemblages that belong to a same assembly motif. \cr

If \code{opt.model = "byelt"},
modelled performances are the average
of mean performances of assemblages
that share a same assembly motif
and that contain the same components
as the assemblage to calibrate.
This procedure corresponds to a linear model within each assembly motif
based on the component occurrence in each assemblage.
If no assemblage contains component belonging to assemblage to calibrate,
performance is the mean performance of all assemblages
as in \code{opt.model = "bymot"}.
}
\seealso{
\code{\link{validate_using_cross_validation}}
   predicts performances of assemblages.
}
\keyword{internal}
