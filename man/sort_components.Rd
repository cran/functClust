% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_fclust.R
\name{sort_components}
\alias{sort_components}
\title{Sort assembly centroids by decreasing or increasing mean performances}
\usage{
sort_components(fobs, mOccur, cols,
                pvalue  = fpvalue(),
                opt.dec = TRUE  )
}
\arguments{
\item{fobs}{numeric vector.
The vector contains observed assemblage performances.}

\item{mOccur}{a matrix of occurrence of components within the assemblages.
\code{dim(mOccur)[1] = length(fobs)}.
\code{dim(mOccur)[2]} equals to the total number of components.}

\item{cols}{an integer vector of \code{length(fobs)}.
The vector contains a colour label for each assemblages.}

\item{pvalue}{a probability,
used as threshold
in the variance analysis. Then \code{pvalue} must be
higher than \code{0} and lower than \code{1}.
Groups significantly different
(at a p-value < \code{pvalue}) are then indicated by differents letters
on the right of boxplots.}

\item{opt.dec}{a logical.
If \code{opt.dec = TRUE},
assembly motifs are sorted by decreasing mean performances.
If \code{opt.dec = FALSE},
assembly motifs are sorted by increasing mean performances.}
}
\value{
Return a table
containing statistical properties of assembly centroids,
sorted by increasing or decreasing mean performances.
}
\description{
Plot
}
\details{
None.
}
\keyword{internal}
