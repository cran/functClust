% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_fclust.R
\name{plot_motifs_infos}
\alias{plot_motifs_infos}
\title{Plot reference graphs for checking that plot of sorted data are right}
\usage{
plot_motifs_infos(res, opt.sort = "performance")
}
\arguments{
\item{res}{the result of a functional clustering
obtained using the function \code{\link{fclust}}.}

\item{opt.sort}{a string,
that specifies the way for sorting the motifs.
It can be \code{"performance"} or \code{"name"},
indicating a sorting by motif performances,
or a sorting by motif names.}
}
\value{
Nothing. It is a procedure.
}
\description{
Plot two reference graphs for checking
that plot of sorted data are right:
  the non-sorted assembly motifs
  and assembly motifs sorted by decreasing mean observed performances.
}
\details{
This function was useful only for setting up the R-code.
It is now useful for users to check the resulst,
and be sure to the goodness of plotting ...
The written values are (form up to bottom):
names of assembly motifs, their effectifs,
the mean observed performance, the order
and the symbols (colour x symbol) systematically associated
with each assembly motif in all plots produced
by the package \code{functClust}.
}
\keyword{internal}
