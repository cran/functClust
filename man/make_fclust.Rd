% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fclust.R
\name{make_fclust}
\alias{make_fclust}
\title{Make the formatted result from a functional clustering}
\usage{
make_fclust(res, opt.na, opt.repeat, affectElt)
}
\arguments{
\item{res}{an object resulting from a functional clustering.}

\item{opt.na}{a logical.
The records for each assemblage can have \code{NA}
in matrix of occurrence or in observed assemblage performances.
If \code{opt.na = FALSE} (by default), an error is returned.
If \code{opt.na = TRUE}, the records with \code{NA} are ignored.}

\item{opt.repeat}{a logical.
in any case, the function looks for
different assemblages with identical elemental composition.
Messages indicate these identical assemblages.
If \code{opt.repeat = FALSE} (by default),
their performances are averaged.
If \code{opt.repeat = TRUE}, nothing is done,
and the data are processed as they are.}

\item{affectElt}{a vector of characters or integers,
as long as the number of components \code{nbElt},
that indicates the labels of different functional clusters
to which each component belongs.
Each functional cluster is labelled as a character or an integer, and
each component must be identified by its name in \code{names(affectElt)}.
The number of functional clusters defined in \code{affectElt}
determines an \emph{a priori} level of component clustering
(\code{level <- length(unique(affectElt))}).\cr

If \code{affectElt = NULL} (by default),
the option \code{opt.method} must be specified.
If \code{affectElt} is specified,
the option \code{opt.method} switchs to \code{apriori}.}
}
\value{
\code{fres}, which is an object
    resulting from a functional clustering.
}
\description{
Check and homogeneize the format of results
   from a functional clustering,
   then make the formatted result.
}
\details{
The results are saved in
  the file \code{"filename.options.csv"}.
  If the file does not exist or is corrupted, the function is stopped.
}
\keyword{internal}
