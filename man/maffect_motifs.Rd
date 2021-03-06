% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/labelling.R
\name{maffect_motifs}
\alias{maffect_motifs}
\title{Labelling of a whole tree of component clustering by assembly motif}
\usage{
maffect_motifs(tree, mOccur)
}
\arguments{
\item{tree}{an integer square-matrix.
The matrix represents a hierarchical tree of species clustering.
A whole tree of component clustering is a set of component affectations.}

\item{mOccur}{a matrix of occurrence (occurrence of components).
Its first dimension equals to \code{length(fct)}. Its second dimension
 equals to the number of components.}
}
\value{
Return a matrix of labels, of \code{dim(tree$aff)}.
}
\description{
A whole tree of component clustering is a set of component affectations.
Take successively each vector of component affectation by clusters,
then look for the assemblages containing components
that belong to the same clusters,
finally label each set of assemblages by a different number.
}
\details{
A hierarchical tree is recorded
as a list of matrix of component affectation
and a vector of coefficients of determination.
The matrix of component affectation is a square matrix.
Its dimension is "number of leaves" = "number of components".
The clusters of components are labelled by integers, an integer by cluster.
The first line of the square matrix corresponds to the trivial cluster.
The trivial cluster get together all the components:
all components are thus labelled by one.
The second line separates the trivial cluster into 2 sub-clusters:
the components are labelled 1 or 2
according to the cluster to which they belong.
The following lines separate the trivial cluster into n sub-clusters,
labelled from 1 to the line numero.
The last line separate each component in a singleton cluster.\cr

A set of \code{s} components generates \code{m=2^s-1 assembly motifs}.
The function only looks for existing assembly motifs.
Each assembly motif is labelled from \code{1} to \code{m},
with \code{m <= s}.
The assembly motifs are first sorted by size (number of clusters),
then labelled gradually:
the first motif identified is labelled \code{1},
the second motif identified is labelled \code{2}, etc...\cr

The function is used to plot tree,
by sorting clusters from the most efficient cluster (on the left),
its size is generally smaller than the least efficient clusters,
towards the least efficient clusters (on the right).
}
\keyword{internal}
