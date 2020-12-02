#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                            LABELLING.R
#
#  Functions to label assemblages and assembly motifs                      ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#' @include
#'  tools.R
#'
NULL



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Compact a vector of indices
#'
#' @description Take a vector \code{v} of indices,
#' then simplify the indexing system
#' and compact it from \code{1} to \code{length(unique(v))}.
#'
#' @usage compact_index(v)
#'
#' @param v a vector of integers, used as indices.
#'
#' @details The indexing vector \code{v} is simplified
#' and compacts from \code{1} to \code{length(unique(v))}.
#' The index order is not changed, only the labels are changed.
#' The function is useful
#' for debugging computations on values sorted by a program.
#'
#' @return Return a vector of integer,
#' that can be browsed by the code \code{seq_len(length(unique(v)))}.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

compact_index <- function(v) {

  res <- v
  set <- sort(unique(v))
  for (i in seq_along(set)) res[v == set[i]] <- i

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Renumber a vector of component affectation
#'
#' @description Take a vector of component affectation
#' and change the cluster label \code{1} as the highest cluster label.
#'
#' @usage shift_affectElt(affectElt)
#'
#' @inheritParams affect_motifs
#'
#' @details A hierarchical tree is recorded
#' as a list of matrix of component affectation
#' and a vector of coefficients of determination.
#' The matrix of component affectation is a square matrix.
#' Its dimension is "number of leaves" = "number of components".
#' The clusters of components are labelled by integers, an integer by cluster.
#' The first line of the square matrix corresponds to the trivial cluster.
#' The trivial cluster get together all the components:
#' all components are thus labelled by one.
#' The second line separates the trivial cluster into 2 sub-clusters:
#' the components are labelled 1 or 2
#' according to the cluster to which they belong.
#' The following lines separate the trivial cluster into n sub-clusters,
#' labelled from 1 to the line numero.
#' The last line separate each component in a singleton cluster.\cr
#'
#' The initial trivial cluster is labelled \code{1} .
#' Consequently, the largest cluster remains labelled \code{1} .
#' The function \code{shift_affectElt} changes the cluster label \code{1}
#' as the highest cluster label.
#' Then decreases each other clusters number by \code{1}.\cr
#'
#' The function is used to plot tree,
#' by sorting clusters from the most efficient cluster (on the left),
#' its size is generally smaller than the least efficient clusters,
#' towards the least efficient clusters (on the right).
#'
#' @return Return a vector of component affectation.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

shift_affectElt <- function(affectElt) {

  affectElt[affectElt == 1] <- max(affectElt) + 1L
  affectElt <- affectElt - 1L

  return(affectElt)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Labelling (by lowercase letters) of components by cluster
#'
#' @description Take a vector of component affectation by clusters,
#' then name the components by cluster.
#'
#' @usage name_clusters(affectElt)
#'
#' @inheritParams affect_motifs
#'
#' @details
#' General:\cr
#' A hierarchical tree is recorded
#' as a list of matrix of component affectation
#' and a vector of coefficients of determination.
#' The matrix of component affectation is a square matrix.
#' Its dimension is "number of leaves" = "number of components".
#' The clusters of components are labelled by integers, an integer by cluster.
#' The first line of the square matrix corresponds to the trivial cluster.
#' The trivial cluster get together all the components:
#' all components are thus labelled by one.
#' The second line separates the trivial cluster into 2 sub-clusters:
#' the components are labelled 1 or 2
#' according to the cluster to which they belong.
#' The following lines separate the trivial cluster into n sub-clusters,
#' labelled from 1 to the line numero.
#' The last line separate each component in a singleton cluster.\cr
#'
#'  A set of \code{s} components generates \code{m = 2^s - 1 assembly motifs}.
#' The clusters of components are previously shifted by \code{shift_affectElt},
#' for sorting clusters from the most efficient cluster
#'  (labelled \code{a}, on the left),
#' its size is generally smaller than the least efficient clusters,
#' towards the least efficient clusters
#'  (labelled \code{m' < m}, on the right).
#'  component clusters are thus labelled: \code{a, b, c, ... m' < m}.
#'
#' Each assembly motif is a combination of clusters.
#' The assembly motifs are simply labelled
#' as the combination of clusters that they contain:
#' \code{a, b, c, ..., ab, ac bc, ..., abc ...}.
#'
#' The order of all vectors of strings is
#' this of the input vector \code{affectElt}.
#'
#' @return Return a vector of strings
#' corresponding to names of components by clusters.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

name_clusters <- function(affectElt) {

  return(fletters()[shift_affectElt(affectElt)])
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Labelling (by lowercase letters) of assemblages by assembly motif
#'
#' @description Take a vector of component affectation by clusters,
#' then look for the assemblages containing components
#' that belong to the same clusters,
#' finally label each set of assemblages by a different number.
#'
#' @usage name_motifs(affectElt, mOccur)
#'
#' @inheritParams affect_motifs
#'
#' @details
#' General:\cr
#' A hierarchical tree is recorded
#' as a list of matrix of component affectation
#' and a vector of coefficients of determination.
#' The matrix of component affectation is a square matrix.
#' Its dimension is "number of leaves" = "number of components".
#' The clusters of components are labelled by integers, an integer by cluster.
#' The first line of the square matrix corresponds to the trivial cluster.
#' The trivial cluster get together all the components:
#' all components are thus labelled by one.
#' The second line separates the trivial cluster into 2 sub-clusters:
#' the components are labelled 1 or 2
#' according to the cluster to which they belong.
#' The following lines separate the trivial cluster into n sub-clusters,
#' labelled from 1 to the line numero.
#' The last line separate each component in a singleton cluster.\cr
#'
#'  A set of \code{s} components generates \code{m = 2^s - 1 assembly motifs}.
#' The clusters of components are previously shifted by \code{shift_affectElt},
#' for sorting clusters from the most efficient cluster
#'  (labelled \code{a}, on the left),
#' its size is generally smaller than the least efficient clusters,
#' towards the least efficient clusters
#'  (labelled \code{m' < m}, on the right).
#'  component clusters are thus labelled: \code{a, b, c, ... m' < m}.
#'
#' Each assembly motif is a combination of clusters.
#' The assembly motifs are simply labelled
#' as the combination of clusters that they contain:
#' \code{a, b, c, ..., ab, ac bc, ..., abc ...}.
#'
#' The order of all vectors of strings is
#' this of the input vector \code{affectElt}.
#'
#' @return Return a list of vector of strings.\cr
#'
#' \itemize{
#'  \item \code{components}: corresponds to \code{colnames(mOccur)},
#'  \item \code{clusters}: corresponds to names of clusters of components,
#'  \item \code{assemblages}: corresponds to \code{rownames(mOccur)},
#'  \item \code{motifs}: corresponds to names of assembly motifs,
#'       termed as combinations of component clusters.
#'  }
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

name_motifs <- function(affectElt, mOccur) {

  # label the component clusters
  nomAffect <- name_clusters(affectElt)

  # label the assembly motifs
  nbAss    <- dim(mOccur)[1]
  namMotif <- character(nbAss)
  for (ass in seq_len(nbAss))
    namMotif[ass] <- paste(sort(unique(nomAffect[as.logical(mOccur[ass, ])])),
                           collapse = "")

  return(namMotif)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Label assemblages by assembly motif
#'
#' @description Take a vector of affectation of components
#' by functional groups,
#' then build the different assembly motifs
#' generated by combining the different functional groups,
#' finally label the assemblages
#' that belong to the same assembly motifs.
#' The labels are integers.
#'
#' @usage affect_motifs(affectElt, mOccur)
#'
#' @param affectElt a vector of integers.
#' The vector contains the labels of the different functional groups
#' to which each component belongs.
#' Each functional group is labelled as an integer.
#'
#' @param mOccur a matrix of occurrence (occurrence of components).
#'  Its first dimension equals to \code{length(fct)}.
#'   Its second dimension
#'   equals to the number of components.
#'
#' @details A set of \code{s} components generates
#' \code{m = 2^s - 1} assembly motifs.
#' The function only looks for existing, observed assembly motifs.
#' Each assembly motif is labelled from \code{1} to \code{m},
#' with \code{m <= s}.
#' The assembly motifs are first sorted by size (number of clusters),
#' then labelled gradually:
#' the first motif identified is labelled \code{1},
#' the second motif identified is labelled \code{2}, etc...
#' The number of observed assembly motifs is lower or equal to \code{m}.
#'
#' @return Return a vector of labels, of \code{dim(mOccur)[2]}.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

affect_motifs <- function(affectElt, mOccur) {

  namMotifs <- name_motifs(affectElt, mOccur)
  setMotifs <- unique(namMotifs)

  assMotifs <- integer(length(namMotifs))
  for (mot in seq_along(setMotifs))
    assMotifs[namMotifs == setMotifs[mot]] <- mot

  return(assMotifs)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Labelling of a whole tree of component clustering by assembly motif
#'
#' @description
#' A whole tree of component clustering is a set of component affectations.
#' Take successively each vector of component affectation by clusters,
#' then look for the assemblages containing components
#' that belong to the same clusters,
#' finally label each set of assemblages by a different number.
#'
#' @usage maffect_motifs(tree, mOccur)
#'
#' @param tree an integer square-matrix.
#'  The matrix represents a hierarchical tree of species clustering.
#'  A whole tree of component clustering is a set of component affectations.
#'
#' @param mOccur a matrix of occurrence (occurrence of components).
#'  Its first dimension equals to \code{length(fct)}. Its second dimension
#'   equals to the number of components.
#'
#' @details
#' A hierarchical tree is recorded
#' as a list of matrix of component affectation
#' and a vector of coefficients of determination.
#' The matrix of component affectation is a square matrix.
#' Its dimension is "number of leaves" = "number of components".
#' The clusters of components are labelled by integers, an integer by cluster.
#' The first line of the square matrix corresponds to the trivial cluster.
#' The trivial cluster get together all the components:
#' all components are thus labelled by one.
#' The second line separates the trivial cluster into 2 sub-clusters:
#' the components are labelled 1 or 2
#' according to the cluster to which they belong.
#' The following lines separate the trivial cluster into n sub-clusters,
#' labelled from 1 to the line numero.
#' The last line separate each component in a singleton cluster.\cr
#'
#' A set of \code{s} components generates \code{m=2^s-1 assembly motifs}.
#' The function only looks for existing assembly motifs.
#' Each assembly motif is labelled from \code{1} to \code{m},
#' with \code{m <= s}.
#' The assembly motifs are first sorted by size (number of clusters),
#' then labelled gradually:
#' the first motif identified is labelled \code{1},
#' the second motif identified is labelled \code{2}, etc...\cr
#'
#' The function is used to plot tree,
#' by sorting clusters from the most efficient cluster (on the left),
#' its size is generally smaller than the least efficient clusters,
#' towards the least efficient clusters (on the right).
#'
#' @return Return a matrix of labels, of \code{dim(tree$aff)}.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

maffect_motifs <- function(tree, mOccur) {

  nbclTot <- length(tree$cor)
  nbclMax <- which(tree$cor == max(tree$cor))[1]

  massMotif <- matrix(0L, nrow = nbclTot, ncol = dim(mOccur)[1],
                       dimnames = list(seq_len(nbclTot), rownames(mOccur)))
  storage.mode(massMotif) <- "integer"

  for (nbcl in seq_len(nbclMax)) {
    affectElt         <- compact_index(tree$aff[nbcl, ])
    massMotif[nbcl, ] <- affect_motifs(affectElt, mOccur)
  }

  if (nbclMax < nbclTot)
    for (nbcl in (nbclMax + 1L):nbclTot)
      massMotif[nbcl, ] <- massMotif[nbclMax, ]

  return(massMotif)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# End of file                                                              ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
