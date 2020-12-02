#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                        VALIDATING_LOO.R
#
# Functions for predicting par leave-one-out                               ####
#      the performances of assemblages
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#' @include
#'  stats.R
#'
NULL


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Arithmetic mean by assembly motif                                        ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic mean
#' by assembly motif
#' using leave-one-out method
#'
#' @description Take a vector \code{fobs} of assembly performances
#' and return a vector of performances
#' predicted as the arithmetic mean
#' of performances of all assemblages. \cr
#'
#' All assemblages share the same assembly motif.
#'
#' @usage amean_bymot_LOO(fobs)
#'
#' @inheritParams validate_using_cross_validation
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the arithmetic mean of performances,
#' by excluding the assemblage to predict.
#'
#' @details Predicted performances are computed using arithmetic mean
#'   (\code{opt.mean = "amean"}) of performances of assemblages
#'   that share a same assembly motif (\code{opt.model = "bymot"})
#'   by excluding the assemblage to predict.
#'
#' @seealso
#' \code{\link{amean_byelt_LOO}} arithmetic mean
#'   by elements occurring within assembly motif \cr
#' \code{\link{gmean_bymot_LOO}} geometric mean by assembly motif \cr
#' \code{\link{gmean_byelt_LOO}} geometric mean
#'   by elements occurring within assembly motif \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

amean_bymot_LOO <- function(fobs) {

  nbass <- length(fobs)

  if (nbass > 1) {
    fprd <- numeric(nbass)
    for (ind in seq_len(nbass)) fprd[ind] <- amean(fobs[-ind])
  } else {
    fprd <- NA
  }

  return(fprd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Predicting the performances
#'  by assembly motif
#'  using leave-one-out method
#'
#' @description Take a vector \code{fobs} of assembly performances
#' and return a vector of performances predicted
#' as the arithmetic mean of performances of other assemblages
#' that share the same assembly motif. \cr
#'
#' Assembly motifs are labelled in the vector \code{assMotif}.
#'
#' @usage validate_amean_bymot_LOO(fobs, assMotif)
#'
#' @inheritParams validate_using_cross_validation
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the arithmetic mean
#' of performances of assemblages
#' that share a same assembly motif,
#' by excluding the assemblage to predict.
#'
#' @details Predicted performances are computed
#' using arithmetic mean (\code{opt.mean = "amean"})
#' of performances of assemblages
#' that share a same assembly motif (\code{opt.model = "bymot"})
#' by excluding the assemblage to predict.
#'
#' @seealso
#' \code{\link{validate_amean_byelt_LOO}} arithmetic mean
#'   by elements occurring within assembly motif \cr
#' \code{\link{validate_gmean_bymot_LOO}} arithmetic mean
#'   by assembly motif \cr
#' \code{\link{validate_gmean_byelt_LOO}} geometric mean
#'   by elements occurring within assembly motif \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

validate_amean_bymot_LOO <- function(fobs, assMotif) {

  fprd <- numeric(length(assMotif))

  setMot <- unique(assMotif)
  for (mot in seq_along(setMot)) {

    indMot       <- which(assMotif == setMot[mot])
    fprd[indMot] <- amean_bymot_LOO(fobs[indMot])
  }

  return(fprd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Predicting the performances
#'  by assembly motif
#'  using leave-one-out method
#'  over several experiments
#'
#' @description Take a vector \code{fobs} of assembly performances
#' over several experiments
#' and return a vector of performances predicted
#' as the arithmetic mean of performances of other assemblages
#' that share the same assembly motif. \cr
#'
#' Assembly motifs are labelled in the vector \code{assMotif}.
#' Experiments are labelled in the vector \code{xpr}.
#'
#' @usage validate_amean_bymot_LOO_xpr(fobs, assMotif, xpr)
#'
#' @inheritParams validate_using_cross_validation
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the arithmetic mean of performances
#' of assemblages that share a same assembly motif,
#' experiment by experiment,
#' by excluding the assemblage to predict.
#'
#' @details Predicted performances are computed
#' using arithmetic mean (\code{opt.mean = "amean"})
#' of performances of assemblages
#' that share a same assembly motif (\code{opt.model = "bymot"}),
#' by excluding the assemblage to predict.
#'
#' @seealso
#' \code{\link{validate_amean_byelt_LOO_xpr}} arithmetic mean
#'   by elements occurring within assembly motif \cr
#' \code{\link{validate_gmean_bymot_LOO_xpr}} geometric mean
#'   by assembly motif \cr
#' \code{\link{validate_gmean_byelt_LOO_xpr}} geometric mean
#'   by elements occurring within assembly motif \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

validate_amean_bymot_LOO_xpr <- function(fobs, assMotif, xpr) {

  fprd <- numeric(length(assMotif))

  setXpr <- unique(names(xpr))
  for (ix in seq_along(setXpr)) {

    indXpr       <- which(names(xpr) == setXpr[ix])
    fprd[indXpr] <- validate_amean_bymot_LOO(fobs[indXpr], assMotif[indXpr])
  }

  return(fprd)
}


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Geometric mean by assembly motif                                         ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric mean
#' by assembly motif
#' using leave-one-out method
#'
#' @description Take a vector \code{fobs} of assembly performances
#' that share a same assembly motif,
#' and return a vector of performances
#' predicted as the geometric mean of performances
#' of all assemblages.
#'
#' @usage gmean_bymot_LOO(fobs)
#'
#' @inheritParams validate_using_cross_validation
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the geometric mean
#' of performances of assemblages
#' that share a same assembly motif,
#' by excluding the assemblage to predict.
#'
#' @details Predicted performances are computed
#' using geometric mean (\code{opt.mean = "gmean"})
#' of performances of assemblages
#' that share a same assembly motif (\code{opt.model = "bymot"}),
#' by excluding the assemblage to predict.
#'
#' @seealso
#' \code{\link{amean_bymot_LOO}} arithmetic mean
#'   by assembly motif \cr
#' \code{\link{amean_byelt_LOO}} arithmetic mean
#'   by elements occurring within assembly motif \cr
#' \code{\link{gmean_byelt_LOO}} geometric mean
#'   by elements occurring within assembly motif \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

gmean_bymot_LOO <- function(fobs) {

  nbass <- length(fobs)

  if (nbass > 1) {
    fprd <- numeric(nbass)
    for (ind in seq_len(nbass)) fprd[ind] <- gmean(fobs[-ind])
  } else {
    fprd <- NA
  }

  return(fprd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Predicting the performances
#' by assembly motif
#' using leave-one-out method
#'
#' @description Take a vector \code{fobs} of assembly performances
#' and return a vector of performances
#' predicted as the geometric mean
#' of performances of all assemblages
#' that share a same assembly motif,
#' by excluding the assemblage to predict.
#'
#' @usage validate_gmean_bymot_LOO(fobs, assMotif)
#'
#' @inheritParams validate_using_cross_validation
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the geometric mean of performances
#' of assemblages that share a same assembly motif.
#'
#' @details Predicted performances are computed
#' using geometric mean (\code{opt.mean = "gmean"})
#' of performances of assemblages
#' that share a same assembly motif (\code{opt.model = "bymot"})
#' by excluding the assemblage to predict.
#'
#' @seealso
#' \code{\link{validate_amean_bymot_LOO}} arithmetic mean
#'   by assembly motif \cr
#' \code{\link{validate_amean_byelt_LOO}} arithmetic mean
#'   by elements occurring within assembly motif \cr
#' \code{\link{validate_gmean_byelt_LOO}} geometric mean
#'   by elements occurring within assembly motif \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

validate_gmean_bymot_LOO <- function(fobs, assMotif) {

  fprd <- numeric(length(assMotif))

  setMot <- unique(assMotif)
  for (mot in seq_along(setMot)) {

    indMot       <- which(assMotif == setMot[mot])
    fprd[indMot] <- gmean_bymot_LOO(fobs[indMot])
  }

  return(fprd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Predicting the performances
#' by assembly motif
#' using leave-one-out method
#' over several experiments
#'
#' @description Take a vector \code{fobs} of assembly performances
#' over several experiments,
#' and return a vector of performances
#' predicted as the geometric mean
#' of performances of all assemblages
#' that share a same assembly motif,
#' by excluding the assemblage to predict.
#'
#' @usage validate_gmean_bymot_LOO_xpr(fobs, assMotif, xpr)
#'
#' @inheritParams validate_using_cross_validation
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the geometric mean
#' of performances of assemblages
#' that share a same assembly motif,
#' experiment by experiment.
#' by excluding the assemblage to predict.
#'
#' @details Predicted performances are computed
#' using geometric mean (\code{opt.mean = "amean"})
#' of performances of assemblages
#' that share a same assembly motif (\code{opt.model = "bymot"}),
#' experiment by experiment,
#' by excluding the assemblage to predict.
#'
#' @seealso
#' \code{\link{validate_amean_bymot_LOO_xpr}} arithmetic mean
#'   by assembly motif \cr
#' \code{\link{validate_amean_byelt_LOO_xpr}} arithmetic mean
#'   by elements occurring within assembly motif \cr
#' \code{\link{validate_gmean_byelt_LOO_xpr}} geometric mean
#'   by elements occurring within assembly motif \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

validate_gmean_bymot_LOO_xpr <- function(fobs, assMotif, xpr) {

  fprd <- numeric(length(assMotif))

  setXpr <- unique(names(xpr))
  for (ix in seq_along(setXpr)) {

    indXpr       <- which(names(xpr) == setXpr[ix])
    fprd[indXpr] <- validate_gmean_bymot_LOO(fobs[indXpr], assMotif[indXpr])
  }

  return(fprd)
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Arithmetic mean by element within each assembly motif                    ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Predicting the performances
#'  by elements occurring within assembly motif
#'  using leave-one-out method
#'
#' @description Take a vector \code{fobs}
#' of performances of assemblages
#' that share a same assembly motif,
#' and return a vector of performances
#' predicted as the arithmetic mean
#' of performances of assemblages
#' that contain the same elements as the assemblage to predict.\cr
#'
#' All assemblages share the same assembly motif.
#'
#' @usage amean_byelt_LOO(fobs, mOccur)
#'
#' @inheritParams validate_using_cross_validation
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the average
#' of mean performances of assemblages
#' that share the same assembly motif
#' and that contain the same elements
#' as the assemblage to predict \code{opt.model = "byelt"} .
#'
#' @details Modelled performances are computed
#' using arithmetic mean \code{opt.mean = "amean"} of performances.
#' Assemblages share a same assembly motif \code{opt.model = "bymot"}.
#' Modelled performances are the average
#' of mean performances of assemblages that contain the same elements
#' as the assemblage to predict,
#' except the assemblage to predict.
#' This procedure corresponds to a linear model with each assembly motif
#' based on the element occurrence in each assemblage.
#'
#' @seealso
#' \code{\link{amean_bymot_LOO}} arithmetic mean
#'   by assembly motif \cr
#' \code{\link{gmean_bymot_LOO}} geometric mean
#'   by assembly motif \cr
#' \code{\link{gmean_byelt_LOO}} geometric mean
#'   by elements occurring within assembly motif \cr
#'
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

amean_byelt_LOO <- function(fobs, mOccur) {

  nbass <- length(fobs)
  fprd  <- numeric(nbass)
  vfct  <- numeric(dim(mOccur)[2])

  for (ind in seq_len(nbass)) {

    vfct[] <- NA
    indOth <- seq_len(nbass)[-ind]
    setElt <- which(mOccur[ind, ] != 0)

    for (elt in seq_along(setElt)) {
      indElt <- which(mOccur[indOth, setElt[elt]] == 1)
      if (length(indElt) > 0)
        vfct[setElt[elt]] <- amean(fobs[indOth[indElt]])

      #      if (length(indElt) > 0) {
      #        vfct[setElt[elt]] <- amean(fobs[indOth[indElt]])
      #    } else {
      #        vfct[setElt[elt]] <- amean(fobs[indOth])
      #    }

    }

    fprd[ind] <- amean(vfct[setElt])
  }

  return(fprd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Predicting the performances
#'  by elements occurring within assembly motif
#'  using leave-one-out method
#'
#' @description Take a vector \code{fobs}
#' of performances of assemblages
#' that share a same assembly motif,
#' and return a vector of performances
#' predicted as the arithmetic mean
#' of performances of assemblages
#' that contain the same elements as the assemblage to predict.\cr
#'
#' All assemblages share the same assembly motif.
#'
#' @usage validate_amean_byelt_LOO(fobs, assMotif, mOccur)
#'
#' @inheritParams validate_using_cross_validation
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the average
#' of mean performances of assemblages
#' that share the same assembly motif
#' and that contain the same elements
#' as the assemblage to predict \code{opt.model = "byelt"} .
#'
#' @details Modelled performances are computed
#' using arithmetic mean \code{opt.mean = "amean"} of performances.
#' Assemblages share a same assembly motif \code{opt.model = "bymot"}.
#' Modelled performances are the average
#' of mean performances of assemblages that contain the same elements
#' as the assemblage to predict,
#' except the assemblage to predict.
#' This procedure corresponds to a linear model with each assembly motif
#' based on the element occurrence in each assemblage.
#'
#' @seealso
#' \code{\link{amean_bymot_LOO}} arithmetic mean
#'   by assembly motif \cr
#' \code{\link{gmean_bymot_LOO}} geometric mean
#'   by assembly motif \cr
#' \code{\link{gmean_byelt_LOO}} geometric mean
#'   by elements occurring within assembly motif \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

validate_amean_byelt_LOO <- function(fobs, assMotif, mOccur) {

  fprd <- numeric(length(assMotif))

  setMot <- unique(assMotif)
  for (mot in seq_along(setMot)) {
    indMot <- which(assMotif == setMot[mot])
    if (length(indMot) > 1) {
      fprd[indMot] <- amean_byelt_LOO(fobs[indMot], mOccur[indMot, ])
    } else {
      fprd[indMot] <- NA
    }
  }

  return(fprd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Predicting the performances
#'  by elements occurring within assembly motif
#'  using leave-one-out method
#'  over several experiments
#'
#' @description Take a vector \code{fobs} of assembly performances
#' over several experiments,
#' and return a vector of performances
#' predicted as the arithmetic mean
#' of performances of all assemblages
#' that share a same assembly motif,
#' by excluding the assemblage to predict.
#'
#' @usage validate_amean_byelt_LOO_xpr(fobs, assMotif, mOccur, xpr)
#'
#' @inheritParams validate_using_cross_validation
#'
#' @return Return a vector of same length as \code{fobs}, of
#' which values are computed as an arithmetic mean.
#'
#' @details Modelled performances are computed
#' using arithmetic mean (\code{opt.mean = "amean"}) of performances.
#' Assemblages share a same assembly motif (\code{opt.model = "bymot"}).
#' Modelled performances are the average
#' of mean performances of assemblages that contain the same elements
#' as the assemblage to predict,
#' experiment by experiment,
#' except the assemblage to predict.
#' This procedure corresponds to a linear model with each assembly motif
#' based on the element occurrence in each assemblage.
#'
#' @seealso
#' \code{\link{validate_amean_bymot_LOO_xpr}} arithmetic mean
#'   by assembly motif \cr
#' \code{\link{validate_gmean_bymot_LOO_xpr}} geometric mean
#'   by assembly motif \cr
#' \code{\link{validate_gmean_byelt_LOO_xpr}} geometric mean
#'   by elements occurring within assembly motif \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

validate_amean_byelt_LOO_xpr <- function(fobs, assMotif, mOccur, xpr) {

  fprd <- numeric(length(assMotif))

  setXpr <- unique(names(xpr))
  for (ix in seq_along(setXpr)) {

    indXpr         <- which(names(xpr) == setXpr[ix])
    fprd[indXpr] <- validate_amean_byelt_LOO(fobs[indXpr],
                                             assMotif[indXpr],
                                             mOccur[indXpr, ] )
  }

  return(fprd)
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Geometric mean by element within each assembly motif                     ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric mean
#'  by elements occurring within assembly motif
#'  using leave-one-out method
#'
#' @description Take a vector \code{fobs}
#' of performances of assemblages
#' that share a same assembly motif,
#' and return a vector of performances
#' predicted as the arithmetic mean
#' of performances of assemblages
#' that contain the same elements as the assemblage to predict.\cr
#'
#' All assemblages share the same assembly motif.
#'
#' @usage gmean_byelt_LOO(fobs, mOccur)
#'
#' @inheritParams validate_using_cross_validation
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the average
#' of mean performances of all assemblages
#' that contain the same elements
#' as the assemblage to predict (\code{opt.model = "byelt"}) .
#'
#' @details Predicted performances are computed
#' using geometric mean (\code{opt.mean = "gmean"})
#' of performances.
#' Assemblages share a same assembly motif (\code{opt.model = "bymot"}).
#' Predicted performances are the average
#' of mean performances of assemblages
#' that contain the same elements
#' as the assemblage to predict,
#' except the assemblage to predict.
#' This procedure corresponds to a linear model with each assembly motif
#' based on the element occurrence in each assemblage.
#'
#' @seealso
#' \code{\link{amean_bymot_LOO}} arithmetic mean
#'   by assembly motif \cr
#' \code{\link{amean_byelt_LOO}} arithmetic mean
#'   by elements occurring within assembly motif \cr
#' \code{\link{gmean_bymot_LOO}} geometric mean
#'   by assembly motif \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

gmean_byelt_LOO <- function(fobs, mOccur) {

  nbass  <- length(fobs)
  fprd <- numeric(nbass)
  vfct   <- numeric(dim(mOccur)[2])

  for (ind in seq_len(nbass)) {

    vfct[] <- NA
    indOth <- seq_len(nbass)[-ind]
    setElt <- which(mOccur[ind, ] != 0)

    for (elt in seq_along(setElt)) {
      indElt <- which(mOccur[indOth, setElt[elt]] == 1)
      if (length(indElt) > 0)
        vfct[setElt[elt]] <- gmean(fobs[indOth[indElt]])
    }

    fprd[ind] <- gmean(vfct[setElt])
  }

  return(fprd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric mean
#'  by elements occurring within assembly motif
#'  using leave-one-out method
#'
#' @description Take a vector \code{fobs}
#' of performances of assemblages
#' that share a same assembly motif,
#' and return a vector of performances
#' predicted as the geometric mean
#' of performances of assemblages
#' that contain the same elements as the assemblage to predict.\cr
#'
#' All assemblages share the same assembly motif.
#'
#' @usage validate_gmean_byelt_LOO(fobs, assMotif, mOccur)
#'
#' @inheritParams validate_using_cross_validation
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the average
#' of mean performances of all assemblages
#' that contain the same elements
#' as the assemblage to predict (\code{opt.model = "byelt"}) .
#'
#' @details Predicted performances are computed
#' using geometric mean (\code{opt.mean = "gmean"})
#' of performances.
#' Assemblages share a same assembly motif (\code{opt.model = "bymot"}).
#' Predicted performances are the average
#' of mean performances of assemblages
#' that contain the same elements
#' as the assemblage to predict,
#' except the assemblage to predict.
#' This procedure corresponds to a linear model with each assembly motif
#' based on the element occurrence in each assemblage.
#'
#' @seealso
#' \code{\link{validate_amean_bymot_LOO}} arithmetic mean
#'   by assembly motif \cr
#' \code{\link{validate_amean_byelt_LOO}} arithmetic mean
#'   by elements occurring within assembly motif \cr
#' \code{\link{validate_gmean_bymot_LOO}} geometric mean
#'   by assembly motif \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

validate_gmean_byelt_LOO <- function(fobs, assMotif, mOccur) {

  fprd <- numeric(length(assMotif))

  setMot <- unique(assMotif)
  for (mot in seq_along(setMot)) {

    indMot <- which(assMotif == setMot[mot])
    if (length(indMot) > 1) {
      fprd[indMot] <- gmean_byelt_LOO(fobs[indMot], mOccur[indMot, ])
    } else {
      fprd[indMot] <- NA
    }
  }

  return(fprd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric mean
#'  by elements occurring within assembly motif
#'  using leave-one-out method
#'  over several experiments
#'
#' @description Take a vector \code{fobs}
#' of performances of assemblages
#' over several experiments
#' that share a same assembly motif,
#' and return a vector of performances
#' predicted as the geometric mean
#' of performances of assemblages
#' that contain the same elements as the assemblage to predict.\cr
#'
#' All assemblages share the same assembly motif.
#'
#' @usage validate_gmean_byelt_LOO_xpr(fobs, assMotif, mOccur, xpr)
#'
#' @inheritParams validate_using_cross_validation
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the average
#' of mean performances of all assemblages,
#' experiment by experiment,
#' that contain the same elements
#' as the assemblage to predict (\code{opt.model = "byelt"}) .
#'
#' @details Predicted performances are computed
#' using geometric mean (\code{opt.mean = "gmean"})
#' of performances.
#' Assemblages share a same assembly motif (\code{opt.model = "bymot"}).
#' Predicted performances are the average,
#' experiment by experiment,
#' of mean performances of assemblages
#' that contain the same elements
#' as the assemblage to predict,
#' except the assemblage to predict.
#' This procedure corresponds to a linear model with each assembly motif
#' based on the element occurrence in each assemblage.
#'
#' @seealso
#'  \code{\link{validate_amean_bymot_LOO_xpr}} arithmetic mean
#'   byassembly motif \cr
#'  \code{\link{validate_amean_byelt_LOO_xpr}} arithmetic mean
#'   by elements occurring within assembly motif \cr
#'  \code{\link{validate_gmean_bymot_LOO_xpr}} geometric mean
#'   by assembly motif \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

validate_gmean_byelt_LOO_xpr <- function(fobs, assMotif, mOccur, xpr) {

  fprd    <- numeric(length(assMotif))

  setXpr    <- unique(names(xpr))
  for (ix in seq_along(setXpr)) {

    index         <- which(names(xpr) == setXpr[ix])
    fprd[index] <- validate_gmean_byelt_LOO(fobs[index],
                                            assMotif[index],
                                            mOccur[index, ] )
  }

  return(fprd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# End of file                                                              ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
