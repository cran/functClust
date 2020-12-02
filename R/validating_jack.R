#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                                PREDICTING_JCK.R
#
# Functions for predicting by cross-validation (jackknife)                 ####
#      the performances of assemblages
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#' @include
#'  stats.R
#'  validating_loo.R
#'
NULL


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Arithmetic mean by assembly motif                                        ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic mean of performances
#' by assembly motif
#' using jackknife method
#'
#' @description Take a vector \code{fobs}
#' of performances of assemblages
#' that share the same assembly motif,
#' and return a vector of performances predicted
#' as the arithmetic mean
#' of performances of other assemblages.
#'
#' @usage amean_bymot_jack(fobs, jack)
#'
#' @inheritParams validate_using_cross_validation
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the arithmetic mean
#' of performances of assemblages
#' that share a same assembly motif,
#' by excluding a subset of assemblages
#' containing the assemblage to predict.
#'
#' @details Predicted performances are computed using arithmetic mean
#'   \code{opt.mean = "amean"} of performances of assemblages
#'   that share a same assembly motif \code{opt.model = "bymot"}.
#'
#' The assemblages belonging to a same assembly motif are divided
#' into \code{jack[2]} subsets of \code{jack[1]} assemblages.
#' Prediction is computed by excluding \code{jack[1]} assemblages,
#' including the assemblage to predict.
#' If the total number of assemblages belonging
#' to the assembly motif is lower than \code{jack[1]*jack[2]},
#' prediction is computed by Leave-One-Out (LOO).
#'
#' @seealso \code{\link{amean_byelt_jack}},
#'          \code{\link{gmean_bymot_jack}},
#'          \code{\link{gmean_byelt_jack}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

amean_bymot_jack <- function(fobs, jack) {

  nbass <- length(fobs)

  if (nbass > jack[1] * jack[2]) {

    fprd  <- numeric(nbass)
    index <- sample.int(nbass)
    size  <- floor(nbass / jack[2])

    for (ind in seq_len(jack[2] - 1)) {
      indjack         <- index[(ind - 1) * size + (1:size)]
      fprd[indjack] <- amean(fobs[-indjack])
    }

    indjack       <- index[(ind * size + 1):nbass]
    fprd[indjack] <- amean(fobs[-indjack])

  } else {

    fprd <- amean_bymot_LOO(fobs)
  }

  return(fprd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Predicting the performances
#'  by assembly motif
#'  using jackknife method
#'
#' @description Take a vector \code{fobs} of assembly performances
#' and return a vector of performances predicted
#' as the arithmetic mean of performances of other assemblages
#' that share the same assembly motif. \cr
#'
#' Assembly motifs are labelled in the vector \code{assMotif}.
#'
#' @usage validate_amean_bymot_jack(fobs, assMotif, jack)
#'
#' @inheritParams validate_using_cross_validation
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the arithmetic mean
#' of performances of assemblages
#' that share a same assembly motif,
#' by excluding a subset of assemblages
#' containing the assemblage to predict.
#'
#' @details Predicted performances are computed
#' using arithmetic mean (\code{opt.mean = "amean"})
#' of performances of assemblages
#'   that share a same assembly motif (\code{opt.model = "bymot"}).
#'
#' The assemblages belonging to a same assembly motif are divided
#'  into \code{jack[2]} subsets of \code{jack[1]} assemblages.
#'   Prediction is computed by excluding \code{jack[1]} assemblages,
#'   including the assemblage to predict.
#'    If the total number of assemblages belonging
#'   to the assembly motif is lower than \code{jack[1]*jack[2]},
#'    prediction is computed by Leave-One-Out (LOO).
#'
#' @seealso \code{\link{validate_amean_byelt_jack_xpr}},
#'          \code{\link{validate_gmean_bymot_jack_xpr}},
#'          \code{\link{validate_gmean_byelt_jack_xpr}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

validate_amean_bymot_jack <- function(fobs, assMotif, jack) {

  fprd <- numeric(length(assMotif))

  setMot <- unique(assMotif)
  for (mot in seq_along(setMot)) {

    indMot       <- which(assMotif == setMot[mot])
    fprd[indMot] <- amean_bymot_jack(fobs[indMot], jack)
  }

  return(fprd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Predicting the performances
#'  by assembly motif
#'  using jackknife method
#'  over several experiments
#'
#' @description Take a vector \code{fobs} of assembly performances
#' and return a vector of performances predicted
#' as the arithmetic mean of performances of other assemblages
#' that share the same assembly motif. \cr
#'
#' Assembly motifs are labelled in the vector \code{assMotif}.
#' Experiments are labelled in the vector \code{xpr}.
#'
#' @usage validate_amean_bymot_jack_xpr(fobs, assMotif, jack, xpr)
#'
#' @inheritParams validate_using_cross_validation
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the arithmetic mean
#' of performances of assemblages
#' that share a same assembly motif,
#' experiment by experiment,
#' by excluding a subset of assemblages
#' containing the assemblage to predict.
#'
#' @details Predicted performances are computed
#' using arithmetic mean (\code{opt.mean = "amean"})
#' of performances of assemblages
#' that share a same assembly motif (\code{opt.model = "bymot"}).
#'
#' The assemblages belonging to a same assembly motif are divided
#'  into \code{jack[2]} subsets of \code{jack[1]} assemblages.
#'   Prediction is computed by excluding \code{jack[1]} assemblages,
#'   including the assemblage to predict.
#'    If the total number of assemblages belonging
#'   to the assembly motif is lower than \code{jack[1]*jack[2]},
#'    prediction is computed by Leave-One-Out (LOO).
#'
#' @seealso \code{\link{validate_amean_byelt_jack_xpr}},
#'          \code{\link{validate_gmean_bymot_jack_xpr}},
#'          \code{\link{validate_gmean_byelt_jack_xpr}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

validate_amean_bymot_jack_xpr <- function(fobs, assMotif, jack, xpr) {

  fprd <- numeric(length(assMotif))

  setXpr <- unique(names(xpr))
  for (ix in seq_along(setXpr)) {

    indXpr       <- which(names(xpr) == setXpr[ix])
    fprd[indXpr] <- validate_amean_bymot_jack(fobs[indXpr],
                                              assMotif[indXpr], jack)
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
#' using jackknife method
#'
#' @description Take a vector \code{fobs}
#' of performances of assemblages
#' that share the same assembly motif,
#' and return a vector of performances predicted
#' as the geometric mean
#' of performances of other assemblages.
#'
#' @usage gmean_bymot_jack(fobs, jack)
#'
#' @inheritParams validate_using_cross_validation
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the geometric mean of performances
#' of assemblages that share a same assembly motif,
#' by excluding a subset of assemblages
#' containing the assemblage to predict.
#'
#' @details Predicted performances are computed
#' using geometric mean (\code{opt.mean = "gmean"})
#' of performances of assemblages
#' that share a same assembly motif (\code{opt.model = "bymot"}).
#'
#' The assemblages belonging to a same assembly motif are divided
#'  into \code{jack[2]} subsets of \code{jack[1]} assemblages.
#'   Prediction is computed by excluding \code{jack[1]} assemblages,
#'   including the assemblage to predict.
#'    If the total number of assemblages belonging
#'   to the assembly motif is lower than \code{jack[1]*jack[2]},
#'    prediction is computed by Leave-One-Out (LOO).
#'
#' @seealso \code{\link{amean_bymot_jack}},
#'          \code{\link{amean_byelt_jack}},
#'          \code{\link{gmean_byelt_jack}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

gmean_bymot_jack <- function(fobs, jack) {

  nbass <- length(fobs)

  if (nbass > jack[1] * jack[2]) {

    fprd <- numeric(nbass)
    index  <- sample.int(nbass)
    size   <- floor(nbass / jack[2])

    for (ind in seq_len(jack[2] - 1)) {
      indjack         <- index[(ind - 1) * size + (1:size)]
      fprd[indjack] <- gmean(fobs[-indjack])
    }

    indjack         <- index[(ind * size + 1):nbass]
    fprd[indjack] <- gmean(fobs[-indjack])

  } else {

    fprd <- gmean_bymot_LOO(fobs)
  }

  return(fprd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Predicting the performances
#' by assembly motif
#' using jackknife method
#'
#' @description Take a vector \code{fobs}
#' of performances of assemblages
#' and return a vector of performances
#' predicted as the geometric mean
#' of performances of other assemblages.
#'
#' @usage validate_gmean_bymot_jack(fobs, assMotif, jack)
#'
#' @inheritParams validate_using_cross_validation
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the geometric mean of performances
#' of assemblages that share a same assembly motif,
#' by excluding a subset of assemblages
#' containing the assemblage to predict.
#'
#' @details Predicted performances are computed
#' using geometric mean (\code{opt.mean = "gmean"})
#'  of performances of assemblages
#'   that share a same assembly motif (\code{opt.model = "bymot"}).
#'
#' The assemblages belonging to a same assembly motif are divided
#'  into \code{jack[2]} subsets of \code{jack[1]} assemblages.
#'   Prediction is computed by excluding \code{jack[1]} assemblages,
#'   including the assemblage to predict.
#'    If the total number of assemblages belonging
#'   to the assembly motif is lower than \code{jack[1]*jack[2]},
#'    prediction is computed by Leave-One-Out (LOO).
#'
#' @seealso \code{\link{validate_amean_bymot_jack}},
#'          \code{\link{validate_amean_byelt_jack}},
#'          \code{\link{validate_gmean_byelt_jack}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

validate_gmean_bymot_jack <- function(fobs, assMotif, jack) {

  fprd <- numeric(length(assMotif))

  setMot <- unique(assMotif)
  for (mot in seq_along(setMot)) {

    indMot         <- which(assMotif == setMot[mot])
    fprd[indMot] <- gmean_bymot_jack(fobs[indMot], jack)
  }

  return(fprd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Predicting the performances
#' by assembly motif
#' using jackknife method
#' over several experiments
#'
#' @description Take a vector \code{fobs}
#' of performances of assemblages
#' over several experiments
#' and return a vector of performances
#' predicted as the geometric mean,
#' experiment by experiment,
#' of performances of other assemblages.
#'
#' @usage validate_gmean_bymot_jack_xpr(fobs, assMotif, jack, xpr)
#'
#' @inheritParams validate_using_cross_validation
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the geometric mean of performances
#' of assemblages that share a same assembly motif,
#' experiment by experiment,
#' by excluding a subset of assemblages
#' containing the assemblage to predict.
#'
#' @details Predicted performances are computed
#' using geometric mean (\code{opt.mean = "gmean"})
#'  of performances of assemblages
#'   that share a same assembly motif (\code{opt.model = "bymot"}).
#'
#' The assemblages belonging to a same assembly motif are divided
#'  into \code{jack[2]} subsets of \code{jack[1]} assemblages.
#'   Prediction is computed by excluding \code{jack[1]} assemblages,
#'   including the assemblage to predict.
#'    If the total number of assemblages belonging
#'   to the assembly motif is lower than \code{jack[1]*jack[2]},
#'    prediction is computed by Leave-One-Out (LOO).
#'
#' @seealso \code{\link{validate_amean_bymot_jack_xpr}},
#'          \code{\link{validate_amean_byelt_jack_xpr}},
#'          \code{\link{validate_gmean_byelt_jack_xpr}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

validate_gmean_bymot_jack_xpr <- function(fobs, assMotif, jack, xpr) {

  fprd <- numeric(length(assMotif))

  setXpr <- unique(names(xpr))
  for (ix in seq_along(setXpr)) {

    indXpr         <- which(names(xpr) == setXpr[ix])
    fprd[indXpr] <- validate_gmean_bymot_jack(assMotif[indXpr],
                                              fobs[indXpr], jack)
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
#' @title Arithmetic mean
#'  by elements occurring within assembly motif
#'  using jackknife method
#'
#' @description Take a numeric vector and return the predicted vector
#' computed as the arithmetic mean
#' of all elements belonging to a same assembly motif.
#'
#' @usage amean_byelt_jack(fobs, mOccur, jack)
#'
#' @inheritParams validate_using_cross_validation
#'
#' @details Modelled performances are computed
#' using arithmetic mean (\code{opt.mean = "amean"}) of performances.
#' Assemblages share a same assembly motif (\code{opt.model = "bymot"}).
#' Modelled performances are the average
#' of mean performances of assemblages that contain the same elements
#' as the assemblage to predict,
#' except a subset of assemblages.
#' This procedure corresponds to a linear model with each assembly motif
#' based on the element occurrence in each assemblage.
#'
#' The assemblages belonging to a same assembly motif are divided
#'  into \code{jack[2]} subsets of \code{jack[1]} assemblages.
#'   Prediction is computed by excluding \code{jack[1]} assemblages,
#'   including the assemblage to predict.
#'    If the total number of assemblages belonging
#'   to the assembly motif is lower than \code{jack[1]*jack[2]},
#'    prediction is computed by leave-one-out (LOO).
#'
#' @seealso \code{\link{amean_bymot_jack}},
#'          \code{\link{gmean_bymot_jack}},
#'          \code{\link{gmean_byelt_jack}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

amean_byelt_jack <- function(fobs, mOccur, jack) {

  nbass  <- length(fobs)
  fprd <- numeric(nbass)

  if (nbass > jack[1] * jack[2]) {

    index <- sample.int(nbass)
    size  <- floor(nbass / jack[2])

    for (ind in seq_len(jack[2] - 1)) {

      indjack <- index[(ind - 1) * size + (1:size)]
      indOth  <- seq_len(nbass)[-indjack]

      tmp             <- mOccur[indOth, ] * fobs[indOth]
      tmp[tmp == 0]   <- NA
      vfct            <- apply(tmp, MARGIN = 2, FUN = amean)

      tmp             <- t(t(mOccur[indjack, ]) * vfct)
      tmp[tmp == 0]   <- NA
      fprd[indjack] <- apply(tmp, MARGIN = 1, FUN = amean)
    }

    indjack <- index[(ind * size + 1):nbass]
    indOth  <- seq_len(nbass)[-indjack]

    tmp             <- mOccur[indOth, ] * fobs[indOth]
    tmp[tmp == 0]   <- NA
    vfct            <- apply(tmp, MARGIN = 2, FUN = amean)

    tmp             <- t(t(mOccur[indjack, ]) * vfct)
    tmp[tmp == 0]   <- NA
    fprd[indjack] <- apply(tmp, MARGIN = 1, FUN = amean)

  } else {

    fprd[ ] <- amean_byelt_LOO(fobs, mOccur)
  }

  return(fprd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Predicting the performances
#'  by elements occurring within assembly motif
#'  using jackknife method
#'
#' @description Take a numeric vector and return the predicted vector
#'  computed as the arithmetic mean of all elements belonging to the same motif.
#'
#' @usage validate_amean_byelt_jack(fobs, assMotif, mOccur, jack)
#'
#' @inheritParams validate_using_cross_validation
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the arithmetic mean
#' of performances of assemblages
#' that share a same assembly motif,
#' by excluding a subset of assemblages
#' containing the assemblage to predict.
#'
#'
#' @details Predicted performances are computed
#' using arithmetic mean (\code{opt.mean = "amean"})
#'  of performances of assemblages
#'   that share a same assembly motif (\code{opt.model = "bymot"}).
#'
#' The assemblages belonging to a same assembly motif are divided
#'  into \code{jack[2]} subsets of \code{jack[1]} assemblages.
#'   Prediction is computed by excluding \code{jack[1]} assemblages,
#'   including the assemblage to predict.
#'    If the total number of assemblages belonging
#'   to the assembly motif is lower than \code{jack[1]*jack[2]},
#'    prediction is computed by Leave-One-Out (LOO).
#'
#' @seealso \code{\link{validate_amean_bymot_jack}},
#'          \code{\link{validate_gmean_bymot_jack}},
#'          \code{\link{validate_gmean_byelt_jack}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

validate_amean_byelt_jack <- function(fobs, assMotif, mOccur, jack) {

  fprd <- numeric(length(assMotif))

  setMot <- unique(assMotif)
  for (mot in seq_along(setMot)) {

    indMot <- which(assMotif == setMot[mot])
    if (length(indMot) > 1) {
      fprd[indMot] <- amean_byelt_jack(fobs[indMot], mOccur[indMot, ], jack)
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
#'  using jackknife method
#'  over several experiments
#'
#' @description Take a numeric vector and return the predicted vector
#' computed as the arithmetic mean of all elements belonging to a same motif.
#'
#' @usage validate_amean_byelt_jack_xpr(fobs, assMotif, mOccur, jack, xpr)
#'
#' @inheritParams validate_using_cross_validation
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the arithmetic mean
#' of performances of assemblages
#' that share a same assembly motif,
#' by excluding a subset of assemblages
#' containing the assemblage to predict.
#'
#' @details Predicted performances are computed
#' using arithmetic mean (\code{opt.mean = "amean"})
#'  of performances of assemblages
#'   that share a same assembly motif (\code{opt.model = "bymot"}).
#'
#' The assemblages belonging to a same assembly motif are divided
#'  into \code{jack[2]} subsets of \code{jack[1]} assemblages.
#'   Prediction is computed by excluding \code{jack[1]} assemblages,
#'   including the assemblage to predict.
#'    If the total number of assemblages belonging
#'   to the assembly motif is lower than \code{jack[1]*jack[2]},
#'    prediction is computed by Leave-One-Out (LOO).
#'
#' @seealso \code{\link{validate_amean_bymot_jack_xpr}},
#'          \code{\link{validate_gmean_bymot_jack_xpr}},
#'          \code{\link{validate_gmean_byelt_jack_xpr}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

validate_amean_byelt_jack_xpr <- function(fobs, assMotif, mOccur, jack, xpr) {

  fprd <- numeric(length(assMotif))

  setXpr <- unique(names(xpr))
  for (ix in seq_along(setXpr)) {

    index         <- which(names(xpr) == setXpr[ix])
    fprd[index] <- validate_amean_byelt_jack(fobs[index],
                                             assMotif[index],
                                             mOccur[index, ],
                                             jack  )
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
#'  using jackknife method
#'
#' @description Take a numeric vector and return the predicted vector
#' computed as the arithmetic mean of all elements belonging to the same motif.
#'
#' @usage gmean_byelt_jack(fobs, mOccur, jack)
#'
#' @inheritParams validate_using_cross_validation
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the geometric mean
#' of performances of assemblages,
#' by excluding a subset of assemblages
#' containing the assemblage to predict.
#'
#' @details Predicted performances are computed
#'  using geometric mean (\code{opt.mean = "gmean"})
#'   of performances of assemblages
#'   that share a same assembly motif (\code{opt.model = "bymot"}).
#'
#' The assemblages belonging to a same assembly motif are divided
#'  into \code{jack[2]} subsets of \code{jack[1]} assemblages.
#'   Prediction is computed by excluding \code{jack[1]} assemblages,
#'   including the assemblage to predict.
#'    If the total number of assemblages belonging
#'   to the assembly motif is lower than \code{jack[1]*jack[2]},
#'    prediction is computed by Leave-One-Out (LOO).
#'
#' @seealso \code{\link{amean_bymot_jack}},
#'          \code{\link{amean_byelt_jack}},
#'          \code{\link{gmean_bymot_jack}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

gmean_byelt_jack <- function(fobs, mOccur, jack) {

  nbass  <- length(fobs)
  fprd <- numeric(nbass)

  if (nbass > jack[1] * jack[2]) {

    index <- sample.int(nbass)
    size  <- floor(nbass / jack[2])

    for (ind in seq_len(jack[2] - 1)) {

      indjack <- index[(ind - 1) * size + (1:size)]
      indOth  <- seq_len(nbass)[-indjack]

      tmp             <- mOccur[indOth, ] * fobs[indOth]
      tmp[tmp == 0]   <- NA
      vfct            <- apply(tmp, MARGIN = 2, FUN = gmean)

      tmp             <- t(t(mOccur[indjack, ]) * vfct)
      tmp[tmp == 0]   <- NA
      fprd[indjack] <- apply(tmp, MARGIN = 1, FUN = gmean)
    }

    indjack <- index[(ind * size + 1):nbass]
    indOth  <- seq_len(nbass)[-indjack]

    tmp             <- mOccur[indOth, ] * fobs[indOth]
    tmp[tmp == 0]   <- NA
    vfct            <- apply(tmp, MARGIN = 2, FUN = gmean)

    tmp             <- t(t(mOccur[indjack, ]) * vfct)
    tmp[tmp == 0]   <- NA
    fprd[indjack] <- apply(tmp, MARGIN = 1, FUN = gmean)

  } else {

    fprd[] <- gmean_byelt_LOO(fobs, mOccur)
  }

  return(fprd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Predicting the performances
#'  by elements occurring within assembly motif
#'  using jackknife method
#'
#' @description Take a numeric vector and return the predicted vector
#'  computed as the arithmetic mean of all elements belonging to the same motif.
#'
#' @usage validate_gmean_byelt_jack(fobs, assMotif, mOccur, jack)
#'
#' @inheritParams validate_using_cross_validation
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the geometric mean
#' of performances of assemblages
#' that share a same assembly motif,
#' by excluding a subset of assemblages
#' containing the assemblage to predict.
#'
#' @details Predicted performances are computed
#' using geometric mean (\code{opt.mean = "gmean"})
#'  of performances of assemblages
#'   that share a same assembly motif (\code{opt.model = "bymot"}).
#'
#' The assemblages belonging to a same assembly motif are divided
#'  into \code{jack[2]} subsets of \code{jack[1]} assemblages.
#'   Prediction is computed by excluding \code{jack[1]} assemblages,
#'   including the assemblage to predict.
#'    If the total number of assemblages belonging
#'   to the assembly motif is lower than \code{jack[1]*jack[2]},
#'    prediction is computed by Leave-One-Out (LOO).
#'
#' @seealso \code{\link{validate_amean_bymot_jack}},
#'          \code{\link{validate_amean_byelt_jack}},
#'          \code{\link{validate_gmean_bymot_jack}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

validate_gmean_byelt_jack <- function(fobs, assMotif, mOccur, jack) {

  fprd <- numeric(length(assMotif))

  setMot <- unique(assMotif)
  for (mot in seq_along(setMot)) {

    indMot <- which(assMotif == setMot[mot])
    if (length(indMot) > 1) {
      fprd[indMot] <- gmean_byelt_jack(fobs[indMot], mOccur[indMot, ], jack)
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
#'  using jackknife method
#'  over several experiments
#'
#' @description Take a vector \code{fobs} of assembly performances
#' and return a vector of performances predicted
#' as the geometric mean of performances of other assemblages
#' that share the same assembly motif.
#' Assembly motifs are described in the vector \code{assMotif}.
#' Experiments are specified in the vector \code{xpr}.
#'
#' @usage validate_gmean_byelt_jack_xpr(fobs, assMotif, mOccur, jack, xpr)
#'
#' @inheritParams validate_using_cross_validation
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the geometric mean of performances
#' of assemblages that share a same assembly motif,
#' by excluding a subset of assemblages
#' containing the assemblage to predict.
#'
#' @details Predicted performances are computed
#' using geometric mean (\code{opt.mean = "gmean"})
#' of performances of assemblages
#'   that share a same assembly motif (\code{opt.model = "bymot"}).
#'
#' The assemblages belonging to a same assembly motif are divided
#'  into \code{jack[2]} subsets of \code{jack[1]} assemblages.
#'   Prediction is computed by excluding \code{jack[1]} assemblages,
#'   including the assemblage to predict.
#'    If the total number of assemblages belonging
#'   to the assembly motif is lower than \code{jack[1]*jack[2]},
#'    prediction is computed by Leave-One-Out (LOO).
#'
#' @seealso \code{\link{validate_amean_bymot_jack_xpr}},
#'          \code{\link{validate_amean_byelt_jack_xpr}},
#'          \code{\link{validate_gmean_bymot_jack_xpr}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

validate_gmean_byelt_jack_xpr <- function(fobs, assMotif, mOccur, jack, xpr) {

  fprd <- numeric(length(assMotif))

  setXpr <- unique(names(xpr))
  for (ix in seq_along(setXpr)) {

    index         <- which(names(xpr) == setXpr[ix])
    fprd[index] <- validate_gmean_byelt_jack(fobs[index],
                                             assMotif[index],
                                             mOccur[index, ],
                                             jack  )
  }

  return(fprd)
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Function for switch on different options                                 ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Predicting by cross-validation of assembly performances
#'
#' @description Take a vector \code{fobs} of assembly performances
#' over several experiments
#' and return a vector of performances
#' predicted as the mean performances of assemblages
#' that share the same assembly motif.\cr

#' Assembly motifs are labelled in the vector \code{assMotif}.
#' Experiments are labelled in the vector \code{xpr}.
#' Modelling options are indicated in \code{opt.mean} and \code{opt.model}.
#' Occurrence matrix \code{mOccur} is used if \code{opt.model =  "byelt"}.
#' Cross-validation is leave-one-out or jackknifesi
#'
#' @usage
#' validate_using_cross_validation(fobs, assMotif, mOccur, xpr,
#'                   opt.mean = "amean", opt.model = "bymot",
#'                   opt.jack = FALSE, jack = c(3,4)  )
#'
#' @param fobs a numeric vector. The vector \code{fobs} contains the
#' quantitative performances of assemblages.
#'
#' @param assMotif a vector of labels of \code{length(fobs)}.
#' The vector \code{assMotif} contains the assembly motifs of assemblages.
#'
#' @param mOccur a matrix of occurrence (occurrence of elements).
#' Its first dimension equals to \code{length(fobs)}. Its second dimension
#'  equals to the number of elements.
#'
#' @param xpr a vector of numerics of \code{length(fobs)}.
#' The vector \code{xpr} contains the weight of each experiment,
#' and the labels (in \code{names(xpr)}) of different experiments.
#' The weigth of each experiment is used
#' in the computation of the Residual Sum of Squares
#' in the function \code{rss_clustering}.
#' The used formula is \code{rss}
#' if each experiment has the same weight.
#' The used formula is \code{wrss}
#' (barycenter of RSS for each experiment)
#' if each experiment has different weights.
#' All assemblages that belong to a given experiment
#' should then have a same weigth.
#' Each experiment is identified by its names (\code{names(xpr)})
#' and the RSS of each experiment is weighted by values of \code{xpr}.
#' The vector \code{} is generated
#' by the function \code{stats::setNames}.
#'
#' @param opt.mean switchs to arithmetic formula \code{opt.mean = "amean"}
#' or geometric formula \code{opt.mean = "gmean"}.
#'
#' @param opt.model  switchs to model type:
#'  simple mean by assembly motif \code{opt.model = "bymot"}
#'  or linear model with assembly motif \code{opt.model = "byelt"}.
#'
#' @param opt.jack a logical,
#' that switchs towards cross-validation method.
#'
#' If \code{opt.jack = FALSE}, a "leave-one-out" is used:
#' predicted performances are computed
#' as the mean of performances of assemblages
#' that share a same assembly motif,
#' experiment by experiment,
#' except the only assemblage to predict. \cr
#'
#' If \code{opt.jack = TRUE}, a jackknife method is used:
#' the set of assemblages belonging to a same assembly motif is divided
#' into \code{jack[2]} subsets of \code{jack[1]} assemblages.
#' Predicted performances are computed,
#' experiment by experiment,
#' by excluding \code{jack[1]} assemblages,
#' including the assemblage to predict.
#' If the total number of assemblages belonging
#' to the assembly motif is lower than \code{jack[1]*jack[2]},
#' predictions are computed by Leave-One-Out method.
#'
#' @param jack an integer vector of length \code{2}.
#' The vector specifies the parameters for jackknife method.
#' The first integer \code{jack[1]} specifies the size of subset,
#' the second integer \code{jack[2]} specifies the number of subsets.
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are predicted
#' according to \code{opt.mean} and \code{opt.model}.
#'
#' @details None.
#'
#' @seealso
#' \code{\link{calibrate_byminrss}} \cr
#' \code{\link{predict_performance}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

validate_using_cross_validation <- function(fobs, assMotif, mOccur, xpr,
                              opt.mean = "amean",
                              opt.model = "bymot",
                              opt.jack = FALSE,   jack = c(3, 4)) {

  optmean <- "amean"
  if (opt.mean == "gmean") optmean <- "gmean"

  optmod <- "bymot"
  if (opt.model == "byelt") optmod <- "byelt"

  optjack <- ""
  if (opt.jack == TRUE) optjack <- "jack"

  optxpr <- ""
  if (length(unique(names(xpr))) != 1) optxpr <- "xpr"

  option <- paste(optmean, optmod, optjack, optxpr, sep = ".")

  return(
    switch(option,
           amean.bymot.. =
             validate_amean_bymot_LOO(fobs, assMotif),
           amean.bymot..xpr =
             validate_amean_bymot_LOO_xpr(fobs, assMotif, xpr),
           amean.bymot.jack. =
             validate_amean_bymot_jack(fobs, assMotif, jack),
           amean.bymot.jack.xpr =
             validate_amean_bymot_jack_xpr(fobs, assMotif, jack, xpr),

           gmean.bymot.. =
             validate_gmean_bymot_LOO(fobs, assMotif),
           gmean.bymot..xpr =
             validate_gmean_bymot_LOO_xpr(fobs, assMotif, xpr),
           gmean.bymot.jack. =
             validate_gmean_bymot_jack(fobs, assMotif, jack),
           gmean.bymot.jack.xpr =
             validate_gmean_bymot_jack_xpr(fobs, assMotif, jack, xpr),

           amean.byelt.. =
             validate_amean_byelt_LOO(fobs, assMotif, mOccur),
           amean.byelt..xpr =
             validate_amean_byelt_LOO_xpr(fobs, assMotif, mOccur, xpr),
           amean.byelt.jack. =
             validate_amean_byelt_jack(fobs, assMotif, mOccur, jack),
           amean.byelt.jack.xpr =
             validate_amean_byelt_jack_xpr(fobs, assMotif, mOccur, jack, xpr),

           gmean.byelt.. =
             validate_gmean_byelt_LOO(fobs, assMotif, mOccur),
           gmean.byelt..xpr =
             validate_gmean_byelt_LOO_xpr(fobs, assMotif, mOccur, xpr),
           gmean.byelt.jack. =
             validate_gmean_byelt_jack(fobs, assMotif, mOccur, jack),
           gmean.byelt.jack.xpr =
             validate_gmean_byelt_jack_xpr(fobs, assMotif, mOccur, jack, xpr) )
  )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# End of file                                                              ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
