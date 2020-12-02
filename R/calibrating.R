#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                             CALIBRATING.R
#
# Functions for modelling the performances of assemblages                  ####
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
#' @title Modelling of performances
#'  by assembly motif
#'
#' @description Take a vector \code{fobs} of assembly performances
#' and return a vector of performances
#' predicted as the arithmetic mean
#' of performances of all assemblages
#' that share a same assembly motif. \cr
#'
#' Assembly motifs are labelled in the vector \code{assMotif}.
#'
#' @usage calibrate_amean_bymot(fobs, assMotif)
#'
#' @inheritParams calibrate_byminrss
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the arithmetic mean
#' of performances of all assemblages
#' that share a same assembly motif.
#'
#' @details Modelled performances are computed
#' using arithmetic mean (\code{opt.mean = "amean"})
#' of performances of assemblages
#' that share a same assembly motif (\code{opt.model = "bymot"})
#' by including all assemblages that belong to a same assembly motif.
#'
#' @seealso
#' \code{\link{calibrate_amean_bymot}} arithmetic mean
#'   by assembly motif. \cr
#' \code{\link{calibrate_amean_byelt}} arithmetic mean
#'   by elements occurring within assembly motif. \cr
#' \code{\link{calibrate_gmean_bymot}} geometric mean
#'   by assembly motif. \cr
#' \code{\link{calibrate_gmean_byelt}} geometric mean
#'   by elements occurring within assembly motif. \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

calibrate_amean_bymot <- function(fobs, assMotif) {

  fprd <- numeric(length(assMotif))

  setMot <- unique(assMotif)
  for (mot in seq_along(setMot)) {

    indMot       <- which(assMotif == setMot[mot])
    fprd[indMot] <- amean(fobs[indMot])
  }

  return(fprd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Modelling of performances
#'  by assembly motif
#'  over several experiments
#'
#' @description Take a vector \code{fobs} of assembly performances
#' over several experiments
#' and return a vector of performances
#' predicted as the arithmetic mean
#' of performances of all assemblages
#' that share the same assembly motif.\cr
#'
#' Assembly motifs are labelled in the vector \code{assMotif}.
#' Experiments are labelled in the vector \code{xpr}.
#'
#' @usage calibrate_amean_bymot_xpr(fobs, assMotif, xpr)
#'
#' @inheritParams calibrate_byminrss
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the arithmetic mean of performances
#' of all assemblages that share a same assembly motif,
#' experiment by experiment.
#'
#' @details Modelled performances are computed
#' using arithmetic mean (\code{opt.mean = "amean"})
#' of performances of assemblages
#' that share a same assembly motif (\code{opt.model = "bymot"})
#' by including all assemblages that belong to a same assembly motif.
#'
#' @seealso
#' \code{\link{calibrate_amean_bymot_xpr}} arithmetic mean
#'   by assembly motif. \cr
#' \code{\link{calibrate_amean_byelt_xpr}} arithmetic mean
#'   by elements occurring within assembly motif. \cr
#' \code{\link{calibrate_gmean_bymot_xpr}} geometric mean
#'   by assembly motif. \cr
#' \code{\link{calibrate_gmean_byelt_xpr}} geometric mean
#'   by elements occurring within assembly motif. \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

calibrate_amean_bymot_xpr <- function(fobs, assMotif, xpr) {

  fprd <- numeric(length(assMotif))

  setXpr <- unique(names(xpr))
  for (ipr in seq_along(setXpr)) {

    indXpr       <- which(names(xpr) == setXpr[ipr])
    fprd[indXpr] <- calibrate_amean_bymot(fobs[indXpr], assMotif[indXpr])
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
#' @title Modelling of performances
#'  by assembly motif
#'
#' @description Take a vector \code{fobs} of assembly performances
#' and return a vector of performances
#' predicted as the geometric mean
#' of performances of all assemblages
#' that share the same assembly motif. \cr
#'
#' Assembly motifs are labelled in the vector \code{assMotif}.
#'
#' @usage calibrate_gmean_bymot(fobs, assMotif)
#'
#' @inheritParams calibrate_byminrss
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the geometric mean of performances
#' of assemblages that share a same assembly motif.
#'
#' @details Modelled performances are computed
#' using arithmetic mean (\code{opt.mean = "amean"})
#' of performances of assemblages
#' that share a same assembly motif (\code{opt.model = "bymot"})
#' by including all assemblages that belong to a same assembly motif.
#'
#' @seealso
#' \code{\link{calibrate_amean_bymot}} arithmetic mean
#'   by assembly motif. \cr
#' \code{\link{calibrate_amean_byelt}} arithmetic mean
#'   by elements occurring within assembly motif. \cr
#' \code{\link{calibrate_gmean_bymot}} geometric mean
#'   by assembly motif. \cr
#' \code{\link{calibrate_gmean_byelt}} geometric mean
#'   by elements occurring within assembly motif. \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

calibrate_gmean_bymot <- function(fobs, assMotif) {

  fprd <- numeric(length(assMotif))

  setMot <- unique(assMotif)
  for (mot in seq_along(setMot)) {

    indMot         <- which(assMotif == setMot[mot])
    fprd[indMot] <- gmean(fobs[indMot])
  }

  return(fprd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Modelling of performances
#' by assembly motif
#' over several experiments
#'
#' @description Take a vector \code{fobs} of assembly performances
#' over several experiments
#' and return a vector of performances
#' predicted as the geometric mean of performances
#' of all assemblages that share the same assembly motif.\cr
#'
#' Assembly motifs are labelled in the vector \code{assMotif}.
#' Experiments are labelled in the vector \code{xpr}.
#'
#' @usage calibrate_gmean_bymot_xpr(fobs, assMotif, xpr)
#'
#' @inheritParams calibrate_byminrss
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the geometric mean
#' of performances of assemblages
#' that share a same assembly motif,
#' experiment by experiment.
#'
#' @details Modelled performances are computed
#' using geometric mean (\code{opt.mean = "gmean"})
#' of performances of assemblages
#' that share a same assembly motif (\code{opt.model = "bymot"}),
#' experiment by experiment,
#' by including all assemblages that belong to a same assembly motif.
#'
#' @seealso
#' \code{\link{calibrate_amean_bymot_xpr}} arithmetic mean
#'   by assembly motif. \cr
#' \code{\link{calibrate_amean_byelt_xpr}} arithmetic mean
#'   by elements occurring within assembly motif. \cr
#' \code{\link{calibrate_gmean_bymot_xpr}} geometric mean
#'   by assembly motif. \cr
#' \code{\link{calibrate_gmean_byelt_xpr}} geometric mean
#'   by elements occurring within assembly motif. \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

calibrate_gmean_bymot_xpr <- function(fobs, assMotif, xpr) {

  fprd    <- numeric(length(assMotif))

  setXpr <- unique(names(xpr))
  for (ipr in seq_along(setXpr)) {

    indXpr         <- which(names(xpr) == setXpr[ipr])
    fprd[indXpr] <- calibrate_gmean_bymot(fobs[indXpr], assMotif[indXpr])
  }

  return(fprd)
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Arithmetic mean by components occurring within an assembly motif         ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic mean
#'  by components occurring within an assembly motif
#'
#' @description Take a vector \code{fobs}
#' of performances of assemblages
#' that share a same assembly motif,
#' and return a vector of performances
#' predicted as the arithmetic mean
#' of performances of assemblages
#' that contain the same components as the assemblage to predict. \cr
#'
#' @usage amean_byelt(fobs, mOccur)
#'
#' @inheritParams calibrate_byminrss
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the average
#' of mean performances of assemblages that contain the same components
#' as the assemblage to calibrate \code{opt.model = "byelt"} .
#'
#' @details Modelled performances are computed
#' using arithmetic mean (\code{opt.mean = "amean"}) of performances.
#' Assemblages share a same assembly motif.
#' Modelled performances are the average
#' of mean performances of assemblages that contain the same components
#' as the assemblage to calibrate (\code{opt.model = "byelt"}).
#' This procedure corresponds to a linear model with each assembly motif
#' based on the component occurrence in each assemblage.
#'
#' @seealso
#'  \code{\link{amean_byelt}} using arithmetic mean.
#'  \code{\link{gmean_byelt}} using geometric mean.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

amean_byelt <- function(fobs, mOccur) {

  fprd <- numeric(length(fobs))

  setElt <- unique(which((mOccur[ , , drop = FALSE] == 1),
                         arr.ind = TRUE)[ , 2])

  for (elt in seq_along(setElt)) {
    indElt         <- which(mOccur[ , setElt[elt]] == 1)
    fprd[indElt] <- fprd[indElt] + amean(fobs[indElt])
  }

  fprd <- fprd / apply(mOccur, MARGIN = 1, FUN = sum)

  return(fprd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Modelling of performances
#'  by components occurring within an assembly motif
#'
#' @description Take a vector \code{fobs}
#' of performances of assemblages
#' and return a vector of performances
#' predicted as the arithmetic mean
#' of performances of assemblages
#' that share a same assembly motif
#' and that contain the same components as the assemblage to calibrate. \cr
#'
#' The assembly motifs are labelled in the vector \code{assMotif}.
#' The component composition of assemblages is specified
#' in the binary matrix \code{mOccur}:
#' \code{0} if the component does not occur,
#' \code{1} if the component occurs.
#'
#' @usage calibrate_amean_byelt(fobs, assMotif, mOccur)
#'
#' @inheritParams calibrate_byminrss
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the average
#' of mean performances of assemblages
#' that share the same assembly motif
#' and that contain the same components
#' as the assemblage to calibrate \code{opt.model = "byelt"} .
#'
#' @details Modelled performances are computed
#' using arithmetic mean (\code{opt.mean = "amean"}) of performances.
#' Modelled performances are the average
#' of mean performances of assemblages
#' that share a same assembly motif and
#' that contain the same components
#' as the assemblage to calibrate (\code{opt.model = "byelt"}).
#' This procedure corresponds to a linear model with each assembly motif
#' based on the component occurrence in each assemblage.
#' If no assemblage contains component belonging to assemblage to calibrate,
#' performance is the mean performance of all assemblages.
#'
#' @seealso
#' \code{\link{calibrate_amean_bymot}} arithmetic mean
#'   by assembly motif. \cr
#' \code{\link{calibrate_amean_byelt}} arithmetic mean
#'   by elements occurring within assembly motif. \cr
#' \code{\link{calibrate_gmean_bymot}} geometric mean
#'   by assembly motif. \cr
#' \code{\link{calibrate_gmean_byelt}} geometric mean
#'   by elements occurring within assembly motif. \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

calibrate_amean_byelt <- function(fobs, assMotif, mOccur) {

  fprd <- numeric(length(assMotif))

  setMot <- unique(assMotif)
  for (mot in seq_along(setMot)) {

    indMot <- which(assMotif == setMot[mot])

    if (length(indMot) > 1) {
      fprd[indMot] <- amean_byelt(fobs[indMot], mOccur[indMot, ])
    } else {
      fprd[indMot] <- fobs[indMot]
    }
  }

  return(fprd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Modelling of performances
#'  by components occurring within an assembly motif
#'  over several experiments
#'
#' @description Take a vector \code{fobs}
#' of performances of assemblages
#' and return a vector of performances
#' predicted as the arithmetic mean
#' of performances of assemblages
#' that share a same assembly motif
#' and that contain the same components as the assemblage to calibrate. \cr
#'
#' The assembly motifs are labelled in the vector \code{assMotif}.\cr
#'
#' The component composition of assemblages is specified
#' in the binary matrix \code{mOccur}:
#' \code{0} if the component does not occur,
#' \code{1} if the component occurs.
#'
#' @usage calibrate_amean_byelt_xpr(fobs, assMotif, mOccur, xpr)
#'
#' @inheritParams calibrate_byminrss
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the average
#' of mean performances of assemblages
#' that share the same assembly motif
#' and that contain the same components
#' as the assemblage to calibrate \code{opt.model = "byelt"} .
#'
#' @details Modelled performances are computed
#' using arithmetic mean (\code{opt.mean = "amean"}) of performances.
#' Assemblages share a same assembly motif.
#' Modelled performances are the average
#' of mean performances of assemblages that contain the same components
#' as the assemblage to calibrate (\code{opt.model = "byelt"}).
#' This procedure corresponds to a linear model with each assembly motif
#' based on the component occurrence in each assemblage.
#' If no assemblage contains component belonging to assemblage to calibrate,
#' performance is the mean performance of all assemblages.
#'
#' @seealso
#' \code{\link{calibrate_amean_bymot_xpr}} arithmetic mean
#'   by assembly motif. \cr
#' \code{\link{calibrate_amean_byelt_xpr}} arithmetic mean
#'   by elements occurring within assembly motif. \cr
#' \code{\link{calibrate_gmean_bymot_xpr}} geometric mean
#'   by assembly motif. \cr
#' \code{\link{calibrate_gmean_byelt_xpr}} geometric mean
#'   by elements occurring within assembly motif. \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

calibrate_amean_byelt_xpr <- function(fobs, assMotif, mOccur, xpr) {

  fprd <- numeric(length(assMotif))

  setXpr <- unique(names(xpr))
  for (ipr in seq_along(setXpr)) {

    indXpr         <- which(names(xpr) == setXpr[ipr])
    fprd[indXpr] <- calibrate_amean_byelt(fobs[indXpr],
                                          assMotif[indXpr],
                                          mOccur[indXpr, ])
  }

  return(fprd)
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Geometric mean by components occurring within an assembly motif          ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric mean
#'  by components occurring within an assembly motif
#'
#' @description Take a vector \code{fobs}
#' of performances of assemblages
#' that share a same assembly motif,
#' and return a vector of performances
#' predicted as the geometric mean
#' of performances of assemblages
#' that contain the same components as the assemblage to predict. \cr
#'
#' All assemblages share the same assembly motif.
#'
#' @usage gmean_byelt(fobs, mOccur)
#'
#' @inheritParams calibrate_byminrss
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the average
#' of mean performances of assemblages that contain the same components
#' as the assemblage to calibrate (\code{opt.model = "byelt"}) .
#'
#' @details Modelled performances are computed
#' using geometric mean (\code{opt.mean = "gmean"}) of performances.
#' Assemblages share a same assembly motif.
#' Modelled performances are the average
#' of mean performances of assemblages that contain the same components
#' as the assemblage to calibrate (\code{opt.model = "byelt"}).
#' This procedure corresponds to a linear model with each assembly motif
#' based on the component occurrence in each assemblage.
#'
#' @seealso
#'  \code{\link{amean_byelt}} using arithmetic mean.
#'  \code{\link{gmean_byelt}} using geometric mean.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

gmean_byelt <- function(fobs, mOccur) {

  fprd    <- numeric(length(fobs))
  fprd[ ] <- 1

  setElt    <- unique(which(mOccur[ , , drop = FALSE] == 1L,
                            arr.ind = TRUE)[ , 2])
  for (elt in seq_along(setElt)) {
    indElt <- which(mOccur[ , setElt[elt]] == 1)

    if (length(indElt) > 0)
      fprd[indElt] <- fprd[indElt] * gmean(fobs[indElt])
  }

  fprd <- fprd ^ (1 / apply(mOccur, MARGIN = 1, FUN = sum))

  return(fprd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Modelling the performances
#'  by components occurring within an assembly motif
#'
#' @description The numeric vector \code{fobs} get together the properties
#' of assemblages belonging to different assembly motifs. The properties
#' \code{fobs} of assemblages belonging to a given assembly motif and
#' containing a given component are separately averaged. The property of
#'  each assemblage is computed as the average of mean values of assemblages
#'   containing the same components as the considered assemblage. The
#'   motif of each vector component is specified in the vector
#'   \code{assMotif}.
#'    The componental composition of each assemblage is specified in the
#'     binary matrix \code{mOccur}: \code{0} if the component does not occur,
#'      \code{1} if the component occurs.
#'
#' @usage calibrate_gmean_byelt(fobs, assMotif, mOccur)
#'
#' @inheritParams calibrate_byminrss
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the average
#' of mean performances of assemblages
#' that share the same assembly motif
#' and that contain the same components
#' as the assemblage to calibrate \code{opt.model = "byelt"} .
#'
#' @details Modelled performances are computed
#' using geometric mean (\code{opt.mean = "gmean"}) of performances.
#' Modelled performances are the average
#' of mean performances of assemblages
#' that share a same assembly motif and
#' that contain the same components
#' as the assemblage to calibrate (\code{opt.model = "byelt"}).
#' This procedure corresponds to a linear model with each assembly motif
#' based on the component occurrence in each assemblage.
#' If no assemblage contains component belonging to assemblage to calibrate,
#' performance is the mean performance of all assemblages.
#'
#' @seealso
#' \code{\link{calibrate_amean_bymot}} arithmetic mean
#'   by assembly motif. \cr
#' \code{\link{calibrate_amean_byelt}} arithmetic mean
#'   by elements occurring within assembly motif. \cr
#' \code{\link{calibrate_gmean_bymot}} geometric mean
#'   by assembly motif. \cr
#' \code{\link{calibrate_gmean_byelt}} geometric mean
#'   by elements occurring within assembly motif. \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

calibrate_gmean_byelt <- function(fobs, assMotif, mOccur) {

  fprd <- numeric(length(assMotif))

  setMot <- unique(assMotif)
  for (mot in seq_along(setMot)) {

    indMot <- which(assMotif == setMot[mot])
    if (length(indMot) > 1) {
      fprd[indMot] <- gmean_byelt(fobs[indMot], mOccur[indMot, ])
    } else {
      fprd[indMot] <- fobs[indMot]
    }
  }

  return(fprd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Modelling the performances
#'  by components occurring within an assembly motif
#'  over several experiments
#'
#' @description Take a vector \code{fobs}
#' of performances of assemblages
#' and return a vector of performances
#' predicted as the geometric mean
#' of performances of assemblages
#' that share a same assembly motif
#' and that contain the same components as the assemblage to calibrate. \cr
#'
#' The assembly motifs are labelled in the vector \code{assMotif}.\cr
#'
#' The component composition of assemblages is specified
#' in the binary matrix \code{mOccur}:
#' \code{0} if the component does not occur,
#' \code{1} if the component occurs.
#'
#' @usage calibrate_gmean_byelt_xpr(fobs, assMotif, mOccur, xpr)
#'
#' @inheritParams calibrate_byminrss
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed as the average
#' of mean performances of assemblages,
#' experiment by experiment,
#' that share the same assembly motif
#' and that contain the same components
#' as the assemblage to calibrate \code{opt.model = "byelt"} .
#'
#' @details Modelled performances are computed
#' using geometric mean (\code{opt.mean = "gmean"}) of performances.
#' Modelled performances are the average,
#' experiment by experiment,
#' of mean performances of assemblages
#' that share a same assembly motif and
#' that contain the same components
#' as the assemblage to calibrate (\code{opt.model = "byelt"}).
#' This procedure corresponds to a linear model with each assembly motif
#' based on the component occurrence in each assemblage.
#' If no assemblage contains component belonging to assemblage to calibrate,
#' performance is the mean performance of all assemblages.
#'
#' @seealso
#' \code{\link{calibrate_amean_bymot_xpr}} arithmetic mean
#'   by assembly motif. \cr
#' \code{\link{calibrate_amean_byelt_xpr}} arithmetic mean
#'   by elements occurring within assembly motif. \cr
#' \code{\link{calibrate_gmean_bymot_xpr}} geometric mean
#'   by assembly motif. \cr
#' \code{\link{calibrate_gmean_byelt_xpr}} geometric mean
#'   by elements occurring within assembly motif. \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

calibrate_gmean_byelt_xpr <- function(fobs, assMotif, mOccur, xpr) {

  fprd <- numeric(length(assMotif))

  setXpr <- unique(names(xpr))
  for (ipr in seq_along(setXpr)) {

    index       <- which(names(xpr) == setXpr[ipr])
    fprd[index] <- calibrate_gmean_byelt( fobs[index],
                                          assMotif[index],
                                          mOccur[index, ] )
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
#' @title Modelling of performances of assemblages
#'
#' @description Take a vector \code{fobs} of assembly performances
#' over several experiments
#' and return a vector of performances
#' predicted as the mean performances of assemblages
#' that share the same assembly motif.\cr
#'
#' Assembly motifs are labelled in the vector \code{assMotif}.
#' Experiments are labelled in the vector \code{xpr}.
#' Modelling options are indicated in \code{opt.mean} and \code{opt.model}.
#' Occurrence matrix \code{mOccur} is used if \code{opt.model =  "byelt"}.
#'
#' @usage
#' calibrate_byminrss(fobs, assMotif, mOccur, xpr,
#' opt.mean = "amean", opt.model = "bymot"  )
#'
#' @param fobs a numeric vector. The vector \code{fobs} contains the
#' quantitative performances of assemblages.
#'
#' @param assMotif a vector of labels of \code{length(fobs)}.
#' The vector \code{assMotif} contains the assembly motifs of assemblages.
#'
#' @param mOccur a matrix of occurrence (occurrence of components).
#' Its first dimension equals to \code{length(fobs)}. Its second dimension
#'  equals to the number of components.
#'
#' @param xpr a vector of numerics of \code{length(fobs)}.
#' The vector \code{xpr} contains the weight of each experiment,
#' and the labels (in \code{names(xpr)}) of different experiments.
#' The weigth of each experiment is used
#' in the computation of the Residual Sum of Squares
#' in the function \code{rss_clustering}.
#' All assemblages that belong to a given experiment
#' should then have a same weigth.
#' Each experiment is identified by its names (\code{names(xpr)})
#' and the RSS of each experiment is weighted by values of \code{xpr}.
#' The vector \code{xpr} is generated
#' by the function \code{stats::setNames}.
#'
#' @param opt.mean a character equals to \code{"amean"} or \code{"gmean"}.
#' Switchs to arithmetic formula if \code{opt.mean = "amean"}.
#' Switchs to geometric formula if \code{opt.mean = "gmean"}.
#'
#' @param opt.model a character equals to \code{"bymot"} or \code{"byelt"}.
#' Switchs to simple mean by assembly motif if \code{opt.model = "bymot"}.
#' Switchs to linear model with assembly motif if \code{opt.model = "byelt"}.
#'
#' @return Return a vector of \code{length(fobs)}.
#' Its values are computed according to \code{opt.mean} and \code{opt.model}.
#'
#' @details Modelled performances are computed
#' using arithmetic mean (\code{opt.mean = "amean"})
#' or geometric mean (\code{opt.mean = "gmean"}).\cr
#'
#' If \code{opt.model = "bymot"},
#' modelled performances are means
#' of performances of assemblages
#' that share a same assembly motif
#' by including all assemblages that belong to a same assembly motif. \cr
#'
#' If \code{opt.model = "byelt"},
#' modelled performances are the average
#' of mean performances of assemblages
#' that share a same assembly motif
#' and that contain the same components
#' as the assemblage to calibrate.
#' This procedure corresponds to a linear model within each assembly motif
#' based on the component occurrence in each assemblage.
#' If no assemblage contains component belonging to assemblage to calibrate,
#' performance is the mean performance of all assemblages
#' as in \code{opt.model = "bymot"}.
#'
#' @seealso
#' \code{\link{validate_using_cross_validation}}
#'    predicts performances of assemblages.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

calibrate_byminrss <- function(fobs, assMotif, mOccur, xpr,
                               opt.mean  = "amean",
                               opt.model = "bymot") {

  optmean <- "amean"
  if (opt.mean == "gmean") optmean <- "gmean"

  optmod <- "bymot"
  if (opt.model == "byelt") optmod <- "byelt"

  optxpr <- ""
  if (length(unique(names(xpr))) != 1) optxpr <- "xpr"

  option <- paste(optmean, optmod, optxpr, sep = ".")

  return(
    switch(option,
           amean.bymot. =
             calibrate_amean_bymot(fobs, assMotif),
           amean.bymot.xpr =
             calibrate_amean_bymot_xpr(fobs, assMotif, xpr),
           gmean.bymot. =
             calibrate_gmean_bymot(fobs, assMotif),
           gmean.bymot.xpr =
             calibrate_gmean_bymot_xpr(fobs, assMotif, xpr),

           amean.byelt. =
             calibrate_amean_byelt(fobs, assMotif, mOccur),
           amean.byelt.xpr =
             calibrate_amean_byelt_xpr(fobs, assMotif, mOccur, xpr),
           gmean.byelt. =
             calibrate_gmean_byelt(fobs, assMotif, mOccur),
           gmean.byelt.xpr =
             calibrate_gmean_byelt_xpr(fobs, assMotif, mOccur, xpr)  )
  )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# End of file                                                              ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
