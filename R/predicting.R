#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                               PREDICTING.R
#
# Predicting performances of new assemblages                               ####
#      by knowing their elemental composition
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



#' @include
#'  stats.R
#'  validating.R
#'
NULL


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Arithmetic mean by Motif                                                 ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Supplementary assemblages to predict                                     ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Prediction of supplementary assemblages by motif
#'
#' @description Takes a numeric f.
#'
#' @usage predict_amean_bymot(appFct, appMotifs, supMotifs)
#'
#' @inheritParams predict_performance
#'
#' @return Return a vector of \code{length(supMotifs)}. The values are
#'  computed using arithmetic mean
#'  of components belonging to \code{appMotifs}
#'  and sharing a same motif.
#'
#' @details Prediction ...
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict_amean_bymot <- function(appFct, appMotifs, supMotifs) {

  supFct   <- numeric(length(supMotifs))
  supFct[] <- NA

  setMot <- unique(supMotifs)
  for (mot in seq_along(setMot)) {

    indSup <- which(supMotifs == setMot[mot])
    indApp <- which(appMotifs == setMot[mot])
    if ( (length(indSup) > 0) & (length(indApp) > 0) )
      supFct[indSup] <- amean(appFct[indApp])
  }

  return(supFct)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Prediction of supplementary assemblages
#'
#' @description Takes a numeric f.
#'
#' @usage predict_gmean_bymot(appFct, appMotifs, supMotifs)
#'
#' @inheritParams predict_performance
#'
#' @return Return a vector of \code{length(supMotifs)}. The values are
#'  computed using arithmetic mean of components belonging to \code{appMotifs}
#'  and sharing a same motif.
#'
#' @details Prediction ...
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict_gmean_bymot <- function(appFct, appMotifs, supMotifs) {

  supFct   <- numeric(length(supMotifs))
  supFct[] <- NA

  setMot <- unique(supMotifs)
  for (mot in seq_along(setMot)) {

    indSup <- which(supMotifs == setMot[mot])
    indApp <- which(appMotifs == setMot[mot])
    if ( (length(indSup) > 0) & (length(indApp) > 0) )
      supFct[indSup] <- gmean(appFct[indApp])
  }

  return(supFct)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Prediction of supplementary assemblages computed
#'
#' @description Takes a numeric f.
#'
#' @usage
#' predict_amean_byelt(appFct, appMotifs, appOccur,
#'                         supMotifs, supOccur )
#'
#' @inheritParams predict_performance
#'
#' @details dd
#'
#' @return  cccc
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict_amean_byelt <- function(appFct, appMotifs, appOccur,
                                supMotifs, supOccur ) {

  setAppMot <- unique(appMotifs)
  setSupMot <- unique(supMotifs)
  setMot    <- sort(union(setAppMot, setSupMot))
  nbMot     <- length(setMot)

  mfct      <- matrix(NA, nrow = nbMot, ncol = dim(appOccur)[2],
                      dimnames = list(setMot, colnames(appOccur)))

  for (mot in seq_along(setAppMot)) {

    motif  <- setAppMot[mot]
    indApp <- which(appMotifs == motif)
    setElt <- unique(which((appOccur[indApp, , drop = FALSE] == 1),
                           arr.ind = TRUE)[ , 2])

    for (elt in seq_along(setElt)) {
      element <- setElt[elt]
      indElt  <- which(appOccur[indApp, element] == 1)
      if (length(indElt) > 0)
        mfct[motif, element] <- amean(appFct[indApp[indElt]])
    }
  }


  supFct  <- numeric(length(supMotifs))
  sizeSup <- apply(supOccur, MARGIN = 1, FUN = sum)

  for (mot in seq_along(setSupMot)) {
    motif     <- setSupMot[mot]
    indSupMot <- which(supMotifs == motif)

    if (length(indSupMot) > 0) {
      setSupElt <- unique(which((supOccur[indSupMot, , drop = FALSE] == 1),
                                arr.ind = TRUE)[ , 2])

      for (elt in seq_along(setSupElt)) {
        element   <- setSupElt[elt]
        indSupElt <- which(supOccur[indSupMot, element] == 1)

        if (length(indSupElt) > 0) {
          index         <- indSupMot[indSupElt]
          supFct[index] <- supFct[index] + mfct[motif, element]
        }
      }

      supFct[indSupMot] <- supFct[indSupMot] / sizeSup[indSupMot]
    }
  }

  return(supFct)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title  Prediction of supplementary assemblages computed
#'      gmean = by using geometric mean
#'      byelt = by motif WITH taking into account species contribution
#'      by including all the assemblages, even the one to predict
#'      for any Function (for instance Fobs)
#'
#' @description Takes a numeric f.
#'
#' @usage
#' predict_gmean_byelt(appFct, appMotifs, appOccur,
#'                         supMotifs, supOccur )
#'
#' @inheritParams predict_performance
#'
#' @details dd
#'
#' @return  cccc
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict_gmean_byelt <- function(appFct, appMotifs, appOccur,
                                supMotifs, supOccur  ) {

  setAppMot <- unique(appMotifs)
  setSupMot <- unique(supMotifs)
  setMot    <- sort(union(setAppMot, setSupMot))
  nbMot     <- length(setMot)

  mfct      <- matrix(NA, nrow = nbMot, ncol = dim(appOccur)[2],
                      dimnames = list(setMot, colnames(appOccur)))

  for (mot in seq_along(setAppMot)) {

    motif  <- setAppMot[mot]
    indApp <- which(appMotifs == motif)
    setElt <- unique(which((appOccur[indApp, , drop = FALSE] == 1),
                           arr.ind = TRUE)[ , 2])

    for (elt in seq_along(setElt)) {
      element <- setElt[elt]
      indElt  <- which(appOccur[indApp, element] == 1)
      if (length(indElt) > 0)
        mfct[motif, element] <- gmean(appFct[indApp[indElt]])
    }
  }


  supFct   <- numeric(length(supMotifs))
  supFct[] <- 1
  sizeSup  <- apply(supOccur, MARGIN = 1, FUN = sum)

  for (mot in seq_along(setSupMot)) {

    motif     <- setSupMot[mot]
    indSupMot <- which(supMotifs == motif)

    if (length(indSupMot) > 0) {
      setSupElt <- unique(which((supOccur[indSupMot, , drop = FALSE] == 1),
                                arr.ind = TRUE)[ , 2])

      for (elt in seq_along(setSupElt)) {
        element   <- setSupElt[elt]
        indSupElt <- which(supOccur[indSupMot, element] == 1)

        if (length(indSupElt) > 0) {
          index         <- indSupMot[indSupElt]
          supFct[index] <- supFct[index] * mfct[motif, element]
        }
      }

      supFct[indSupMot] <- supFct[indSupMot] ^ (1/sizeSup[indSupMot])
    }
  }

  return(supFct)
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Functions for switch on different options                                ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Predicting performances of assemblages
#' by only knowing their elemental composition
#'
#' @description Takes a vector \code{fct} of assembly performances
#' over several experiments
#' and returns a vector of performances
#' predicted as the mean performances of assemblages
#' that share the same assembly motif. \cr
#'
#' Assembly motifs are labelled in the vector \code{assMotif}.
#' Experiments are labelled in the vector \code{xpr}.
#' Modelling options are indicated in \code{opt.mean} and \code{opt.model}.
#' Occurrence matrix \code{mOccur} is used if \code{opt.model =  "byelt"}.
#' Cross-validation is leave-one-out or jackknifesi
#'
#' @usage
#' predict_performance(appFct, appMotifs, appOccur,
#'             supMotifs, supOccur,
#'             opt.mean = "amean",
#'             opt.model  = "bymot"  )
#'
#' @param appFct a vector of numeric values (assembly properties).
#'
#' @param appMotifs a vector of labels of \code{length(fct)} (assembly motifs).
#'
#' @param appOccur a matrix of occurrence (occurrence of components).
#' Its first dimension equals to \code{length(fct)}. Its second dimension
#'  equals to the number of components.
#'
#' @param supMotifs a vector of labels of \code{length(fct)} (assembly motifs).
#'
#' @param supOccur a matrix of occurrence (occurrence of components).
#' Its first dimension equals to \code{length(fct)}. Its second dimension
#'  equals to the number of components.
#'
#' @param opt.mean equal to \code{"amean"} (by default) or \code{"gmean"}.
#'
#' @param opt.model  equal to \code{"bymot"} (by default) or \code{"byelt"}.
#'
#' @return Return the arithmetic mean of a vector, as standard \code{mean}
#'  function.
#'
#' @details Prediction is computed using arithmetic mean \code{amean} by motif
#' \code{bymot} in a whole (WITHOUT taking into account species contribution).
#' The components belonging to a same motif are divided into \code{jack[2]}
#' subsets of \code{jack[1]} components. Prediction is computed by excluding
#' \code{jack[1]} components, of which the component to predict. If the total
#' number of components belonging to the motif is lower than
#' \code{jack[1]*jack[2]}, prediction is computed by Leave-One-Out (LOO).
#'
#' @seealso
#' \code{\link{calibrate_byminrss}} \cr
#' \code{\link{validate_using_cross_validation}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict_performance <- function(appFct, appMotifs, appOccur,
                                supMotifs, supOccur,

                                opt.mean = "amean",
                                opt.model  = "bymot"  ) {

  option <- paste(opt.mean, opt.model, sep = ".")

  return(
    switch(option,
           amean.bymot =
             predict_amean_bymot(appFct, appMotifs, supMotifs) ,
           gmean.bymot =
             predict_gmean_bymot(appFct, appMotifs, supMotifs) ,

           amean.byelt =
             predict_amean_byelt(appFct, appMotifs, appOccur,
                                 supMotifs, supOccur) ,
           gmean.byelt =
             predict_gmean_byelt(appFct, appMotifs, appOccur,
                                 supMotifs, supOccur)
    )
  )
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Main function for predicting performances of new assemblages             ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

###############################################################################
#'
#' @title nnnn
#' @description vvv
#'
#' @param tree vv
#' @param res.prd ff
#' @param supOccur ff
#' @param appOccur ff
#' @param appFct ff
#' @param opt.mean ff
#' @param opt.model ff
#'
#' @details ccc
#' @return vv
#' @keywords internal
#'
###############################################################################

fclust_predict <- function(tree, res.prd,
                               supOccur, appOccur, appFct,
                               opt.mean   = "amean",
                               opt.model  = "byelt"    ) {

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  #  Main calculations of Cal, Prd, R2cal and R2prd
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #  Cas 1 où on suppose le nombre d'éléments identiques
  #      dans apprentissage et supplémentaires
  #  Prévoir le cas 2 où le nombre d'éléments est différent
  #      dans apprentissage et supplémentaire

  nbElt      <- dim(appOccur)[2]    # ????

  nbAppAss   <- dim(appOccur)[1]
  nbSupAss   <- dim(supOccur)[1]

  mAssMotifs <- maffect_motifs(tree, rbind(appOccur, supOccur))
  mAppMotifs <- mAssMotifs[ , 1:nbAppAss]
  mSupMotifs <- mAssMotifs[ , (nbAppAss + 1):(nbAppAss + nbSupAss)]

  # Compute the raw predictions
  mSup <- mSd <- tNbcl <-
    matrix(NA, nrow = nbElt, ncol = nbSupAss,
           dimnames = list(seq_len(nbElt), rownames(supOccur)))
  storage.mode(tNbcl) <- "integer"

  for (nbcl in seq_len(nbElt)) {

    mSup[nbcl, ] <- predict_performance(appFct, mAppMotifs[nbcl, ], appOccur,
                                        mSupMotifs[nbcl, ], supOccur,
                                        opt.mean, opt.model)

    mSd[nbcl, ] <- res.prd$uRmse[[nbcl]][mSupMotifs[nbcl, ]]

    index <- which(mSup[nbcl, ] > appFct)
    mSd[nbcl, index] <- -mSd[nbcl, index]
  }


  # Compute the associated statistiques
  tNbcl[1, ] <- 1L
  for (nbcl in 2:nbElt) {
    tNbcl[nbcl, ] <- nbcl
    tNbcl[nbcl, is.na(mSup[nbcl, ])] <- tNbcl[(nbcl - 1L), is.na(mSup[nbcl, ])]
  }

  tSup <- mSup
  for (nbcl in 2:nbElt)
    tSup[nbcl, is.na(mSup[nbcl, ])] <- tSup[(nbcl - 1L), is.na(mSup[nbcl, ])]


  # Compute the associated statistiques
  nbK <- integer(nbElt)
  for (nbcl in seq_len(nbElt))
    nbK[nbcl] <- length(unique(mSupMotifs[nbcl, ]))

  mStats <- compute_fit_stats(mSup, mSup, appFct, nbK)


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  #  Outputs
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  #                  "tError" A vérifier

  res        <- list(rownames(supOccur), appFct, opt.mean, opt.model,
                     mSupMotifs, tSup, mSd, mStats, tNbcl)
  names(res) <- c("names", "fct", "opt.mean", "opt.model",
                  "mMotifs", "tSup", "tSd", "tStats", "tNbcl")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# End of file                                                              ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
