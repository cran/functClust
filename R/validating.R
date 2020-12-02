#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                             PREDICTING_COMPUT.R
#
# Set of functions for computing
#         Predictions and associated Statistics                            ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#' @include
#'   stats.R
#'   labelling.R
#'   calibrating.R
#'   validating_loo.R
#'   validating_jack.R
#'
NULL


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Functions for computing statistics                                       ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Statistics of model goodness-of-fit
#'
#' @description Taks a matrix of calibrations, a matrix of predictions,
#'   the vector of observed performances,
#'   the number of observed assembly motifs,
#'   and return a matrix of statistics for model goodness-of-fit.
#'
#' @usage compute_fit_stats(mCal, mPrd, fobs, xpr, nbK)
#'
#' @inheritParams validate_ftree
#'
#' @param mCal a numeric matrix.
#' This matrix is the matrix of performances predicted by the tree model.
#'
#' @param mPrd a numeric matrix.
#' This matrix is the matrix of performances predicted by cross-validation.
#'
#' @param nbK an integer.
#' This integer corresponds to the number of observed assembly motifs.
#'
#' @details Be careful, the matrix order is not ramdon.
#' The first argument \code{mCal} is matrix of modelled values.
#' The second argument \code{mPrd} is matrix of values
#' predicted by cross-validation.
#' The third argument \code{fobs} is the vector of observed values.
#'
#' @return Return statistics of model goodness-of-fit.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

compute_fit_stats <- function(mCal, mPrd, fobs, xpr, nbK) {

  nbClu <- dim(mCal)[1]

  ntmp   <- c("missing", "R2cal", "R2prd", "AIC", "AICc")
  mStats <- matrix(0, nrow = nbClu, ncol = length(ntmp),
                   dimnames = list(seq_len(nbClu), ntmp))

  setXpr <- unique(names(xpr))
  wmp    <- numeric(length(setXpr))
  tmp    <- matrix(0, nrow = length(setXpr), ncol = length(ntmp),
                   dimnames = list(setXpr, ntmp))

  for (nbcl in seq_len(nbClu)) {

    mStats[nbcl, "missing"] <- sum(is.na(mPrd[nbcl, ])) / length(mCal[nbcl, ])
    mStats[nbcl, "R2cal"]   <- R2mse(mCal[nbcl, ], fobs)
    mStats[nbcl, "R2prd"]   <- R2mse(mPrd[nbcl, ], fobs)

    # Calculation of AIC as the median of AIC of each performance.
    # Because each performance is a pseudoreplication

    for (ipr in seq_along(setXpr)) {
      index            <- which(names(xpr) == setXpr[ipr])
      wmp[ipr]         <- unique(xpr[index])
      tmp[ipr, "AIC"]  <- AIC_( mCal[nbcl, index], fobs[index], nbK[nbcl])
      tmp[ipr, "AICc"] <- AICc( mCal[nbcl, index], fobs[index], nbK[nbcl])
    }
    mStats[nbcl, "AIC"]  <- median(tmp[ , "AIC"],  each = wmp)
    mStats[nbcl, "AICc"] <- median(tmp[ , "AICc"], each = wmp)
  }

  return(mStats)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#' @title Valid hierarchical tree and Statistics of model goodness-of-fit
#'
#' @description Take a matrix of calibrations, a matrix of predictions,
#'   the vector of observed performances,
#'   the number of observed assembly motifs,
#'   and return a matrix of statistics for model goodness-of-fit.
#'
#' @usage compute_ftree_stats(mCal, mPrd, mStats, fobs, xpr, nbK)
#'
#' @inheritParams validate_ftree
#'
#' @param mCal a numeric matrix.
#' This matrix is the matrix of performances predicted by the model.
#'
#' @param mPrd a numeric matrix.
#' This matrix is the matrix of performances predicted by cross-validation.
#'
#' @param mStats a numeric matrix.
#' This matrix is the matrix of statistics of model goodness-of-fit.
#'
#' @param fobs a numeric vector.
#' This vector is the vector of observed performances.
#'
#' @param nbK an integer.
#' This integer corresponds to the number of observed assembly motifs.
#'
#' @details Be careful, the matrix order is not ramdon.
#' The first argument \code{mCal} is matrix of modelled values.
#' The second argument \code{mPrd} is matrix of values
#' predicted by cross-validation.
#' The fourth argument \code{fobs} is the vector of observed values.
#'
#' @return
#'
#' \itemize{
#'  \item \code{tCal}: a matrix of the valid part of hierarchical tree,
#'     that is the part of tree that increases predictive ability of model,\cr
#'  \item \code{tCal} and \code{tPrd}: the valid part of hierarchical tree,
#'     that is the part of tree that increases predictive ability of model,\cr
#'  \item \code{tStats}: statistics of tree model goodness-of-fit,\cr
#'  \item \code{tNbcl}: the number of clusters used or
#'  computing each performance.
#'}
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

compute_ftree_stats <- function(mCal, mPrd, mStats, fobs, xpr, nbK) {

  nbClu <- dim(mCal)[1]
  nbAss <- dim(mCal)[2]

  tCal <- tPrd <- tNbcl <-
    matrix(NA, nrow = nbClu, ncol = nbAss,
           dimnames = list(seq_len(nbClu), names(fobs)))

  for (nbcl in seq_len(nbClu)) {
    for (ass in seq_len(nbAss)) {

      last  <- min(sum(!is.na(mPrd[1:nbcl, ass])),
                   sum(!is.na(mStats[1:nbcl, "R2prd"])),
                   nbcl)
      index <- first_argmax(mStats[1:last, "R2prd"])

      tCal[nbcl, ass]  <- mCal[index, ass]
      tPrd[nbcl, ass]  <- mPrd[index, ass]
      tNbcl[nbcl, ass] <- index
    }
  }
  tStats <- compute_fit_stats(tCal, tPrd, fobs, xpr, nbK)


  if (dim(tStats)[1] > 1)
    for (nbcl in 2:nbClu)
      if (tStats[nbcl, "R2prd"] <= tStats[nbcl - 1, "R2prd"]) {

        tCal[nbcl, ]  <- tCal[nbcl - 1, ]
        tPrd[nbcl, ]  <- tPrd[nbcl - 1, ]
        tNbcl[nbcl, ] <- tNbcl[nbcl - 1, ]
      }
  tStats <- compute_fit_stats(tCal, tPrd, fobs, xpr, nbK)

  res        <- list(tCal, tPrd, tStats, tNbcl)
  names(res) <- c("tCal", "tPrd", "tStats", "tNbcl")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#' @title Statistics of assembly motifs
#'
#' @description Take a matrix of tree calibrations,
#'  a matrix of tree predictions,
#'  a matrix of assembly motifs,
#'  and return statistics of each observed assembly motif.
#'
#' @usage compute_motif_stats(tCal, tPrd, mMotifs)
#'
#' @param tCal a numeric matrix.
#' This matrix is the matrix of performances predicted by the valid tree model.
#'
#' @param tPrd a numeric matrix.
#' This matrix is the matrix of performances predicted by cross-validation.
#'
#' @param mMotifs a numeric matrix.
#' This matrix is the matrix of assembly motifs.
#'
#' @details The different assembly motifs have different length:
#' the motif set can be treated as a list.
#' Each assembly motif is separately analysed.
#'
#' \itemize{
#'  \item \code{uTab}: the assemblages that belong to the assembly motif,
#'  \item \code{uMean}: the arithmetic mean of motif performances,
#'  \item \code{uSd}: the standard deviation of motif performances,
#'  \item \code{uRmse}: the Root Mean Square Error of motif performances,
#'  \item \code{uR2}: the Coefficient of determination of motif performances,
#'  \item \code{uSlope}:
#'  the slope of linear regression with motif performances.
#'}
#'
#' @return Returns the statistics for each assembly motif.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

compute_motif_stats <- function(tCal, tPrd, mMotifs) {

  nbClu <- dim(tCal)[1]

  # build the general structure of motif_stats

  uTab <- list()
  for (nbcl in 1:nbClu) uTab[[nbcl]] <- table(mMotifs[nbcl, ])

  # compute the values of motif_Stats

  uMean <- uSd <- uRmse <- uR2 <- uSlope <- uTab

  for (nbcl in 1:nbClu) {

    for (motif in names(uTab[[nbcl]])) {

      index <- which(mMotifs[nbcl, ] == motif, arr.ind = TRUE)
      uMean[[nbcl]][motif] <- amean(tPrd[nbcl, index])

      if (length(index) > 1) {

        uSd[[nbcl]][motif]    <- asd(tPrd[nbcl, index])
        uRmse[[nbcl]][motif]  <- rmse(tPrd[nbcl, index], tCal[nbcl, index])
        uR2[[nbcl]][motif]    <- R2mse(tPrd[nbcl, index], tCal[nbcl, index])

#        if (length(!is.na(tPrd[nbcl, index])) > 1) {
#          uSlope[[nbcl]][motif] <-
#            stats::lm(tPrd[nbcl, index] ~ tCal[nbcl, index])$coef[2]
##        } else {
#"          uSlope[[nbcl]][motif] <- uSlope[[nbcl - 1]][oldMotif]
#        }

      } else {

        oldMotif <- mMotifs[nbcl - 1, index]

#        uSd[[nbcl]][motif]    <- uSd[[nbcl - 1]][oldMotif]
#        uRmse[[nbcl]][motif]  <- uRmse[[nbcl - 1]][oldMotif]
#        uR2[[nbcl]][motif]    <- uR2[[nbcl - 1]][oldMotif]
#        uSlope[[nbcl]][motif] <- uSlope[[nbcl - 1]][oldMotif]
      }
    }
  }

  res        <- list(uTab, uMean, uSd, uRmse, uR2, uSlope)
  names(res) <- c("uTab", "uMean", "uSd", "uRmse", "uR2", "uSlope")

  return(res)
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Main functions                                                           ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Predictions of assembly performances
#'        using a species clustering tree
#'
#' @description Take a hierarchical tree of species clustering,
#'   a matrix of occurrency and the corresponding vector of performances,
#'   and return the predictions, statistics and other informations.
#'
#' @usage
#' validate_ftree(tree.I, fobs, mOccur,
#'               xpr = stats::setNames(rep(1, length(fobs)),
#'                                     rep("a", length(fobs))),
#'               opt.method = "divisive", opt.mean = "amean",
#'               opt.model = "byelt",
#'               opt.jack   = FALSE,
#'               jack       = as.integer(c(3, 4)),
#'               opt.nbMax  = dim(mOccur)[2])
#'
#' @param tree.I an integer square-matrix.
#'  The matrix represents a hierarchical tree of species clustering.
#'
#' @param fobs a numeric vector. The vector \code{fobs} contains the
#' quantitative performances of assemblages.
#'
#' @param mOccur a matrix of occurrence (occurrence of elements).
#'  Its first dimension equals to \code{length(fobs)}. Its second dimension
#'   equals to the number of elements.
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
#' The vector \code{xpr} is generated
#' by the function \code{stats::setNames}.
#'
#' @param opt.method a string that specifies the method to use.
#' \code{opt.method = c("sort", "divisive", "agglomerative", "cluster")}.
#' The three first methods generate hierarchical trees.
#' Each tree is complete, running from a unique trunk
#' to as many leaves as components.
#' The last method generates, at each level of the tree,
#' a clustering of components into a given, predefined number of clusters.
#' Because it is repeated from the trunk until to leaves,
#' by increasing the number of clusters,
#' the method generates a non-hierarchical tree. \cr
#'
#' If \code{opt.method = "sort"}, the components are sorted
#' by their effect of assemblage performances. \cr
#'
#' If \code{opt.method = "divisive"}, the components are clustered
#' according to a hierarchical process
#' by using a divisive method,
#' from the trivial cluster where all components are together,
#' towards the clustering where each component is a cluster. \cr
#'
#' If \code{opt.method = "agglomerative"}, the components are clustered
#' according to a hierarchical process
#' by using an agglomerative method,
#' from the trivial clustering where each component is a clsuter,
#' towards the cluster where all components are together.
#' The method that gives the best result is \code{opt.method = "divisive"}. \cr
#'
#' If \code{opt.method = "cluster"}, the components are clustered
#' according to a non-hierarchical process
#' by using the method of McNaughton-Smith et al., 1964.
#' In this case, one must specify the number of wished clusters.
#'
#' Recall that, if \code{affectElt} is specified,
#' the option \code{opt.method} does not need to be filled out.
#' \code{affectElt} determines a level of component clustering,
#' and a tree is built:
#' \emph{(i)} by using \code{opt.method =  "divisive"}
#' from the defined level in tree towards as many leaves as components;
#' \emph{(ii)} by using \code{opt.method =  "agglomerative"}
#' from the defined level in tree towards the trunk of tree.
#'
#' @param opt.mean a character equals to \code{"amean"} or \code{"gmean"}.
#' Switch to arithmetic formula if \code{opt.mean = "amean"}.
#' Switch to geometric formula if \code{opt.mean = "gmean"}.
#'
#' @param opt.model a character equals to \code{"bymot"} or \code{"byelt"}.
#' Switch to simple mean by assembly motif if \code{opt.model = "bymot"}.
#' Switch to linear model with assembly motif if \code{opt.model = "byelt"}.
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
#' @param opt.nbMax an integer, that indicates the maximum number
#' of tree levels to cluster.
#' By default, \code{opt.nbMax = nbElt} for clustering components
#' all along the tree, from the trunc to the leaves, to be able to determine
#' the optimum number of component functional groups.
#' However, in \code{ftest_*} and \code{fboot_*} functions,
#' there is no point in cluster the components
#' beyond the optimum number of functional groups. In these functions,
#' \code{opt.nbMax = } optimum number of functional groups, by default.
#'
#' @details None.
#'
#' @return Return a list containing predictions of assembly performances
#'       and statistics computed by using a species clustering tree.
#'
#'  Recall of inputs:
#'  \itemize{
#'  \item \code{nbElt, nbAss}{: the numbers of components, of assemblages}
#'  \item \code{opt.method}{: the method used to cluster components,}
#'  \item \code{opt.mean}{: the option for mean values computing,}
#'  \item \code{opt.model}{: the option for prediction modelling,}
#'  \item \code{opt.jack}{: the option for method of cross-validation,}
#'  \item \code{jack}{: the parameters for jackknife,}
#'  \item \code{fobs}{: the vector of observed performances of assemblages,}
#'  \item \code{mOccur}{: the matrix of component occurrence,}
#'  \item \code{xpr}{: the vector of labels of different experiments.}
#'  }
#'
#'  Primary and secondary trees of element clustering:
#'  \itemize{
#'  \item \code{tree.I}{: the primary tree of component clustering,}
#'  \item \code{tree.II}{: the validated secondary
#'   tree of component clustering,}
#'  \item \code{nbOpt}{: the optimum number of clusters,}
#'  }
#'
#'  Matrices of calibration and prediction using tree.I
#'    and associated statistics:
#'  \itemize{
#'  \item \code{mCal}{: the matrix of modelled values,}
#'  \item \code{mPrd}{: the matrix of values predicted by cross-validation,}
#'  \item \code{mMotifs}{: the matrix of labels of assembly motifs,}
#'  \item \code{mStats}{: the matrix of associated statistics.}
#'  }
#'
#'  Matrices of calibaration and prediction using tree.II
#'    and associated statistics:
#'  \itemize{
#'  \item \code{tCal}{: the matrix of values modelled
#'    using the valid part of tree,}
#'  \item \code{tPrd}{: the matrix of values predicted
#'    using the valid part of tree,}
#'  \item \code{tStats}{: statistics of valid tree model goodness-of-fit,}
#'  \item \code{tNbcl}{: the number of clusters used
#'    or computing each performance.}
#'  }
#'
#' @importFrom stats setNames
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

validate_ftree <- function(tree.I, fobs, mOccur,
                           xpr = stats::setNames(rep(1, length(fobs)),
                                                 rep("a", length(fobs))),
                           opt.method = "divisive",
                           opt.mean   = "amean",
                           opt.model  = "byelt",
                           opt.jack   = FALSE,
                           jack       = as.integer(c(3, 4)),
                           opt.nbMax  = dim(mOccur)[2]) {

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  #  Main computations of Cal, Prd, R2cal and R2prd
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  if (getOption("verbose") == TRUE)
    cat("\rPruning of primary tree => validated secondary tree\n")

  nbElt <- nbClu <- as.integer(dim(mOccur)[2])
  nbAss <- as.integer(length(unique(rownames(mOccur))))
  nbXpr <- as.integer(length(unique(names(xpr))))

  # Compute the predictions by cross-validation
  nbK  <- integer(nbClu)
  mCal <- mPrd <- mMot <-
    matrix(NA, nrow = nbClu, ncol = nbXpr * nbAss,
           dimnames = list(seq_len(nbClu), rownames(mOccur)))
  storage.mode(mMot) <- "integer"

  nbcl <- 1
  for (nbcl in seq_len(opt.nbMax)) {

    mMot[nbcl, ] <- affect_motifs(tree.I$aff[nbcl, ], mOccur)
    nbK[nbcl]    <- length(unique(mMot[nbcl, ]))

    mCal[nbcl, ] <- calibrate_byminrss(fobs, mMot[nbcl, ], mOccur, xpr,
                                       opt.mean, opt.model)
    mPrd[nbcl, ] <- validate_using_cross_validation(fobs, mMot[nbcl, ],
                                                    mOccur, xpr,
                                                    opt.mean, opt.model,
                                                    opt.jack, jack)
  }

  if (opt.nbMax < nbClu)
    for (nbcl in (opt.nbMax + 1):nbClu) {
      mMot[nbcl, ] <- mMot[opt.nbMax, ]
      nbK[nbcl]    <- nbK[opt.nbMax]
      mCal[nbcl, ] <- mCal[opt.nbMax, ]
      mPrd[nbcl, ] <- mPrd[opt.nbMax, ]
    }


  # Compute the associated statistiques: 1. of primary, fitted tree
  #                                      2. of secondary, validated tree
  #                                      3. for each assembly motif
  mStats  <- compute_fit_stats(mCal, mPrd, fobs, xpr, nbK)
  res.prd <- compute_ftree_stats(mCal, mPrd, mStats, fobs, xpr, nbK)
  nbOpt   <- as.integer(first_argmin(res.prd$tStats[ , "AICc"]))

  tree.II     <- tree.I
  tree.II$cor <- res.prd$tStats[ , "R2cal"]


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  #  Outputs
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  res <- list(nbElt, nbAss, nbXpr,
              opt.method, opt.mean, opt.model, opt.jack, jack,
              fobs, mOccur, xpr,
              tree.I, tree.II, nbOpt,
              mCal, mPrd, mMot, mStats,
              res.prd$tCal, res.prd$tPrd, res.prd$tNbcl, res.prd$tStats)

  names(res) <- c("nbElt", "nbAss", "nbXpr",
                  "opt.method", "opt.mean", "opt.model", "opt.jack", "jack",
                  "fobs", "mOccur", "xpr",
                  "tree.I", "tree.II", "nbOpt",
                  "mCal", "mPrd", "mMotifs", "mStats",
                  "tCal", "tPrd", "tNbcl", "tStats")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# End of file                                                              ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
