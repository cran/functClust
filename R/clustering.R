#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                               CLUSTERING.R
#
#  Functions for clustering components or assembly motifs                  ####
#     on the basis of the performances of assemblages of components
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#' @include
#'  stats.R
#'  labelling.R
#'  calibrating.R
#'
NULL


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Miscellaneous functions                                                 ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Number of Stirling of second kind
#'
#' @description Take an integer and return the number
#'  of Stirling of second kind, that is the number
#'  of all possible partitions into \code{k} clusters
#'   among \code{n} components.
#'
#' @usage stirling(n)
#'
#' @param n a positive integer.
#' This corresponds to the size of set,
#' that is the number of components that belong to the set.
#'
#' @return Return a vector of integer. The \code{names} are the number
#'  \code{k} of clusters. The \code{values} are the number of partitions
#'   into \code{k} clusters possibly combinated with \code{n} components.
#'
#' @details None.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

stirling <- function(n) {

  res        <- integer(n)
  names(res) <- seq_len(n)

  for (k in seq_len(n)) {
    tmp <- 0
    for (j in 0:k) tmp <- tmp + (-1) ^ (k - j) * choose(k, j) * j ^ n
    res[k] <- tmp / factorial(k)
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#' @title Maximal number of clustering models to test
#'
#' @description Take the number of componens,
#' then compute the maximal number of clustering models to test.
#'
#' @usage nb_tests(n)
#'
#' @param n an interger, that is the number of components
#' that belong to the system.
#'
#' @details None.
#'
#' @return Return the maximal number of clustering models to test.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

nb_tests <- function(n) {

  return( n * (n + 1) * (2 * n + 1) / 12 + n * (n + 1) / 4 - n)
}




#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Component clustering                                                    ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Message for following functional clustering computation
#'
#' @description Print a message indicating
#'   where the functional clustering computation is done.
#'
#' @usage
#' notify_fclust(where, all)
#'
#' @param where an integer, that indicates the last tree level done.
#'
#' @param all an integer, that indicates the number of tree levels to do.
#'
#' @details None.
#'
#' @return Noting. It is a procedure.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

notify_fclust <- function(where, all) {

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Return TRUE if the element is odd
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  is_odd <- function(x) {
    return( ifelse(x %% 2 == 1, TRUE, FALSE) )
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  if (getOption("verbose") == TRUE) {

    main  <- "Functional clustering of components at level"
    rline <- 76 - nchar(paste(main, paste0(all, "/", all, ":"), sep = " "))

    if (all <= rline) {
      nbzero  <- where
      nbpoint <- all - nbzero
    } else {
      nbzero  <- floor(rline * where / all)
      if ( (is_odd(where) == TRUE) & (nbzero < rline) ) nbzero <- nbzero + 1
      nbpoint <- rline - nbzero
    }

    cat("\r", main, " ",
        paste0(where, "/", all, ": "),
        paste0(paste0(rep("o", nbzero),  collapse = ""),
               paste0(rep(".", nbpoint), collapse = "")), sep = "")

    if (where == all) cat("\n")
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title  Residual Sum of Squares of a given clustering model
#'
#' @description Take performance, occurrence matrix
#' and a clustering model,
#' then compute the corresponding Residual Sum of Squares.
#'
#' @usage
#' rss_clustering(fobs, affectElt, mOccur, xpr, opt.mean, opt.model )
#'
#' @inheritParams fit_ftree
#'
#' @details None.
#'
#' @return Return the Residual Sum of Squares
#' of a given clustering model.
#' Its value is computed according to \code{opt.mean} and \code{opt.model}.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

rss_clustering <- function(fobs, affectElt, mOccur, xpr,
                           opt.mean, opt.model) {

  prd <- calibrate_byminrss(fobs,
                            affect_motifs(affectElt, mOccur),
                            mOccur, xpr,
                            opt.mean, opt.model)

  setXpr <- unique(names(xpr))
  tmp <- wmp <- numeric(length(setXpr))

  for (ipr in seq_along(setXpr)) {
    index    <- which(names(xpr) == setXpr[ipr])
    wmp[ipr] <- unique(xpr[index])
    tmp[ipr] <- rss(fobs[index], prd[index])
  }

  res <- sum(wmp * tmp)

  return(res)
}




#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title  Total Residual Sum of Squares of observed performances
#'
#' @description Take performance,
#' then compute the total Residual Sum of Squares.
#'
#' @usage
#' rss_total(fobs, xpr, opt.mean)
#'
#' @inheritParams fit_ftree
#'
#' @details None.
#'
#' @return Return the Total Residual Sum of Squares
#' of observed performances of assemblages.
#' Its value is computed according to \code{opt.mean}.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

rss_total <- function(fobs, xpr, opt.mean) {

  setXpr <- unique(names(xpr))
  tmp <- wmp <- numeric(length(setXpr))

  for (ipr in seq_along(setXpr)) {
    index    <- which(names(xpr) == setXpr[ipr])
    wmp[ipr] <- unique(xpr[index])
    tmp[ipr] <- rss(fobs[index], mean_fct(fobs[index], opt.mean))
  }

  res <- sum(wmp * tmp)

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Check a matrix of component affectation to functional groups
#'
#' @description Take a matrix of component affectation,
#' then check that the numero of each affectation
#' increases regularly,
#' from \code{1} until to the number of components,
#' with the number of functional groups.
#'
#' @usage
#' check_ftree(mat)
#'
#' @param mat a square-matrix of integer.
#'
#' @details None.
#'
#' @return Return a regular matrix.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

check_ftree <- function(mat) {

  nbElt <- dim(mat)[1]

  for (elt in 1:nbElt) {

    index  <- max(mat[elt, ])
    index1 <- which(mat == index, arr.ind = TRUE)
    index2 <- which(mat == elt, arr.ind = TRUE)

    for (i in 1:dim(index1)[1])
      mat[index1[i, "row"], index1[i, "col"]] <- elt
    for (i in 1:dim(index2)[1])
      mat[index2[i, "row"], index2[i, "col"]] <- index
  }

  return(mat)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Hierarchical divisive clustering of components
#'
#' @description We proceed by division,
#' varying the number of functional groups of components
#' from 1 to the number of components.
#' All components are initially regrouped
#' into a single, large, trivial functional group.
#' At each step, one of the functional groups is split
#' into two new functional groups: the new functional groups selected are
#' those that minimize the Residual Sum of Squares of the clustering.
#' The process stops when each component is isolated in a singleton,
#' that is when there are so many clsyters as components.
#' As a whole, the process generates a hierarchical divisive tree
#' of component clustering, whose RSS decreases monotonically
#'  with the number of functional groups.
#'
#' @details At each hierarchical level of the divisive tree,
#' the division of the existing functional groups
#' into new functional groups proceeds as follows.
#' Each existing functional group is successively split
#' into two new functional groups. To do that, each component
#' of the functional group is isolated into a singleton:
#' the singleton-component that minimizes RSS is selected
#' as the nucleus of the new functional group.
#' Each of the other components belonging to the existing functional group
#' is successively moved towards the new functional group:
#' the component clustering that minimizes RSS is kept.
#' Moving component into the new functional group continues
#' as long as the new component clustering decreases RSS.
#'
#' @usage
#' divisive_ftree(fobs, mOccur, xpr, opt.mean, opt.model, opt.nbMax)
#'
#' @inheritParams fit_ftree
#'
#' @return Return an object "tree",
#' that is a list containing
#' \emph{(i)} \code{tree$aff}: an integer square-matrix of
#' component affectation to functional groups,
#' \emph{(ii)} \code{tree$cor}: a numeric vector of
#' coefficient of determination.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

divisive_ftree <- function(fobs, mOccur, xpr,
                           opt.mean, opt.model, opt.nbMax) {

  nbElt <- dim(mOccur)[2]
  res   <- matrix(Inf, nrow = nbElt, ncol = nbElt,
                  dimnames = list(seq_len(nbElt), colnames(mOccur)))

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # First line : Trunk
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  nbclElt   <- 1L
  affectElt <- as.integer(rep(1L, nbElt))
  RES.aff   <- matrix(affectElt, nrow = nbElt, ncol = nbElt, byrow = TRUE,
                      dimnames = list(seq_len(nbElt), colnames(mOccur)))
  storage.mode(RES.aff) <- "integer"

  RES.rss    <- RES.ord <- numeric(nbElt)
  RES.rss[1] <- oldRes <- rss_clustering(fobs, affectElt, mOccur,
                                         xpr, opt.mean, opt.model)

  notify_fclust(nbclElt, nbElt)

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Following lines : Leaves
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  oldRes.affectElt <- affectElt

  if (nbElt > 1L) for (nbclElt in 2:nbElt) {

    # nbclElt <- nbclElt + 1L
    oldRes <- res[ , ] <- Inf

    # One tests split of each leaves of tree
    for (clElt in 1:(nbclElt - 1L)) {
      affectElt <- RES.aff[nbclElt - 1L, ]
      setElt    <- which(affectElt == clElt)

      test2     <- (length(setElt) > 1)
      while (test2) {
        # One tests components one-after-one
        for (elt in setElt) {
          oldAffect       <- affectElt
          affectElt[elt]  <- nbclElt
          res[clElt, elt] <- rss_clustering(fobs, affectElt, mOccur,
                                            xpr, opt.mean, opt.model)
          affectElt       <- oldAffect
        }

        # Decision: one keeps the best local fit
        minRes <- min(res, na.rm = TRUE)
        if (minRes < oldRes) {
          coord                  <- which(res == minRes, arr.ind = TRUE)
          affectElt[coord[1, 2]] <- nbclElt
          oldRes                 <- minRes
          oldRes.affectElt       <- affectElt
        } else {
          test2 <- FALSE
        }
      }
      # END of WHILE on test
    }
    # END of LOOP on LEAVE

    #     with calibrate_byminrss(), oldRes is always < Inf,
    # but with calibrate_byminrss(), oldRes is == Inf when all Pred == NA
    if (minRes < Inf) {
      affectElt          <- oldRes.affectElt
      RES.aff[nbclElt, ] <- affectElt
      RES.rss[nbclElt]   <- rss_clustering(fobs, affectElt, mOccur,
                                           xpr, opt.mean, opt.model)
    } else {
      nbclElt <- nbclElt - 1L
    }

    # to avoid unuseful time-consuming computations...
    if ( (minRes == Inf) | (RES.rss[nbclElt] == 0) |
         (RES.rss[nbclElt] == RES.rss[nbclElt - 1L]) |
         (nbclElt > opt.nbMax) ) {

      tmp   <- table(RES.aff[nbclElt, ])
      index <- as.integer(names(tmp))[tmp > 1]
      while (length(index) > 0) {
        nbclElt <- nbclElt + 1L
        affectElt[ which(affectElt == index[1])[1] ] <- nbclElt

        RES.aff[nbclElt, ] <- affectElt
        RES.rss[nbclElt]   <- RES.rss[nbclElt - 1L]

        tmp     <- table(RES.aff[nbclElt, ])
        index   <- as.integer(names(tmp))[tmp > 1]
      }
      break
    }

    notify_fclust(nbclElt, nbElt)
  }
  # END of the LOOP

  notify_fclust(nbElt, nbElt)

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Build and format the tree
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  RES.rss <- 1.0 - RES.rss / rss_total(fobs, xpr, opt.mean)
  names(RES.rss) <- seq_len(nbElt)

  tree        <- list(RES.aff, RES.rss)
  names(tree) <- c("aff", "cor")

  return(tree)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Hierarchical agglomerative clustering of components
#'
#' @description We proceed by grouping,
#' varying the number of functional groups of components
#' from the number of components until to 1.
#' All components are initially dispersed
#' into a singleton, as many singletons as components.
#' At each step, one of the functional groups is grouped
#' with another functional group: the new functional groups selected are
#' those that minimize the Residual Sum of Squares of the clustering.
#' The process stops when all components are grouped
#' into a large, unique functional group.
#' As a whole, the process generates a hierarchical aggloimerative tree
#' of component clustering, whose RSS decreases monotonically
#' with the number of functional groups.
#'
#' @details At each hierarchical level of the agglomerative tree,
#' the clustering of the existing functional groups
#' into new functional groups proceeds as follows.
#' Each existing functional group is successively grouped
#' with other functional groups.
#' The component clustering that minimizes RSS is kept.
#'
#' @usage
#' agglomerative_ftree(fobs, mOccur, xpr, opt.mean, opt.model)
#'
#' @inheritParams fit_ftree
#'
#' @return Return an object "tree",
#' that is a list containing
#' \emph{(i)} \code{tree$aff}: an integer square-matrix of
#' component affectation to functional groups,
#' \emph{(ii)} \code{tree$cor}: a numeric vector of
#' coefficient of determination.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

agglomerative_ftree <- function(fobs, mOccur, xpr, opt.mean, opt.model) {

  nbElt <- dim(mOccur)[2]
  res   <- matrix(Inf, nrow = nbElt, ncol = nbElt,
                  dimnames = list(seq_len(nbElt), colnames(mOccur)))

  nbclElt <- nbElt
  RES.aff <- matrix(as.integer(0), nrow = nbElt, ncol = nbElt,
                    dimnames = list(seq_len(nbElt), colnames(mOccur)))
  storage.mode(RES.aff) <- "integer"
  RES.aff[nbElt, ] <- affectElt <- as.integer(seq(1, nbElt))
  names(affectElt) <- colnames(mOccur)

  RES.rss <- numeric(nbElt)
  RES.rss[nbElt] <- oldRes <- rss_clustering(fobs, affectElt, mOccur,
                                             xpr, opt.mean, opt.model)

  notify_fclust(nbclElt, nbElt)

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Following lines : branchs and trunk
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  # Test adding on a component to each assembly class
  if (nbElt > 1) for (nbclElt in (nbElt - 1):1) {

    # nbclElt  <- nbclElt - 1
    res[ , ] <- 0.0
    setElt   <- sort(unique(affectElt))

    set1 <- setdiff(setElt, max(setElt))
    for (elt1 in set1) {
      set2 <- setElt[setElt > elt1]

      # Test the components one-after-one
      for (elt2 in set2) {
        oldAffect       <- affectElt
        affectElt[affectElt == elt2] <- elt1
        res[elt1, elt2] <- rss_clustering(fobs, affectElt, mOccur,
                                          xpr, opt.mean, opt.model)
        affectElt       <- oldAffect
      }
    }

    # Decision: one keeps the best fit
    affectElt <- RES.aff[nbclElt + 1, ]

    maxRes <- min(res[res != 0])
    coord  <- which(res == maxRes, arr.ind = TRUE)
#    if (dim(coord)[1] > 1) cat(" ")

    affectElt[affectElt == coord[1, 2]] <- coord[1, 1]
    RES.aff[nbclElt, ] <- affectElt
    RES.rss[nbclElt]   <- rss_clustering(fobs, affectElt, mOccur,
                                         xpr, opt.mean, opt.model)
    # print to follow the computation
    notify_fclust(nbclElt, nbElt)
  }
  # END of LOOP

  notify_fclust(nbElt, nbElt)

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Build and format the tree
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  RES.aff <- check_ftree(RES.aff)

  RES.rss        <- 1.0 - RES.rss
  names(RES.rss) <- seq_len(nbElt)

  tree        <- list(RES.aff, RES.rss)
  names(tree) <- c("aff", "cor")

  return(tree)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Hierarchical clustering of components
#' from an \emph{a priori} component clustering
#'
#' @description An \emph{a priori} clustering of components
#' is given, therefore coerced by the user.
#' We generate a hierarchical tree coerced
#' by the \emph{a priori} component clustering.
#' We proceeds in two steps:
#' \emph{(i)} by division from the coerced tree-level towards the leaves,
#' \emph{(i)} by grouping from the coerced tree-level towards the trunk. \cr
#'
#' @details
#' \code{"divisive"}: We proceed by division,
#' varying the number of functional groups of components
#' from 1 to the number of components.
#' All components are initially regrouped
#' into a single, large, trivial functional group.
#' At each step, one of the functional groups is split
#' into two new functional groups: the new functional groups selected are
#' those that minimize the Residual Sum of Squares of the clustering.
#' The process stops when each component is isolated in a singleton,
#' that is when there are so many clsyters as components.
#' As a whole, the process generates a hierarchical divisive tree
#' of component clustering, whose RSS decreases monotonically
#' with the number of functional groups. \cr
#'
#' At each hierarchical level of the divisive tree,
#' the division of the existing functional groups
#' into new functional groups proceeds as follows.
#' Each existing functional group is successively split
#' into two new functional groups. To do that, each component
#' of the functional group is isolated into a singleton:
#' the singleton-component that minimizes RSS is selected
#' as the nucleus of the new functional group.
#' Each of the other components belonging to the existing functional group
#' is successively moved towards the new functional group:
#' the component clustering that minimizes RSS is kept.
#' Moving component into the new functional group continues
#' as long as the new component clustering decreases RSS.
#'
#' \code{"agglomerative"}: We proceed by grouping,
#' varying the number of functional groups of components
#' from the number of components until to 1.
#' All components are initially dispersed
#' into a singleton, as many singletons as components.
#' At each step, one of the functional groups is grouped
#' with another functional group: the new functional groups selected are
#' those that minimize the Residual Sum of Squares of the clustering.
#' The process stops when all components are grouped
#' into a large, unique functional group.
#' As a whole, the process generates a hierarchical aggloimerative tree
#' of component clustering, whose RSS decreases monotonically
#' with the number of functional groups. \cr
#'
#' At each hierarchical level of the agglomerative tree,
#' the clustering of the existing functional groups
#' into new functional groups proceeds as follows.
#' Each existing functional group is successively grouped
#' with other functional groups.
#' The component clustering that minimizes RSS is kept.
#'
#' @usage
#' complete_ftree(fobs, mOccur, xpr, affectElt, opt.mean, opt.model,
#'               opt.nbMax = dim(mOccur)[2] )
#'
#' @inheritParams fit_ftree
#'
#' @return Return an object "tree",
#' that is a list containing
#' \emph{(i)} \code{tree$aff}: an integer square-matrix of
#' component affectation to functional groups,
#' \emph{(ii)} \code{tree$cor}: a numeric vector of
#' coefficient of determination.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

complete_ftree <- function(fobs, mOccur, xpr, affectElt,
                           opt.mean, opt.model, opt.nbMax = dim(mOccur)[2]) {

  nbElt <- dim(mOccur)[2]
  res   <- matrix(Inf, nrow = nbElt, ncol = nbElt,
                  dimnames = list(seq_len(nbElt), colnames(mOccur)))

  nbclElt <- length(unique(affectElt))
  RES.aff <- matrix(as.integer(0), nrow = nbElt, ncol = nbElt,
                    dimnames = list(seq_len(nbElt), colnames(mOccur)))
  storage.mode(RES.aff) <- "integer"
  RES.aff[nbclElt, ] <- affectElt

  RES.rss <- numeric(nbElt)
  RES.rss[nbclElt] <- oldRes <- rss_clustering(fobs, affectElt, mOccur,
                                               xpr, opt.mean, opt.model)

  notify_fclust(nbclElt, nbElt)

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Agglomerative component clustering from "a priori" tree-level towards trunk
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  #     Following lines : branchs and trunk
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  nbcl <- nbclElt

  # Test adding on a component to each assembly class
  if ( (nbElt > 1) & (nbclElt > 1) )
    for (nbcl in (nbclElt - 1):1) {

      # nbcl <- nbcl - 1
      res[ , ] <- 0.0
      setElt   <- sort(unique(affectElt))

      set1 <- setdiff(setElt, max(setElt))
      for (elt1 in set1) {
        set2 <- setElt[setElt > elt1]

        # Test the components one-after-one
        for (elt2 in set2) {
          oldAffect       <- affectElt
          affectElt[affectElt == elt2] <- elt1
          res[elt1, elt2] <- rss_clustering(fobs, affectElt, mOccur,
                                            xpr, opt.mean, opt.model)
          affectElt       <- oldAffect
        }
      }

      # Decision: one keeps the best fit
      affectElt <- RES.aff[nbcl + 1, ]

      maxRes <- min(res[res != 0])
      coord  <- which(res == maxRes, arr.ind = TRUE)

      affectElt[affectElt == coord[1, 2]] <- coord[1, 1]
      RES.aff[nbcl, ] <- affectElt
      RES.rss[nbcl]   <- rss_clustering(fobs, affectElt, mOccur,
                                        xpr, opt.mean, opt.model)
      # print to follow the computation
      notify_fclust(nbcl, nbElt)

    }   # END of LOOP

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Build and format the tree
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  RES.aff[1:nbclElt, ] <- check_ftree(RES.aff[1:nbclElt, , drop = FALSE])



  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Divisive component clustering from "a priori" tree-level towards leaves
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Following lines : Leaves
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  oldRes.affectElt <- affectElt <- RES.aff[nbclElt, ]

  if ((nbElt > 1) & (nbclElt < opt.nbMax)) {
    for (nbcl in (nbclElt + 1):nbElt) {

      # nbcl <- nbcl + 1
      oldRes <- res[ , ] <- Inf

      # One tests split of each leaves of tree
      for (clElt in 1:(nbcl - 1)) {
        affectElt <- RES.aff[nbcl - 1, ]
        setElt    <- which(affectElt == clElt)

        test2     <- (length(setElt) > 1)
        while (test2) {
          # One tests components one-after-one
          for (elt in setElt) {
            oldAffect       <- affectElt
            affectElt[elt]  <- nbcl
            res[clElt, elt] <- rss_clustering(fobs, affectElt, mOccur,
                                              xpr, opt.mean, opt.model)
            affectElt       <- oldAffect
          }

          # Decision: one keeps the best local fit
          minRes <- min(res, na.rm = TRUE)
          if (minRes < oldRes) {
            coord                 <- which(res == minRes, arr.ind = TRUE)
            affectElt[coord[1, 2]] <- nbcl
            oldRes                 <- minRes
            oldRes.affectElt       <- affectElt
          } else {
            test2 <- FALSE
          }
        }
        # END of WHILE on test
      }
      # END of LOOP on LEAVE

      #     with calibrate_byminrss(), oldRes is always < Inf,
      # but with calibrate_byminrss(), oldRes is == Inf when all Pred == NA
      if (minRes < Inf) {
        affectElt       <- oldRes.affectElt
        RES.aff[nbcl, ] <- affectElt
        RES.rss[nbcl]   <- rss_clustering(fobs, affectElt, mOccur,
                                          xpr, opt.mean, opt.model)
      } else {
        nbcl <- nbcl - 1
      }

      # to avoid unuseful time-consuming computations...
      if ( (minRes == Inf) | (RES.rss[nbcl] == 0) |
           (RES.rss[nbcl] == RES.rss[nbcl - 1]) |
           (nbcl == opt.nbMax) )  {

        tmp   <- table(RES.aff[nbcl, ])
        index <- as.integer(names(tmp))[tmp > 1]
        while (length(index) > 0) {
          nbcl <- nbcl + 1
          affectElt[ which(affectElt == index[1])[1] ] <- nbcl

          RES.aff[nbcl, ] <- affectElt
          RES.rss[nbcl]   <- RES.rss[nbcl - 1]

          tmp     <- table(RES.aff[nbcl, ])
          index   <- as.integer(names(tmp))[tmp > 1]
        }
        break
      }

      # print to follow the computation
      notify_fclust(nbcl, nbElt)

    }  # END of the LOOP

  } else {

    # to fill up the upper part of the tree
    nbcl  <- nbclElt
    tmp   <- table(RES.aff[nbcl, ])
    index <- as.integer(names(tmp))[tmp > 1]
    while (length(index) > 0) {
      nbcl <- nbcl + 1
      affectElt[ which(affectElt == index[1])[1] ] <- nbcl

      RES.aff[nbcl, ] <- affectElt
      RES.rss[nbcl]   <- RES.rss[nbcl - 1]

      tmp     <- table(RES.aff[nbcl, ])
      index   <- as.integer(names(tmp))[tmp > 1]
    }
  }

  notify_fclust(nbElt, nbElt)

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Build and format the tree
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  RES.aff <- check_ftree(RES.aff)

  RES.rss        <- 1.0 - RES.rss
  names(RES.rss) <- seq_len(nbElt)

  tree        <- list(RES.aff, RES.rss)
  names(tree) <- c("aff", "cor")

  return(tree)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Main function                                                           ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Clustering of components for the performances of assemblages
#'
#' @description The function checks the inputs,
#' then switch on different methods of component clustering
#' according to the option \code{opt.method}.
#' Three methods (\code{opt.method = c("sort", "divisive", "agglomerative")})
#' generate a hierarchical tree.
#' The last method (\code{opt.method = "cluster"})
#' generates a non-hierarchical tree of clustering.
#'
#' @usage
#' fit_ftree(fobs, mOccur,
#'          xpr = stats::setNames(rep(1, length(fobs)),rep("a", length(fobs))),
#'          affectElt  = rep(1, dim(mOccur)[2]),
#'          opt.method = "divisive",
#'          opt.mean   = "amean",
#'          opt.model  = "byelt",
#'          opt.nbMax  = dim(mOccur)[2] )
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
#' @param affectElt a vector of integers
#' of \code{length(affectElt) == dim(mOccur)[1]},
#' that is the number of components.
#' The vector contains the labels of different functional clusters
#' to which each component belongs.
#' Each functional cluster is labelled as an integer, and
#' each component must be identified by its name in \code{names(affectElt)}.
#' The number of functional clusters defined in \code{affectElt}
#' determines an \emph{a priori} level of component clustering
#' (\code{level <- length(unique(affectElt))}).\cr
#'
#' If \code{affectElt = NULL} (value by default),
#' the option \code{opt.method} must be filled out.
#' A tree is built,
#' from a unique trunk to as many leaves as components
#' by using the specified method. \cr
#'
#' If \code{affectElt} is specified,
#' the option \code{opt.method} does not need to be filled out.
#' \code{affectElt} determines an \emph{a priori}
#' level of component clustering,
#' and a tree is built:
#' \emph{(i)} by using \code{opt.method =  "divisive"}
#' from the \emph{a priori} defined level in tree towards
#' as many leaves as components;
#' \emph{(ii)} by using \code{opt.method =  "agglomerative"}
#' from the \emph{a priori} defined level in tree towards the tree trunk
#' (all components are together withi a trivial singleton).
#'
#' @param opt.method a string that specifies the method to use.
#' \code{opt.method = c("divisive", "agglomerative", "complete")}.
#' All the methods generate a hierarchical tree.
#' Each tree is complete, running from a unique trunk
#' to as many leaves as components. \cr
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
#' If \code{opt.method = "complete"}, \code{affectElt} must be specified.
#' \code{affectElt} determines a level of component clustering,
#' and a tree is built:
#' \emph{(i)} by using \code{opt.method =  "divisive"}
#' from the defined level in tree towards as many leaves as components;
#' \emph{(ii)} by using \code{opt.method =  "agglomerative"}
#' from the defined level in tree towards the trunk of tree.
#'
#' @param opt.mean a character equals to \code{"amean"} or \code{"gmean"}.
#' Switchs to arithmetic formula if \code{opt.mean = "amean"}.
#' Switchs to geometric formula if \code{opt.mean = "gmean"}. \cr
#'
#' Modelled performances are computed
#' using arithmetic mean (\code{opt.mean = "amean"})
#' or geometric mean (\code{opt.mean = "gmean"})
#' according to \code{opt.model}.
#'
#' @param opt.model a character equals to \code{"bymot"} or \code{"byelt"}.
#' Switchs to simple mean by assembly motif if \code{opt.model = "bymot"}.
#' Switchs to linear model with assembly motif if \code{opt.model = "byelt"}.
#'   \cr
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
#' as the assemblage to predict.
#' This procedure corresponds to a linear model within each assembly motif
#' based on the component occurrence in each assemblage.
#' If no assemblage contains component belonging to assemblage to predict,
#' performance is the mean performance of all assemblages
#' as in \code{opt.model = "bymot"}.
#'
#' @param opt.nbMax an integer, comprizes between 1 and nbElt,
#' that indicates the last level of hierarchical tree to compute.
#' This option is very useful to shorten computing-time
#' in the test-functions
#' \code{\link{ftest_components}}, \code{\link{ftest_assemblages}},
#' \code{\link{ftest_performances}}, \code{\link{fboot_assemblages}},
#' \code{\link{fboot_performances}} or \code{\link{ftest}}
#' where the function \code{\link{fit_ftree}} is run very numerous times.
#'
#' @details None.
#'
#' @return Return a primary tree of clustering of components
#' for their effect on the performances of component assemblages.
#'
#' @importFrom stats setNames
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

fit_ftree <- function(fobs, mOccur,
                      xpr = stats::setNames(rep(1, length(fobs)),
                                            rep("a", length(fobs))),
                      affectElt  = rep(1, dim(mOccur)[2]),
                      opt.method = "divisive",
                      opt.mean   = "amean",
                      opt.model  = "byelt",
                      opt.nbMax  = dim(mOccur)[2]    )  {

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Check the quantitative inputs
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  if (!is.matrix(mOccur)) mOccur <- as.matrix(mOccur)
  storage.mode(mOccur) <- "integer"

  tmp <- names(table(mOccur))
  if ( (sum(is.na(mOccur)) != 0) | (length(tmp) != 2) )
    stop("The matrix of occurence should be a binary matrix")

  if (((tmp[1] != "0") & (tmp[1] != "FALSE")) ||
      ((tmp[2] != "1") & (tmp[2] != "TRUE")))
    stop("The matrix of occurrence should be a binary matrix")

  nbAss <- dim(mOccur)[1]
  if ((length(fobs) != nbAss) | (sum(!is.na(fobs)) != nbAss))
    stop("length(fobs) == dim(mOccur)[1]")

  if ((length(names(fobs)) != 0) | (length(rownames(mOccur)) != 0)) {
    if (length(names(fobs)) != 0) rownames(mOccur) <- names(fobs)
  } else {
    stop("names(fobs) or rownames(mOccur) must be informed")
  }

  if (length(xpr) != length(fobs)) stop("length(xpr) == length(fobs)")

  setXpr <- unique(names(xpr))
  if (length(setXpr) > 1) for (ipr in seq_along(setXpr))
    if (length(unique(xpr[names(xpr) == setXpr[ipr]])) > 1)
      stop("xpr: all assemblages of a same experiment must have a same weight")



  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Check the options
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  all <- c("divisive", "agglomerative", "apriori")
  if (is.na(pmatch(opt.method, all)))
    stop("'opt.method' should be ", list_in_quote(all))
  opt.method <- match.arg(opt.method, all)

  if (length(affectElt) == dim(mOccur)[2]) {
    if (is.character(affectElt) == TRUE) {
      affectElt <- char_to_int(affectElt)
    } else {
      affectElt <- compact_index(affectElt)
    }
    names(affectElt) <- colnames(mOccur)
  } else {
    stop("'length(affectElt)' should be equal to 'nbElt'")
  }

  all <- c("amean", "gmean")
  if (is.na(pmatch(opt.mean, all)))
    stop("'opt.mean' should be ", list_in_quote(all))
  opt.mean <- match.arg(opt.mean, all)

  all <- c("bymot", "byelt")
  if (is.na(pmatch(opt.model, all)))
    stop("'opt.model' should be ", list_in_quote(all))
  opt.model <- match.arg(opt.model, all)

  if ( (opt.nbMax < 1) | (opt.nbMax > dim(mOccur)[2]) )
    stop("opt.nbMax is comprised between 1 and nbElt")


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Cluster the components
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  res <-
    switch(opt.method,
           divisive      = divisive_ftree(fobs, mOccur, xpr,
                                          opt.mean, opt.model, opt.nbMax),
           agglomerative = agglomerative_ftree(fobs, mOccur, xpr,
                                               opt.mean, opt.model),
           apriori       = complete_ftree(fobs, mOccur, xpr, affectElt,
                                          opt.mean, opt.model, opt.nbMax)
    )

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# End of file                                                              ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
