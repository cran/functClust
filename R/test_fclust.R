#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                         TEST_CLUSTERING.R
#
#  Functions for testing robustess of functional clustering                ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#' @include
#'  clustering.R
#'  validating.R
#'  plot_fclust.R
#'
NULL



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Graphical tools                                                         ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Plot a point with x- and y-error bars
#'
#' @description Take two vectors,
#'   compute their mean and standard deviation,
#'   then plot the mean point with error bars along x- and y-axis.
#'
#' @usage
#' points_sd(x, y, opt = c("mean", "mean"), fig, col, bg, cex = 2, name = "")
#'
#' @param x,y two numeric vectors, x in abscissa and y in coordinate.
#'
#' @param opt two strings, equal to \code{"mean"} or \code{"median"}.
#' Indicate if the point should be plotted at \code{mean(x), mean(y)},
#' \code{mean(x), median(y)}, \code{median(x), mean(y)}
#' or \code{median(x), median(y)}.
#'
#' @param fig,col two integers, indicating the symbol and the colour to use.
#'
#' @param bg a string, indicating the colour background to use.
#'
#' @param cex a numeric, indicating the size of symbol to plot.
#'
#' @param name a string, indicating the label to add near the point.
#'
#' @details None.
#'
#' @importFrom graphics arrows points text
#'
#' @return Nothing. It is a procedure.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

points_sd <- function(x, y, opt = c("mean", "mean"),
                      fig, col, bg = "white", cex = 2, name = "") {

  mx <- ifelse(opt[1] == "mean", mean(x), median(x))
  my <- ifelse(opt[2] == "mean", mean(y), median(y))

  # plot the error-bars
  if (length(x) > 1) {

    dx <- asd(x)
    if (dx > fepsilon())
      graphics::arrows(x0 = mx - dx, y0 = my,
                       x1 = mx + dx, y1 = my,
                       length = 0.1, angle = 90, code = 3,
                       col = col, lwd = 1, lty = "solid")

    dy <- asd(y)
    if (dy > fepsilon())
      graphics::arrows(x0 = mx, y0 = my - dy,
                       x1 = mx, y1 = my + dy,
                       length = 0.1, angle = 90, code = 3,
                       col = col, lwd = 1, lty = "solid")
  }

  # plot the mean point
  graphics::points(x = mx, y = my,
                   type = "p", pch = fig, col = col,
                   bg = bg, cex = cex)
  graphics::text(x = mx, y = my,
                 labels = name, pos = 4, col = col, bg = "white")
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title sort a matrix
#'
#' @description Take a matrix, then sorted the colums by increasing values
#' of 1st, 2nd, 3rd,..    n lines of matrix.
#'
#' @usage sort_matrix(mat)
#'
#' @param mat a numeric matrix
#'
#' @details None.
#'
#' @return A vector of indices for properly sort the columns of matrix
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

sort_matrix <- function(mat) {

  if (is.matrix(mat)) {

    lst <- vector(mode = "list", length(dim(mat)[2]))
    for (i in seq_len(dim(mat)[2])) lst[[i]] <- mat[i, ]
    res <- do.call(what = order, args = lst)

  } else {

    res <- 1
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Mathematical mark of an ordered string vector
#'
#' @description Take two vectors: a ordered vector of numeric
#' and the corresponding of strings. Then compare each term
#' of the ordered suit of numeric and add "=", ">" or "<"
#' at the end of each string.
#'
#' @usage mark_string(vNum, vStr)
#'
#' @param vNum an ordered (increasing or decreasing) vector of numerics.
#'
#' @param vStr the corresponding ordered vector of strings.
#'
#' @details None.
#'
#' @return The vector of strings with added mathematical symbol.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

mark_string <- function(vNum, vStr) {

  if (length(vStr) > 1)
    for (i in seq_len(length(vStr) - 1)) {
      if (vNum[i] == vNum[i + 1]) {
        vStr[i] <- paste0(vStr[i], " = ")
      } else {
        if (vNum[i] > vNum[i + 1]) {
          vStr[i] <- paste0(vStr[i], " > ")
        } else {
          vStr[i] <- paste0(vStr[i], " < ")
        }
      }
    }

  return(vStr)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Test significant components                                             ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'a
#' @title Evaluate the weight of each component on functional clustering
#'
#' @description Evaluate by cross-validation (leave-one-out)
#'  the effect induced by each component
#'  on the result of a functional clustering.
#'
#' @usage
#' ftest_components(fres, opt.nbMax = fres$nbOpt,
#'                  opt.R2 = FALSE, opt.plot = FALSE)
#'
#' @inheritParams ftest
#'
#' @details Each component of the interactive system in consideration
#'  is successively removed from the dataset,
#'  the remaining components are functionally clustered,
#'  then indices of distance between clustering trees
#'  with and without the component are computed.
#'  The components can then be hierarchised depending on
#'  the distance induced by their removing from dataset.
#'  The used distance criteria are :
#'   "Czekanowski_Dice", "Folkes_Mallows", "Jaccard", "Kulczynski",
#'   "Precision", "Rand", "Recall", "Rogers_Tanimoto", "Russel_Rao",
#'   "Sokal_Sneath1" and "Sokal_Sneath2" index.
#'   For more informations, see the notice of R-package \code{clusterCrit}.
#'   The test is time-consuming.
#'
#' @references
#' Package "clusterCrit": Clustering Indices,
#'   by Bernard Desgraupes (University of Paris Ouest - Lab Modal'X)
#'
#' @importFrom clusterCrit getCriteriaNames extCriteria
#'
#' @importFrom graphics text
#'
#' @return a list containing a matrix by clustering index.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

ftest_components <- function(fres, opt.nbMax = fres$nbOpt,
                             opt.R2 = FALSE, opt.plot = FALSE) {


  # organize the inputs

  index      <- which(names(fres$xpr) == names(fres$xpr)[1])
  fobsBase   <- fres$fobs[index]
  mOccurBase <- fres$mOccur[index, ]

  setElt        <- seq_len(fres$nbElt)
  names(setElt) <- colnames(fres$tree.I$aff)

  lAffectElt <- vector(mode = "list", length = fres$nbElt)
  for (nb in seq_len(fres$nbElt))
    lAffectElt[[nb]] <- cut_ftree(fres$tree.I, nb)


  # plot the reference tree if opt.plot = TRUE

  if ((getOption("verbose") == TRUE) && (opt.plot == TRUE)) {
    cols <- fcolours()[shift_affectElt(cut_ftree(fres$tree.II, fres$nbOpt))]
    plot_ftree(fres$tree.II, cols = cols)
    graphics::text(x = 1, y = 0, labels = "Reference tree", pos = 4)
  }


  # Begin the loop on components

  indCrit   <- c(1:2,4:5,8:14)
  critNames <- clusterCrit::getCriteriaNames(FALSE)[indCrit]    #external crits
  vCriteria <- numeric(length(critNames))
  names(vCriteria) <- c(critNames)

  mCrit <- array(0, dim = c(fres$nbElt, fres$nbElt, length(critNames)))
  dimnames(mCrit) <- list(seq_len(fres$nbElt),
                          colnames(fres$mOccur), critNames)

  if (opt.R2 == TRUE) {
    mStats <- array(0, dim = c(fres$nbElt, fres$nbElt, 2))
    dimnames(mStats) <- list(seq_len(fres$nbElt),
                             colnames(fres$mOccur), c("R2", "E"))
  }

  ind.tree  <- sort_ftree(fres$tree.I, index.return = TRUE)$ix

#  for (elt in seq_len(fres$nbElt)) {
  for (elt in ind.tree) {

    # compute the primary tree

    element <- setElt[elt]

    main <- paste0("The component '", names(element), "' is removed")
    if (getOption("verbose") == TRUE) cat(main, "\n")

    indElt <- setdiff(setElt, element)
    tmp    <- check_repeat(fobsBase, mOccurBase[ , indElt])
    index  <- which(names(fres$fobs) %in% rownames(tmp$mOccur))

    tree.I <- fit_ftree(fres$fobs[index], fres$mOccur[index, indElt],
                        fres$xpr[index],  fres$affectElt[indElt],
                        fres$opt.method,  fres$opt.mean, fres$opt.model,
                        opt.nbMax )

    for (nb in seq_len(fres$nbElt - 1)) {
      vCriteria[] <- as.vector(unlist(clusterCrit::extCriteria(
        as.integer(cut_ftree(tree.I, nb)),
        as.integer(lAffectElt[[nb]][indElt]),
        crit = "all")))[indCrit]
      mCrit[nb, element, critNames] <- vCriteria[critNames]
    }

    # compute the secondary tree

    if (opt.R2 == TRUE) {
      res <- validate_ftree(tree.I,
                            fres$fobs[index], fres$mOccur[index, indElt],
                            fres$xpr[index],
                            fres$opt.method,  fres$opt.mean, fres$opt.model,
                            fres$opt.jack, fres$jack,
                            opt.nbMax )

      mStats[seq_len(fres$nbElt - 1), element, "R2"] <- res$tStats[ , "R2cal"]
      mStats[seq_len(fres$nbElt - 1), element, "E"]  <- res$tStats[ , "R2prd"]
    }

    # plot the resulting tree if opt.plot = TRUE

    if ((getOption("verbose") == TRUE) && (opt.plot == TRUE)) {
      if (opt.R2 == TRUE) {
        plot_ftree(res$tree.II, cols = cols[indElt])
      } else {
        plot_ftree(tree.I, cols = cols[indElt])
      }
      graphics::text(x = 1, y = 0, labels = main, pos = 4)
    }
  }

  coord <- which(mCrit == "NaN", arr.ind = TRUE)
  mCrit[coord] <- ifelse(amean(mCrit[fres$nbElt, , ]) < 0.1, 0, 1)


  # Organize the results

  rtest <- vector(mode = "list", length = length(critNames))
  names(rtest) <- critNames
  for (ind in seq_along(critNames))
    rtest[[ind]] <- as.matrix(mCrit[ , , critNames[ind]])

  if (opt.R2 == TRUE) {
    mStats[fres$nbElt, , "R2"] <- mStats[fres$nbElt - 1, , "R2"]
    mStats[fres$nbElt, , "E"]  <- mStats[fres$nbElt - 1, , "E"]
    rtest$R2 <- as.matrix(mStats[ , , "R2"])
    rtest$E  <- as.matrix(mStats[ , , "E"])
  }


  return(rtest)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Plot test of the weight of each component on functional clustering
#'
#' @description Evaluate by cross-validation (leave-one-out)
#'  the effect induced by leaving out each component
#'  on the result of a functional clustering.
#'
#' @usage
#' ftest_plot_components(fres, rtest, main = "Title", opt.crit = "Jaccard",
#'                       opt.comp = list("sorted.tree"))
#'
#' @inheritParams ftest_plot
#'
#' @details None.
#'
#' @references
#' Package "clusterCrit": Clustering Indices,
#'   by Bernard Desgraupes (University of Paris Ouest - Lab Modal'X)
#'
#' @importFrom clusterCrit getCriteriaNames extCriteria
#'
#' @importFrom graphics plot lines points axis
#'
#' @return a list containing a matrix by clustering index.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

ftest_plot_components <- function(fres, rtest,
                                  main = "Title", opt.crit = "Jaccard",
                                  opt.comp = list("sorted.tree")) {


  affectElt <- cut_ftree(fres$tree.II, fres$nbOpt)
  setNoms   <- unique(name_clusters(affectElt))
  setClu    <- unique(shift_affectElt(affectElt))

  cols      <- extend_vector(fcolours(), fres$nbElt)[shift_affectElt(affectElt)]
  pchs      <- extend_vector(fsymbols(), fres$nbElt)[shift_affectElt(affectElt)]
  setCol    <- unique(cols)
  setPch    <- unique(pchs)

  indClu    <- seq_len(fres$nbOpt)
  setCrit   <- names(rtest)
  indCrit   <- which(setCrit %in% opt.crit)


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Check options
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  all.comp <- list("all.together", "fgroups.together", "comps.together",
                   "fgroups.byfg", "comps.byfg", "sorted.tree", "sorted.leg",
                   "all")
  names(all.comp) <- all.comp

  opt.comp <- check_plot_options(opt.comp, all.comp, main = "opt.components",
                                 affectElt = affectElt)

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  affectElt <- shift_affectElt(affectElt)

  crit <- 5
  for (crit in indCrit) {

    mCrit    <- as.matrix(rtest[[crit]])
    mainCrit <- paste(main, "comp", setCrit[crit], sep = " / ")

    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    # plot general median
    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    if ("all.together" %in% opt.comp) {

      mainLoc <- paste(mainCrit, "general median index", sep = " / ")

      graphics::plot(x = mCrit[indClu, 1], xlim = c(0, 1),
                     y = fres$tree.II$cor[indClu], ylim = c(0, 1),
                     xlab = paste0(setCrit[crit], " index"), ylab = "R2",
                     main = mainLoc,
                     type = "n", tck = 0.02, las = 2)

      for (lev in seq_len(fres$nbOpt - 1))
        graphics::lines(x = c(median(mCrit[lev, ]), median(mCrit[lev + 1, ])),
                        y = c(fres$tree.II$cor[lev],
                              fres$tree.II$cor[lev + 1]),
                        lty = "dashed", col = fcolours()[1])

      for (lev in seq_len(fres$nbOpt))
        points_sd(x = mCrit[lev, ], y = fres$tree.II$cor[lev],
                  opt = c("median", "mean"),
                  fig = fsymbols()[1], col = fcolours()[1], bg = fcolours()[1])
    }


    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    # plot median component clusters, all together
    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    if ("fgroups.together" %in% opt.comp) {

      mainLoc <- paste(mainCrit, "median index of fgroups", sep = " / ")

      graphics::plot(x = mCrit[indClu, 1], xlim = c(0, 1),
                     y = fres$tree.II$cor[indClu], ylim = c(0, 1),
                     xlab = paste0(setCrit[crit], " index"), ylab = "R2",
                     main = mainLoc,
                     type = "n", tck = 0.02, las = 2)

      for (clu in setClu) {

        setElt <- which(affectElt == setClu[clu])
        vtmp   <- apply(mCrit[indClu, setElt, drop = FALSE],
                        MARGIN = 1, FUN = median)

        graphics::points(x = vtmp,
                         y = fres$tree.II$cor[indClu],
                         type = "b", lty = "solid",
                         pch = setPch[clu], col = setCol[clu],
                         bg = setCol[clu], cex = 2)

        graphics::axis(labels = setNoms[clu],
                       tick = FALSE, las = 1, side = 3,
                       at  = vtmp[fres$nbOpt],
                       pos = fres$tree.II$cor[fres$nbOpt] + 0.02)
      }
    }


    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    # plot mean component clusters, functional group by functional group
    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    if ("fgroups.byfg" %in% opt.comp) {

      for (clu in setClu) {

        mainLoc  <- paste(mainCrit,
                          paste0("median index for fgroup '",
                                 setNoms[clu], "'"),
                          sep = " / ")
        setElt <- which(affectElt == setClu[clu])

        graphics::plot(x = mCrit[indClu, 1, drop = FALSE], xlim = c(0, 1),
                       y = fres$tree.II$cor[indClu], ylim = c(0, 1),
                       xlab = paste0(setCrit[crit], " index"), ylab = "R2",
                       main = mainLoc,
                       type = "n", tck = 0.02, las = 2)

        for (lev in seq_len(fres$nbOpt - 1))
          graphics::lines(x = c(median(mCrit[lev, setElt]),
                                median(mCrit[lev + 1, setElt])),
                          y = c(fres$tree.II$cor[lev],
                                fres$tree.II$cor[lev + 1]),
                          lty = "dashed", col = setCol[clu])

        for (lev in seq_len(fres$nbOpt))
          points_sd(x = mCrit[lev, setElt], y = fres$tree.II$cor[lev],
                    opt = c("median", "mean"),
                    fig = setPch[clu], col = setCol[clu], bg = setCol[clu],
                    cex = 2)

        graphics::axis(labels = setNoms[clu],
                       tick = FALSE, las = 1, side = 3,
                       at  = median(mCrit[fres$nbOpt, setElt]),
                       pos = fres$tree.II$cor[fres$nbOpt] + 0.02)
      }
    }


    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    # plot a graph by component, functional group by functional group
    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    if ("comps.byfg" %in% opt.comp) {

      for (clu in setClu) {

        mainLoc <- paste(mainCrit,
                         paste0("component indices within fgroup '",
                                setNoms[clu], "'"),
                         sep = " / ")
        setElt <- which(affectElt == setClu[clu])

        graphics::plot(x = mCrit[indClu, 1, drop = FALSE], xlim = c(0, 1),
                       y = fres$tree.II$cor[indClu], ylim = c(0, 1),
                       xlab = paste0(setCrit[crit], " index"), ylab = "R2",
                       main = mainLoc,
                       type = "n", tck = 0.02, las = 2)

        for (elt in seq_along(setElt)) {

          for (lev in seq_len(fres$nbOpt - 1))
            graphics::lines(x = c(mCrit[lev, setElt[elt]],
                                  mCrit[lev + 1, setElt[elt]]),
                            y = c(fres$tree.II$cor[lev],
                                  fres$tree.II$cor[lev + 1]),
                            lty = "dashed", col = setCol[clu])

          for (lev in seq_len(fres$nbOpt))
            points_sd(x = mCrit[lev, setElt[elt]],  y = fres$tree.II$cor[lev],
                      opt = c("median", "mean"),
                      fig = extend_vector(fsymbols(), length(setElt))[elt],
                      col = setCol[clu],
                      cex = 2)

          graphics::axis(labels = colnames(fres$mOccur)[setElt[elt]],
                         tick = FALSE, las = 2, side = 3,
                         at  = mCrit[fres$nbOpt, setElt[elt]],
                         pos = fres$tree.II$cor[fres$nbOpt] + 0.02)
        }
      }
    }


    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    # plot a graph by component
    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    if ("comps.together" %in% opt.comp) {

      mainLoc <- paste(mainCrit, "all component indices", sep = " / ")

      graphics::plot(x = mCrit[indClu, 1], xlim = c(0, 1),
                     y = fres$tree.II$cor[indClu], ylim = c(0, 1),
                     xlab = paste0(setCrit[crit], " index"), ylab = "R2",
                     main = mainLoc,
                     type = "n", tck = 0.02, las = 2)

      for (clu in seq_along(setClu)) {

        setElt <- which(affectElt == setClu[clu])

        for (elt in seq_along(setElt)) {

          graphics::points(x = mCrit[indClu, setElt[elt]],
                           y = fres$tree.II$cor[indClu],
                           type = "b", lty = "solid",
                           col = cols[setElt[elt]],
                           pch = extend_vector(fsymbols(),
                                               length(setElt))[elt],
                           cex = 2)

          graphics::axis(labels = colnames(fres$mOccur)[setElt[elt]],
                         tick = FALSE, las = 2, side = 3,
                         at  = mCrit[fres$nbOpt, setElt[elt]],
                         pos = fres$tree.II$cor[fres$nbOpt] + 0.02)
        }
      }
    }


    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    # plot a tree with components sorted by importance and/or the legend
    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    if ( ("sorted.tree" %in% opt.comp) | ("sorted.leg" %in% opt.comp) )  {

      stmp  <- sort_ftree(fres$tree.II, index.return = TRUE)

      stree <- ttree <- stmp$x
      for (i in seq_len(fres$nbElt)) stree[i, ] <- shift_affectElt(stree[i, ])

      sCrit <- mCrit[ , stmp$ix]

      sindex  <- seq_len(fres$nbElt)
      for (clu in 1:fres$nbOpt) {

        index <- which(stree[fres$nbOpt, ] == clu)

        lCrit <- list()
        for (i in 1:fres$nbOpt) lCrit[[i]] <- sCrit[i, index]
        ssindex <- do.call(what = base::order, args = lCrit)

        sindex[index] <- sindex[index][ssindex]
      }

      sCrit <- sCrit[ , sindex]
      ttree <- ttree[ , sindex]

      tree.III <- list(ttree, fres$tree.II$cor)
      names(tree.III) <- c("aff", "cor")

      saffectElt <- cut_ftree(tree.III, fres$nbOpt)

      # "sorted.tree"

      if ("sorted.tree" %in% opt.comp) {

        mainLoc <- paste(mainCrit, "sorted component tree", sep = " / ")

        plot_ftree(tree.III,
                   cols = fcolours()[shift_affectElt(saffectElt)],
                   main = mainLoc, opt.sort = FALSE)
        graphics::text(x = 1, y = 0, labels = "Sorted reference tree", pos = 4)
      }

      # "sorted.leg"

      if ("sorted.leg" %in% opt.comp) {

        ssetCmp <- colnames(sCrit)
        for (clu in unique(saffectElt)) {

          index <- which(saffectElt == clu)

          v <- sCrit[fres$nbOpt, index]
          if (length(v) > 1)
            for (i in seq_len(length(v) - 1))
              if (v[i] == v[i + 1]) {
                ssetCmp[index][i] <- paste0(ssetCmp[index][i], " = ")
              } else {
                ssetCmp[index][i] <- paste0(ssetCmp[index][i], " > ")
              }
        }

        names(saffectElt) <- ssetCmp
        plot_clusters_content(saffectElt, opt.sort = FALSE, sepChar = "")
      }
    }
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Test significant assemblages                                            ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Evaluate the weight of each assemblage on functional clustering
#'
#' @description Evaluate by cross-validation (leave-one-out)
#'  the effect induced by leaving out each assemblage
#'  on the result of a functional clustering.
#'
#' @usage
#' ftest_assemblages(fres, opt.nbMax = fres$nbOpt,
#'                   opt.R2 = FALSE, opt.plot = FALSE)
#'
#' @inheritParams ftest
#'
#' @details Each assemblage of the dataset is successively removed,
#'  the remaining assemblage collection is functionally clustered,
#'  then indices of distance between clustering trees
#'  with and without the assemblage are computed.
#'  The assemblages can then be hierarchised depending on
#'  the distance induced by their removing from dataset.
#'  The used distance criteria are :
#'   "Czekanowski_Dice", "Folkes_Mallows", "Jaccard", "Kulczynski",
#'   "Precision", "Rand", "Recall", "Rogers_Tanimoto", "Russel_Rao",
#'   "Sokal_Sneath1" and "Sokal_Sneath2" index.
#'   For more informations, see the notice of R-package \code{clusterCrit}.
#'   The test is time-consuming.
#'
#' @references
#' Package "clusterCrit": Clustering Indices,
#'   by Bernard Desgraupes (University of Paris Ouest - Lab Modal'X)
#'
#' @importFrom clusterCrit getCriteriaNames extCriteria
#'
#' @importFrom graphics text
#'
#' @return a list containing a matrix by clustering index.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

ftest_assemblages <- function(fres, opt.nbMax = fres$nbOpt,
                              opt.R2 = FALSE, opt.plot = FALSE) {


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  nbAss  <- fres$nbAss
  setXpr <- unique(names(fres$xpr))
  setAss <- which(names(fres$xpr) == setXpr[1])
  names(setAss) <- rownames(fres$mOccur[setAss, ])

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  lAffectElt <- vector(mode = "list", length = fres$nbElt)
  for (nb in seq_len(fres$nbElt))
    lAffectElt[[nb]] <- cut_ftree(fres$tree.I, nb)


  # plot the reference tree if opt.plot = TRUE

  if ((getOption("verbose") == TRUE) && (opt.plot == TRUE)) {
    cols <- fcolours()[cut_ftree(fres$tree.II, fres$nbOpt)]
    plot_ftree(tree = fres$tree.II, cols = cols)
    graphics::text(x = 1, y = 0, labels = "Reference tree", pos = 4)
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  indCrit   <- c(1:2,4:5,8:14)
  critNames <- clusterCrit::getCriteriaNames(FALSE)[indCrit]    #external crits
  vCriteria <- numeric(length(critNames))
  names(vCriteria) <- c(critNames)

  mCrit <- array(0, dim = c(fres$nbElt, fres$nbAss, length(critNames)))
  dimnames(mCrit) <- list(seq_len(fres$nbElt),
                          unique(rownames(fres$mOccur)), critNames)

  if (opt.R2 == TRUE) {
    mStats <- array(0, dim = c(fres$nbElt, fres$nbAss, 2))
    dimnames(mStats) <- list(seq_len(fres$nbElt),
                             unique(rownames(fres$mOccur)), c("R2", "E"))
  }

  for (ass in seq_along(setAss)) {

    # compute the primary tree

    assemblage <- setAss[ass]

    main <- paste0("The assemblage '", names(assemblage), "' is removed")
    if (getOption("verbose") == TRUE) cat(main, "\n")

    indAss <- setdiff(setAss, ass)
    index  <- which(names(fres$fobs) %in% names(setAss)[indAss])

    tree.I <- fit_ftree(fres$fobs[index], fres$mOccur[index, ],
                        fres$xpr[index], fres$affectElt,
                        fres$opt.method, fres$opt.mean, fres$opt.model,
                        opt.nbMax  )

    for (nb in seq_len(fres$nbElt)) {
      vCriteria[] <- as.vector(unlist(clusterCrit::extCriteria(
        as.integer(cut_ftree(tree.I, nb) ),
        as.integer(lAffectElt[[nb]]),
        crit = "all")))[indCrit]
      mCrit[nb, assemblage, ] <- vCriteria
    }

    if (sum(is.na(mCrit[fres$nbElt, assemblage, ])) == fres$nbAss)
      mCrit[fres$nbElt, assemblage, ] <- rep(0, fres$nbAss)

    # compute the secondary tree

    if (opt.R2 == TRUE) {
      res <- validate_ftree(tree.I,
                            fres$fobs[index], fres$mOccur[index, ],
                            fres$xpr[index],
                            fres$opt.method,  fres$opt.mean, fres$opt.model,
                            fres$opt.jack, fres$jack,
                            opt.nbMax  )

      mStats[ , assemblage, "R2"] <- res$tStats[ , "R2cal"]
      mStats[ , assemblage, "E"]  <- res$tStats[ , "R2prd"]
    }

    # plot the resulting tree if opt.plot = TRUE

    if ((getOption("verbose") == TRUE) && (opt.plot == TRUE)) {
      if (opt.R2 == TRUE) {
        plot_ftree(res$tree.II, cols = cols)
      } else {
        plot_ftree(tree.I, cols = cols)
      }
      graphics::text(x = 1, y = 0, labels = main, pos = 4)
    }
  }

  coord <- which(mCrit == "NaN", arr.ind = TRUE)
  mCrit[coord] <- ifelse(amean(mCrit[fres$nbElt, , ]) < 0.1, 0, 1)


  # Organize the results

  rtest <- vector(mode = "list", length = length(critNames))
  names(rtest) <- critNames
  for (ind in seq_along(critNames))
    rtest[[ind]] <- as.matrix(mCrit[ , , critNames[ind]])

  if (opt.R2 == TRUE) {
    mStats[fres$nbElt, , "R2"] <- mStats[fres$nbElt - 1, , "R2"]
    mStats[fres$nbElt, , "E"]  <- mStats[fres$nbElt - 1, , "E"]
    rtest$R2 <- as.matrix(mStats[ , , "R2"])
    rtest$E  <- as.matrix(mStats[ , , "E"])
  }


  return(rtest)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Plot the evaluation of weight of each assemblage
#'   on functional clustering
#'
#' @description Evaluate by cross-validation (leave-one-out)
#'  the effect induced by leaving out each assemblage
#'  on the result of a functional clustering.
#'
#' @usage
#' ftest_plot_assemblages(fres, rtest,
#'                        main = "Title", opt.crit = "Jaccard", opt.ass = NULL)
#'
#' @inheritParams ftest_plot
#'
#' @details The trees obtained by leaving out each assemblage
#' are compared to the reference tree obtained with all assemblages
#' using different criteria :
#' "Czekanowski_Dice", "Folkes_Mallows", "Jaccard", "Kulczynski",
#'  "Precision", "Rand", "Recall", "Rogers_Tanimoto", "Russel_Rao",
#'  "Sokal_Sneath1" and "Sokal_Sneath2" index.
#'  For more informations, see the notice of R-package \code{clusterCrit}.
#'
#' @references
#' Package "clusterCrit": Clustering Indices,
#'   by Bernard Desgraupes (University of Paris Ouest - Lab Modal'X)
#'
#' @importFrom graphics plot lines axis
#'
#' @return Nothing. It is a procedure.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

ftest_plot_assemblages <- function(fres, rtest,
                                   main = "Title", opt.crit = "Jaccard",
                                   opt.ass = NULL) {


  affectElt <- cut_ftree(fres$tree.II, fres$nbOpt)
  indClu    <- seq_len(fres$nbOpt)

  setCrit <- names(rtest)
  indCrit <- which(setCrit %in% opt.crit)

  assMot  <- fres$mMotifs[fres$nbOpt, ][seq_len(fres$nbAss)]
  setMot  <- unique(assMot)

  assNam  <- name_motifs(cut_ftree(
    fres$tree.II, fres$nbOpt), fres$mOccur)[seq_len(fres$nbAss)]
  nameMot <- unique(assNam)


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Check options
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  all.ass <- list("all.together", "motifs.together", "assemblages.together",
                  "motifs.bymot", "assemblages.bymot", "sorted.leg",
                  "all")
  names(all.ass) <- all.ass

  opt.ass <- check_plot_options(opt.ass, all.ass, main = "opt.assemblages",
                                affectElt = affectElt)


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  cols    <- extend_vector(fcolours(), fres$nbAss)
  pchs    <- extend_vector(fsymbols(), fres$nbAss)

  crit <- 5
  for (crit in indCrit) {

    mCrit    <- as.matrix(rtest[[crit]])
    mainCrit <- paste(main, "ass", setCrit[crit], sep = " / ")

    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    # plot general median
    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    if ("all.together" %in% opt.ass) {

      mainLoc <- paste(mainCrit, "general median index", sep = " / ")

      graphics::plot(x = mCrit[indClu, 1], xlim = c(0, 1),
                     y = fres$tree.II$cor[indClu], ylim = c(0, 1),
                     xlab = paste0(setCrit[crit], " index"), ylab = "R2",
                     main = mainLoc,
                     type = "n", tck = 0.02, las = 2)

      for (clu in seq_len(fres$nbOpt - 1))
        graphics::lines(x = c(median(mCrit[clu, ]), median(mCrit[clu + 1, ])),
                        y = c(fres$tree.II$cor[clu],
                              fres$tree.II$cor[clu + 1]),
                        lty = "dashed", col = fcolours()[1])

      for (clu in seq_len(fres$nbOpt))
        points_sd(x = mCrit[clu, ],  y = fres$tree.II$cor[clu],
                  opt = c("median", "mean"),
                  fig = fsymbols()[1], col = fcolours()[1], bg = fcolours()[1])
    }



    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    # plot median component clusters, all together
    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    if ("motifs.together" %in% opt.ass) {

      mainLoc <- paste(mainCrit, "median index by assembly motifs", sep = " / ")

      graphics::plot(x = mCrit[indClu, 1], xlim = c(0, 1),
                     y = fres$tree.II$cor[indClu], ylim = c(0, 1),
                     xlab = paste0(setCrit[crit], " index"), ylab = "R2",
                     main = mainLoc,
                     type = "n", tck = 0.02, las = 2)

      for (mot in seq_along(setMot)) {

        setAss <- which(fres$mMotifs[fres$nbOpt,
                                     seq_len(fres$nbAss)] == setMot[mot])

        vtmp   <- apply(mCrit[indClu, setAss, drop = FALSE],
                        MARGIN = 1, FUN = median)

        graphics::points(x = vtmp,
                         y = fres$tree.II$cor[indClu],
                         type = "b", lty = "solid",
                         pch = pchs[setMot[mot]], col = cols[setMot[mot]],
                         cex = 2)

        graphics::axis(labels = nameMot[mot],
                       tick = FALSE, las = 2, side = 3,
                       at  = vtmp[fres$nbOpt],
                       pos = fres$tree.II$cor[fres$nbOpt] + 0.02)
      }
    }


    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    # plot median assemblages, assembly motif by assembly motif
    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    if ("motifs.bymot" %in% opt.ass) {

      for (mot in seq_along(setMot)) {

        mainLoc <- paste(mainCrit,
                         paste0("median index for assembly motif '",
                                nameMot[mot], "'"),
                         sep = " / ")

        setAss <- which(fres$mMotifs[fres$nbOpt,
                                     seq_len(fres$nbAss)] == setMot[mot])

        graphics::plot(x = mCrit[indClu, 1, drop = FALSE], xlim = c(0, 1),
                       y = fres$tree.II$cor[indClu], ylim = c(0, 1),
                       xlab = paste0(setCrit[crit], " index"), ylab = "R2",
                       main = mainLoc,
                       type = "n", tck = 0.02, las = 2)

        for (lev in seq_len(fres$nbOpt - 1))
          graphics::lines(x = c(median(mCrit[lev, setAss]),
                                median(mCrit[lev + 1, setAss])),
                          y = c(fres$tree.II$cor[lev],
                                fres$tree.II$cor[lev + 1]),
                          lty = "dashed", col = cols[setMot[mot]])

        for (lev in seq_len(fres$nbOpt))
          points_sd(x = mCrit[lev, setAss], y = fres$tree.II$cor[lev],
                    opt = c("median", "mean"),
                    fig = pchs[setMot[mot]], col = cols[setMot[mot]],
                    cex = 2)

        graphics::axis(labels = nameMot[mot],
                       tick = FALSE, las = 2, side = 3,
                       at  = median(mCrit[fres$nbOpt, setAss]),
                       pos = fres$tree.II$cor[fres$nbOpt] + 0.02)
      }
    }


    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    # plot a graph by assemblage, assembly motif by assembly motif
    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    if ("assemblages.bymot" %in% opt.ass) {

      for (mot in seq_along(setMot)) {

        mainLoc  <- paste(mainCrit,
                          paste0("assemblage index within assembly motif '",
                                 nameMot[mot], "'"),
                          sep = " / ")

        setAss <- which(fres$mMotifs[fres$nbOpt,
                                     seq_len(fres$nbAss)] == setMot[mot])

        graphics::plot(x = mCrit[indClu, 1, drop = FALSE], xlim = c(0, 1),
                       y = fres$tree.II$cor[indClu], ylim = c(0, 1),
                       xlab = paste0(setCrit[crit], " index"), ylab = "R2",
                       main = mainLoc,
                       type = "n", tck = 0.02, las = 2)

        for (ass in seq_along(setAss)) {

          graphics::points(x = mCrit[indClu, setAss[ass]],
                           y = fres$tree.II$cor[indClu],
                           type = "b", lty = "solid",
                           col = cols[setMot[mot]],
                           pch = extend_vector(fsymbols(),
                                               length(setAss))[ass],
                           cex = 2)

          graphics::axis(labels = unique(rownames(fres$mOccur))[setAss[ass]],
                         tick = FALSE, las = 2, side = 3,
                         at  = mCrit[fres$nbOpt, setAss[ass]],
                         pos = fres$tree.II$cor[fres$nbOpt] + 0.02)
        }
      }
    }


    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    # plot a graph by assemblage
    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    if ("assemblages.together" %in% opt.ass) {

      mainLoc <- paste(mainCrit, "all assemblages", sep = " / ")

      graphics::plot(x = mCrit[indClu, 1], xlim = c(0, 1),
                     y = fres$tree.II$cor[indClu], ylim = c(0, 1),
                     xlab = paste0(setCrit[crit], " index"), ylab = "R2",
                     main = mainLoc,
                     type = "n", tck = 0.02, las = 2)

      for (mot in seq_along(setMot)) {

        setAss <- which(fres$mMotifs[fres$nbOpt,
                                     seq_len(fres$nbAss)] == setMot[mot])

        for (ass in seq_along(setAss)) {

          graphics::points(x = mCrit[indClu, setAss[ass]],
                           y = fres$tree.II$cor[indClu],
                           type = "b", lty = "solid",
                           col = cols[setMot[mot]],
                           pch = extend_vector(fsymbols(),
                                               length(setAss))[ass],
                           cex = 2)

          graphics::axis(labels = unique(rownames(fres$mOccur))[setAss[ass]],
                         tick = FALSE, las = 2, side = 3,
                         at  = mCrit[fres$nbOpt, setAss[ass]],
                         pos = fres$tree.II$cor[fres$nbOpt] + 0.02)
        }
      }
    }


    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    # plot a tree with components sorted by importance and/or the legend
    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    if ("sorted.leg" %in% opt.ass)  {

      sindex <- seq_len(fres$nbAss)
      for (mot in seq_along(setMot)) {

        index <- which(assMot == setMot[mot])

        lCrit <- list()
        for (i in 1:fres$nbOpt) lCrit[[i]] <- mCrit[i, index]
        ssindex <- do.call(what = order, args = lCrit)

        sindex[index] <- sindex[index][ssindex]
      }

      ssetCmp <- names(mCrit[fres$nbOpt, ])[sindex]
      for (mot in seq_along(setMot)) {

        index <- which(fres$mMotifs[fres$nbOpt,
                                    seq_len(fres$nbAss)] == setMot[mot])
        v <- mCrit[fres$nbOpt, ][sindex][index]
        if (length(v) > 1)
          for (i in seq_len(length(v) - 1))
            if (v[i] == v[i + 1]) {
              ssetCmp[index][i] <- paste0(ssetCmp[index][i], " = ")
            } else {
              ssetCmp[index][i] <- paste0(ssetCmp[index][i], " > ")
            }
      }

      tmp <- fres$fobs
      names(tmp)[sindex] <- ssetCmp

      plot_motifs_content(fobs     = tmp[sindex],
                          assMotif = assMot[sindex],
                          assNames = assNam[sindex],
                          opt.sort = FALSE,
                          sepChar  = "")

    }
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Test significant performances                                           ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Evaluate the weight of each performance on functional clustering
#'
#' @description Evaluate by cross-validation (leave-one-out)
#'  the effect induced by leaving out each performance
#'  on the result of a functional clustering.
#'
#' @usage
#' ftest_performances(fres, opt.nbMax = fres$nbOpt,
#'                    opt.R2 = FALSE, opt.plot = FALSE)
#'
#' @inheritParams ftest
#'
#' @details Each performance of the dataset is successively removed,
#'  the remaining performance collection is functionally analysed,
#'  then indices of distance between clustering trees
#'  with and without the performance are computed.
#'  The performances can then be hierarchised depending on
#'  the distance induced by their removing from dataset.
#'  The used distance criteria are :
#'   "Czekanowski_Dice", "Folkes_Mallows", "Jaccard", "Kulczynski",
#'   "Precision", "Rand", "Recall", "Rogers_Tanimoto", "Russel_Rao",
#'   "Sokal_Sneath1" and "Sokal_Sneath2" index.
#'   For more informations, see the notice of R-package \code{clusterCrit}.
#'   The test is time-consuming.
#'
#' @references
#' Package "clusterCrit": Clustering Indices,
#'   by Bernard Desgraupes (University of Paris Ouest - Lab Modal'X)
#'
#' @importFrom clusterCrit getCriteriaNames extCriteria
#'
#' @importFrom graphics text
#'
#' @return a list containing a matrix by clustering index.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

ftest_performances <- function(fres, opt.nbMax = fres$nbOpt,
                               opt.R2 = FALSE, opt.plot = FALSE) {


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  setXpr <- unique(names(fres$xpr))

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  lAffectElt <- vector(mode = "list", length = fres$nbElt)
  for (nb in seq_len(fres$nbElt))
    lAffectElt[[nb]] <- cut_ftree(fres$tree.I, nb)

  # plot the reference tree if opt.plot = TRUE

  if ((getOption("verbose") == TRUE) && (opt.plot == TRUE)) {
    cols <- fcolours()[cut_ftree(fres$tree.II, fres$nbOpt)]
    plot_ftree(tree = fres$tree.II, cols = cols)
    graphics::text(x = 1, y = 0, labels = "Reference tree", pos = 4)
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  indCrit   <- c(1:2,4:5,8:14)
  critNames <- clusterCrit::getCriteriaNames(FALSE)[indCrit]    #external crits
  vCriteria <- numeric(length(critNames))
  names(vCriteria) <- c(critNames)

  mCrit <- array(0, dim = c(fres$nbElt, fres$nbXpr, length(critNames)))
  dimnames(mCrit) <- list(seq_len(fres$nbElt), setXpr, critNames)

  if (opt.R2 == TRUE) {
    mStats <- array(0, dim = c(fres$nbElt, fres$nbXpr, 2))
    dimnames(mStats) <- list(seq_len(fres$nbElt), setXpr, c("R2", "E"))
  }

  ipr <- 1
  for (ipr in seq_along(setXpr)) {

    performance <- setXpr[ipr]

    main <- paste0("The performance '", setXpr[ipr], "' is removed")
    if (getOption("verbose") == TRUE) cat(main, "\n")

    indXpr <- setdiff(setXpr, performance)
    index  <- which(names(fres$xpr) %in% indXpr)

    tree.I <- fit_ftree(fres$fobs[index], fres$mOccur[index, ],
                        fres$xpr[index],  fres$affectElt,
                        fres$opt.method,  fres$opt.mean,
                        fres$opt.model,
                        opt.nbMax  )

    for (nb in seq_len(fres$nbElt)) {
      vCriteria[] <- as.vector(unlist(clusterCrit::extCriteria(
        as.integer(cut_ftree(tree.I, nb) ),
        as.integer(lAffectElt[[nb]]),
        crit = "all")))[indCrit]
      mCrit[nb, performance, ] <- vCriteria
    }

    # compute the secondary tree

    if (opt.R2 == TRUE) {
      res <- validate_ftree(tree.I,
                            fres$fobs[index], fres$mOccur[index, ],
                            fres$xpr[index],
                            fres$opt.method,  fres$opt.mean, fres$opt.model,
                            fres$opt.jack, fres$jack,
                            opt.nbMax  )

      mStats[seq_len(fres$nbElt), performance, "R2"] <- res$tStats[ , "R2cal"]
      mStats[seq_len(fres$nbElt), performance, "E"]  <- res$tStats[ , "R2prd"]
    }

    # plot the resulting tree if opt.plot = TRUE

    if ((getOption("verbose") == TRUE) && (opt.plot == TRUE)) {
      if (opt.R2 == TRUE) {
        plot_ftree(res$tree.II, cols = cols)
      } else {
        plot_ftree(tree.I, cols = cols)
      }
      graphics::text(x = 1, y = 0, labels = main, pos = 4)
    }
  }

  coord <- which(mCrit == "NaN", arr.ind = TRUE)
  mCrit[coord] <- ifelse(amean(mCrit[fres$nbElt, , ]) < 0.1, 0, 1)


  # Organize the results

  rtest <- vector(mode = "list", length = length(critNames))
  names(rtest) <- critNames
  for (ind in seq_along(critNames))
    rtest[[ind]] <- as.matrix(mCrit[ , , critNames[ind]])

  if (opt.R2 == TRUE) {
    mStats[fres$nbElt, , "R2"] <- mStats[fres$nbElt - 1, , "R2"]
    mStats[fres$nbElt, , "E"]  <- mStats[fres$nbElt - 1, , "E"]
    rtest$R2 <- as.matrix(mStats[ , , "R2"])
    rtest$E  <- as.matrix(mStats[ , , "E"])
  }


  return(rtest)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Plot the evaluation of weight of each performance
#'   on functional clustering
#'
#' @description Evaluate by cross-validation (leave-one-out)
#'  the effect induced by leaving out each performance
#'  on the result of a functional clustering.
#'
#' @usage
#' ftest_plot_performances(fres, rtest,
#'                         main = "Title", opt.crit = "Jaccard", opt.perf = NULL)
#'
#' @inheritParams ftest_plot
#'
#' @details The trees obtained by leaving out each performance
#' are compared to the reference tree obtained with all performances
#' using different criteria :
#' "Czekanowski_Dice", "Folkes_Mallows", "Jaccard", "Kulczynski",
#'  "Precision", "Rand", "Recall", "Rogers_Tanimoto", "Russel_Rao",
#'  "Sokal_Sneath1" and "Sokal_Sneath2" index.
#'  For more informations, see the notice of R-package \code{clusterCrit}.
#'
#' @references
#' Package "clusterCrit": Clustering Indices,
#'   by Bernard Desgraupes (University of Paris Ouest - Lab Modal'X)
#'
#' @importFrom graphics plot lines axis
#'
#' @return Nothing. It is a procedure.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

ftest_plot_performances <- function(fres, rtest,
                                    main = "Title", opt.crit = "Jaccard",
                                    opt.perf = NULL) {

  setCrit <- names(rtest)[names(rtest) %in% opt.crit]
  indCrit <- which(setCrit %in% opt.crit)

  indClu  <- seq_len(fres$nbOpt)

  cols    <- extend_vector(fcolours(), fres$nbXpr)
  pchs    <- extend_vector(fsymbols(), fres$nbXpr)

  setXpr  <- unique(names(fres$xpr))

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Check options
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  all.perf <- list("all.together", "performances.together", "sorted.leg",
                   "all")
  names(all.perf) <- all.perf

  opt.perf <- check_plot_options(opt.perf, all.perf, main = "opt.performances",
                                 affectElt = fres$affectElt)


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  crit <- 1
  for (crit in indCrit) {

    mCrit    <- as.matrix(rtest[[setCrit[crit]]])
    mainCrit <- paste(main, "perf", setCrit[crit], sep = " / ")

    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    # plot general median
    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    if ("all.together" %in% opt.perf) {

      mainLoc <- paste(mainCrit, "general median", sep = " / ")

      graphics::plot(x = mCrit[indClu, 1], xlim = c(0, 1),
                     y = fres$tree.II$cor[indClu], ylim = c(0, 1),
                     xlab = paste0(setCrit[crit], " index"), ylab = "R2",
                     main = mainLoc,
                     type = "n", tck = 0.02, las = 2)

      for (clu in seq_len(fres$nbOpt - 1))
        graphics::lines(x = c(median(mCrit[clu, ]), median(mCrit[clu + 1, ])),
                        y = c(fres$tree.II$cor[clu],
                              fres$tree.II$cor[clu + 1]),
                        lty = "dashed", col = fcolours()[1])

      for (clu in seq_len(fres$nbOpt))
        points_sd(x = mCrit[clu, ], y = fres$tree.II$cor[clu],
                  opt = c("median", "mean"),
                  fig = fsymbols()[1], col = fcolours()[1], bg = fcolours()[1])
    }


    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    # plot a graph by performance
    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    if ("performances.together" %in% opt.perf) {

      mainLoc <- paste(mainCrit, "all performances", sep = " / ")

      graphics::plot(x = mCrit[indClu, 1], xlim = c(0, 1),
                     y = fres$tree.II$cor[indClu], ylim = c(0, 1),
                     xlab = paste0(setCrit[crit], " index"), ylab = "R2",
                     main = mainLoc,
                     type = "n", tck = 0.02, las = 2)

      for (ipr in seq_along(setXpr)) {

        graphics::points(x = mCrit[indClu, setXpr[ipr]],
                         y = fres$tree.II$cor[indClu],
                         type = "b", lty = "solid",
                         col = cols[ipr],
                         pch = pchs[ipr],
                         cex = 2)

        graphics::axis(labels = setXpr[ipr],
                       tick = FALSE, las = 2, side = 3,
                       at  = mCrit[fres$nbOpt, setXpr[ipr]],
                       pos = fres$tree.II$cor[fres$nbOpt] + 0.02)
      }
    }


    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    # plot a tree with components sorted by importance and/or the legend
    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    if ("sorted.leg" %in% opt.perf)  {

      lCrit <- list()
      for (i in 1:fres$nbOpt) lCrit[[i]] <- mCrit[i, ]
      sindex  <- do.call(what = order, args = lCrit)
      ssetXpr <- setXpr[sindex]

      v <- mCrit[fres$nbOpt, ][sindex]
      for (i in seq_len(fres$nbXpr - 1))
        if (v[i] == v[i + 1]) {
          ssetXpr[i] <- paste0(ssetXpr[i], " = ")
        } else {
          ssetXpr[i] <- paste0(ssetXpr[i], " > ")
        }

      label <- paste("Decreasing order of performance effects for",
                     setCrit[crit], "index:", sep = " ")
      tab <- c(label, "", concat_by_line(ssetXpr, sepChar = ""))

      plot_by_page(tab)
    }
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# End of file                                                              ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
