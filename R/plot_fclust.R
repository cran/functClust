#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                                 PLOT_FCLUST.R
#
#                set of functions for plotting the results
#                of components or motifs sorting and clustering
#
#                         Benoit JAILLARD, summer 2020
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#' @include
#' stats.R
#' tools.R
#' labelling.R
#'
NULL



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Graphic tools                                                           ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Plot Option.tree                                                        ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Simplify a tree by keeping only significant components
#'
#' @description Take a tree,
#' then keep only the components identified in \code{signifElt}
#' by removing all other components.
#'
#' @usage
#' simplify_ftree(tree, signifElt)
#'
#' @inheritParams plot_ftree
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
#' The vector of coefficients of determination
#' determines the level of each split of a cluster into two new clusters.
#'
#' @return Return a tree simplified,
#'  that is without less of components.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

simplify_ftree <- function(tree, signifElt) {

  tree$aff <- tree$aff[ , signifElt]

  fct   <- function(x) { length(unique(x)) }
  laff  <- apply(tree$aff, MARGIN = 1, FUN = fct)
  ind1  <- which(table(laff) > 1)

  index <- NULL
  for (elt in seq_along(ind1)) index <- c(index, which(laff == ind1[elt])[-1])

  if (length(index) > 0) {

    tree$cor <- tree$cor[-index]
    tree$aff <- tree$aff[-index, ]

    for (elt in 2:dim(tree$aff)[2]) {

      l1 <- tree$aff[elt - 1, ]
      l2 <- tree$aff[elt, ]

      for (clu in 1:elt) {

        if (length(which(l2 == clu)) == 0) {

          index <- which(l1 == clu)
          for (i in seq_along(index)) {
            elt2  <- elt
            while ((elt2 <= dim(tree$aff)[2])
                   && (tree$aff[elt2, index[i]] == l2[index[i]]) ) {
              tree$aff[elt2, index[i]] <- clu
              elt2 <- elt2 + 1
            }
          }
        }
        tree$aff[elt, ] <- compact_index(tree$aff[elt, ])
      }
    }

    rownames(tree$aff) <- seq_len(sum(signifElt != 0))
  }

  return(tree)
}



##xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Write the components belonging to each functional group
#'
#' @description Write the components belonging to each functional group,
#' on several pages if need.
#'
#' @usage
#' plot_clusters_content(affectElt,
#'                       opt.sort = TRUE,
#'                       sepChar  = ", ",
#'                       nbline   = 25)
#'
#' @param affectElt a vector of integers.
#' The vector contains the labels of the different functional groups
#' to which each component belongs.
#' Each functional group is labelled as an integer.
#'
#' @param opt.sort a logical.
#' If \code{opt.sort = TRUE}, the names of assemblages
#' are sorted alphabetically inside each motif.
#'
#' @param sepChar a string.
#' The string is a separator between each element of list of strings.
#'
#' @param nbline an integer,
#' that corresponds to a number of lines to plot by page.
#'
#' @details None.
#'
#' @return Nothing. It is a procedure.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot_clusters_content <- function(affectElt,
                                  opt.sort = TRUE,
                                  sepChar  = ", ",
                                  nbline   = 25 ) {

  noms   <- name_clusters(affectElt)

  setClu <- table(noms)
  nbrow  <- length(setClu)

  # compute the number of necessary columns
  tmp    <- list()

  for (i in seq_len(nbrow))
    tmp[[i]] <- names(affectElt)[which(noms == names(setClu)[i])]

  if (opt.sort == TRUE)
    for (i in seq_len(nbrow)) tmp[[i]] <- sort(tmp[[i]])

  # adjust the length of rownames
  lmax <- max(nchar(names(setClu)))
  for (i in seq_len(nbrow)) {
    a <- names(setClu)[i]
    while (nchar(a) < lmax) a <- paste(" ", a, sep = "")
    names(tmp[[i]])[1] <- a
  }

  blanc <- NULL
  for (k in 1:(lmax + 5)) blanc <- paste(blanc, " ", sep = "")

  # concat the matrix by row
  tab <- c("Cluster content:", "")
  for (i in seq_len(nbrow)) {

    tpp    <- concat_by_line(tmp[[i]], sepChar = sepChar)
    tpp[1] <- paste0(names(tmp[[i]])[1], " = { ", tpp[1])
    if (length(tpp) > 1)
      for (j in 2:length(tpp)) tpp[j] <- paste0(blanc, tpp[j])
    tpp[length(tpp)] <- paste0(tpp[length(tpp)], " }")

    tab  <- c(tab, tpp)
  }

  # plot the legend
  plot_by_page(tab, nbline)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Cut a tree at a given level
#'
#' @description Take a tree and cut the matrix of component affectation
#' at a given level.
#'
#' @usage cut_ftree(tree, nbcl)
#'
#' @inheritParams plot_ftree
#'
#' @param nbcl an integer. \code{1 <= nbcl <= dim(tree$aff)[1]}.
#'
#' @details A hierarchical tree is recorded
#' as a list of matrix of component affectation
#' and a vector of coefficients of determination.\cr
#'
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
#' The function \code{cut_ftree} samples the line \code{nbcl}
#' of the square matrix of component affectation.
#'
#' @return Return a vector of component affectation.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

cut_ftree <- function(tree, nbcl) {

  return( compact_index(tree$aff[nbcl, ]) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Sort the resulting file of tree
#'
#' @description Take a tree, group the components by cluster,
#' and sort the clusters by increasing label.
#' Increasing label corresponds to increasing efficiency of the clusters.
#'
#' @usage sort_ftree(tree, index.return = FALSE)
#'
#' @inheritParams plot_ftree
#'
#' @param index.return a logical.
#' If \code{index.return = TRUE},
#'  \code{sort_ftree} only returns the sorted matrix of component affectation.
#' If \code{index.return = FALSE},
#'  \code{sort_ftree} returns the sorted matrix of component affectation,
#'  the sorted indices and the sorted names of components.
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
#' The last line separate each component in a singleton cluster.
#'
#' @return
#' \item{index.return = TRUE}{\code{sort_ftree} returns
#' the sorted matrix of component affectation.
#' }
#' \item{index.return = FALSE}{\code{sort_ftree} returns
#'  the sorted matrix of component affectation,
#'  the sorted indices and the sorted names of components.
#' }
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

sort_ftree <- function(tree, index.return = FALSE) {

  X      <- shift_affectElt(tree$aff)

  nbline <- dim(X)[1]
  nbitem <- dim(X)[2]


  # Separate first the sorted from unsorted components
  mask <- ordre <- seq_len(nbitem)
  for (lin in 2:nbline) {
    index       <- sort(x = X[lin, mask], index.return = TRUE)
    index       <- index$ix + nbitem - length(mask)
    X[ , mask]  <- X[ , index]
    ordre[mask] <- ordre[index]
    mask        <- mask[X[lin, mask] == max(X)]
  }

  # Sort the leaves of tree
  new <- old <- seq_len(nbitem)
  for (lin in seq_len(nbline)) {
    for (item in seq_len(nbitem))
      old[item] <- which(unique(X[lin,]) == X[lin, item])
    new             <- sort(old, index.return = TRUE)
    X[ , seq_len(nbitem)]  <- X[ , new$ix]
    ordre[seq_len(nbitem)] <- ordre[new$ix]
  }

  # Re-numeration of leaves
  old <- integer(nbline)
  for (lin in seq_len(nbline)) {
    v  <- setdiff(unique(X[lin, ]), unique(old))
    if (length(v) == 1) old[lin] <- v
  }

  newX <- matrix(0, nrow = nbline, ncol = nbitem)
  for (lin in seq_len(nbline))
    for (j in seq_len(lin)) newX[lin, which(old[j] == X[lin, ])] <- j

  colnames(newX) <- colnames(X)[ordre]

  if (index.return == TRUE) {
    res        <- list(newX, ordre, colnames(newX))
    names(res) <- c("x", "ix", "noms")
  } else {
    res <- newX
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Plot a hierarchical tree
#'
#' @description Take a tree and plot its representation,
#'  trunk at bottom and leaves at top.
#'
#' @usage
#' plot_ftree(tree, signifElt = rep(TRUE, length(tree$cor)), cols = "black",
#'           main = "", opt.sort = TRUE)
#'
#' @param tree a list of a matrix of component affectation
#' and a vector of coefficients of determination.
#'
#' @param signifElt a vector of logical.
#' \code{TRUE} means
#' that deleting the component changes the component clustering.
#' \code{FALSE} means
#' that deleting the component does not change the clustering.
#'
#' @param cols a vector of colours, sorted as the components of tree.
#'
#' @param main a string, that appears as main title of the plot.
#'
#' @param opt.sort a logical. \code{opt.sort} = TRUE by default.
#' If \code{opt.sort} = FALSE, the tree is not re-organized.
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
#' The vector of coefficients of determination
#' determines the level of each split of a cluster into two new clusters.
#'
#' @return Nothing. It is a procedure that plots a hierarchical tree.
#'
#' @importFrom graphics plot points lines axis title
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot_ftree <- function(tree,
                       signifElt = rep(TRUE, length(tree$cor)),
                       cols      = "black", main = "",
                       opt.sort  = TRUE) {

  vpch            <- numeric(dim(tree$aff)[2])
  vpch[ ]         <- 1
  vpch[signifElt] <- 19

  if (length(tree$cor) > 1) {
    if (opt.sort == TRUE) {
      tmpx <- sort_ftree(tree, index.return = TRUE)
      X    <- tmpx$x
      vpch <- vpch[tmpx$ix]
      if (length(cols) > 1) cols <- cols[tmpx$ix]
    } else {
      X    <- tree$aff
    }
  } else {
    X    <- shift_affectElt(tree$aff)
    vpch <- vpch[1]
    if (length(cols) > 1) cols <- cols[1]
  }

  Y   <- tree$cor
  Y[Y < 0] <- 0


  nbline <- dim(X)[1]
  nbitem <- dim(X)[2] + 1
  YY     <- c(0,Y)

  xx     <- matrix(0, nrow = nbline, ncol = nbline)

  # plot the framework
  graphics::plot(x = seq_len(nbitem), y = seq(1/nbitem, 1, 1/nbitem),
                 xlab = "component",  ylab = "R2-value",
                 xlim = c(1, nbitem), ylim = c(0, 1),
                 type = "n", tck = 0.02, las = 1)

  # plot the vertical lines
  for (lin in 1:nbline)
    for (leave in 1:lin) {
      xx[lin, leave] <- 1 + mean(which(X[lin, ] == leave))
      graphics::lines(x = rep(xx[lin, leave], 2), y = c(YY[lin], YY[lin + 1]),
                      lty = "solid")
    }

  # plot the horizontal jonction lines
  if (nbline > 1) for (lin in 2:nbline) {
    XX <- which(xx[lin, ] != xx[lin - 1, ])
    if (length(XX) == 1) XX <- rep(XX, 2)
    graphics::lines(x = xx[lin, XX], y = rep(YY[lin], 2), lty = "solid" )
  }

  # plot the final horizontal jonction lines
  setLeave <- sort(unique(X[nbline, ]))
  for (leave in seq_along(setLeave)) {
    tmp <- which(X[nbline, ] == setLeave[leave])
    graphics::lines(x = c(min(tmp), max(tmp)) + 1,
                    y = rep(YY[nbline + 1], 2),
                    lty = "solid" )
  }

  # plot the symbols of the most likely partition
  xpos <- 1 + c(1:(nbitem - 1))
  ypos <- YY[length(YY)] + 0.02

  if (length(cols) > 1) {
    graphics::points(x = xpos, y = rep(ypos, (nbitem - 1)), pch = vpch,
                     cex = 2, col = cols, bg = "white")
  } else {
    graphics::points(x = xpos, y = rep(ypos, (nbitem - 1)), pch = vpch,
                     bg = "white")
  }

  #  plot the names of components
  graphics::axis(side = 3, at = 2:nbitem, labels = colnames(X),
                 tick = FALSE, las = 2, pos = ypos)

  graphics::title(main)
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Plot Option.performances                                                 ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Add assemblage names on a plot
#'
#' @description Plot the name of each assemblage near its performance
#' on a graph plotted using the functions
#' \code{\link{plot_prediction_simple}} or
#' \code{\link{plot_prediction_LOO}}.
#'
#' @usage
#' add_ass_names(Fx, Fy, assMotif, cols)
#'
#' @param Fx,Fy two numeric vectors of \code{length(Fx)}.
#' The first vector contains reference, observed assemblage performances.
#' The second vector contains assemblage performances
#' modelled by components clustering.
#'
#' @param assMotif an integer vector of \code{length(Fx)}.
#' The vector contains the labels of each assemblages.
#' The labels are only used to plot the assemblage performances
#' that share a same assembly motif
#' with a same symbol and color.
#' If it is omitted, all points are plotted using a same symbol and color.
#' The default symbol and color is \code{"red circle"}.
#'
#' @param cols an integer vector of \code{length(Fx)}.
#' It contains the colours specific to each assembly motif
#' recorded in \code{assMotif}.
#'
#' @details
#' A given assemblage is always plotted
#' using the same symbol and color.
#'
#' @return Nothing. It is a procedure.
#'
#' @importFrom graphics text
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

add_ass_names <- function(Fx, Fy, assMotif, cols) {

  index <- which(!is.na(Fy) == TRUE)
  if (length(index) == 0) {
    stop("The vector Fy cannot be null")
  } else {
    Fx        <- Fx[index]
    Fy        <- Fy[index]
    assMotif <- assMotif[index]
  }

  thre   <- max(Fx) - (max(Fx) - min(Fx)) / 5
  index1 <- which(Fx <= thre)
  index2 <- which(Fx >  thre)

  # to the right of
  if (length(index1) != 0)
    graphics::text(x = Fx[index1], y = Fy[index1],
                   labels = names(Fx)[index1],
                   col = cols[assMotif[index1]], pos = 4)

  # to the left of
  if (length(index2) != 0)
    graphics::text(x = Fx[index2], y = Fy[index2],
                   labels = names(Fx)[index2],
                   col = cols[assMotif[index2]], pos = 2)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#' @title  Plot statistics of combinatorial analysis
#'
#' @description Take a matrix of statistical results
#'  and plot Coefficient of Determination (\code{R2}),
#'  Efficiency (\code{E})
#'  proportion of missing values,
#'  and \code{AIc} and \code{AICc},
#'  \emph{versus} the number of clusters of components.
#'
#' @usage plot_stats(mStats, main = "")
#'
#' @param mStats a matrix
#'  that contains all statistics for combinatorial analysis.
#'  The matrix is generated by the function \code{\link{compute_fit_stats}}.
#'
#' @param main a string. The string is used as graph title.
#'
#' @details The colnames of matrix are:
#' \code{"missing"}, \code{"R2cal"}, \code{"R2prd"}, \code{"pcal"},
#' \code{"pprd"}, \code{"AIC"}, \code{"AICc"}.
#'
#' @return Nothing. It is a procedure.
#'
#' @importFrom graphics plot points text abline
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot_stats <- function(mStats, main = "") {

  nbclMax <- dim(mStats)[1]

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # plot a first graph with R2cal and R2prd
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  graphics::plot(
    y    = mStats[ ,"R2cal"], ylim = c(0, 1),
    x    = seq_len(nbclMax),  xlim = c(1, nbclMax),
    type = "n", pch = 1, cex = 2, tck = 0.02, las = 1,
    bg   = "white", col = "black",
    ylab = "R2 (in black),   E (in red)",
    xlab = "number of component clusters",
    main = main )

  # plot R2 calibration
  graphics::points(
    y    = mStats[ ,"R2cal"],
    x    = seq_len(nbclMax),
    type = "b", pch = 1, cex = 2, bg = "white", col = "black")

  # plot R2 prediction with pvalues
  graphics::points(
    y    = mStats[ ,"R2prd"],
    x    = seq_len(nbclMax),
    type = "b", pch = 1, cex = 2, bg = "white", col = "red3")


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # plot a second graph with "Predicting ratio"
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  graphics::plot(
    y    = mStats[ ,"R2cal"], ylim = c(0, 1),
    x    = seq_len(nbclMax),  xlim = c(1, nbclMax),
    type = "n", pch = 1, cex = 2, tck = 0.02, las = 1,
    bg   = "white", col = "black",
    ylab = "Predicting ratio",
    xlab = "number of component clusters",
    main = main )

  # plot predictiong ratio
  graphics::points(
    y    = 1 - mStats[ ,"missing"],
    x    = seq_len(nbclMax),
    type = "b", pch = 0, cex = 2, bg = "white", col = "blue" )


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # plot a third graph with AIC and AICc
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  ylim <- c(min(mStats[ , "AIC"], mStats[ , "AICc"], na.rm = TRUE),
            max(mStats[1, "AIC"], mStats[1, "AICc"], na.rm = TRUE))

  # plot AIC
  graphics::plot(
    y    = mStats[ ,"AIC"],  ylim = ylim,
    x    = seq_len(nbclMax), xlim = c(1, nbclMax),
    type = "n", pch = 1, cex = 2, tck = 0.02, las = 1,
    bg   = "white", col = "black",
    ylab = "AIC", xlab = "number of component clusters",
    main = main )

  fct <- mStats[!is.na(mStats[ , "AIC"]), "AIC"]
  graphics::abline(v = first_argmin(fct), col = "green3", lty = "dashed")
  graphics::points(y = fct, x = seq_along(fct),
                   type = "b", pch = 1, cex = 2, bg = "white", col = "blue3" )

  # plot AICc
  graphics::plot(
    y    = mStats[ ,"AIC"],  ylim = ylim,
    x    = seq_len(nbclMax), xlim = c(1, nbclMax),
    type = "n", pch = 1, cex = 2, tck = 0.02, las = 1,
    bg   = "white", col = "black",
    ylab = "AICc", xlab = "number of component clusters",
    main = main )

  fct <- mStats[!is.na(mStats[ , "AICc"]), "AICc"]
  graphics::abline(v = first_argmin(fct), col = "green3", lty = "dashed")
  graphics::points(y = fct, x = seq_along(fct),
                   type = "b", pch = 1, cex = 2, bg = "white", col = "red3")

}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Plot Simulated or Predicted \emph{vs} Observed performances
#'
#' @description Take two vectors
#' corresponding to assemblage performances
#' modelled by component clustering model,
#' or assemblage performances predicted by cross-validation,
#' and reference, observed assemblage performances,
#' then plot modelled assemblage performances
#' \emph{versus} observed assemblage performances.
#'
#' @usage
#' plot_prediction_simple(Fx, Fy,
#'                        assMotif = rep(1, length(Fx)),
#'                        xylab     = c("Observed performance",
#'                                      "Modelled performance"),
#'                        xylim     = range(Fx),
#'                        figs      = rep(fsymbols()[1], length(Fx)),
#'                        cols      = rep(fcolours()[1], length(Fx)),
#'                        nbcl      = 1,
#'                        main      = "",
#'                        opt.mean  = "amean",
#'                        pvalue    = fpvalue())
#'
#' @param Fx,Fy two numeric vector of \code{length(Fx)}.
#' The first vector contains assemblage performances to plot on x-axis.
#' The second vector contains assemblage performances to plot on y-axis.
#'
#' @param assMotif an integer vector of \code{length(Fx)}.
#' The vector contains the labels of each assemblages.
#' The labels are only used to plot the assemblage performances
#' that share a same assembly motif
#' with a same symbol and color.
#' If it is omitted, all points are plotted using a same symbol and color.
#' The default symbol and color is \code{"red circle"}.
#'
#' @param xylab a vector of two strings.
#' The strings are the labels of \code{(x,y)} axis.
#'
#' @param xylim a vector of two numerics.
#'  The numerics are extrem values of \code{(x,y)} axis.
#'
#' @param figs,cols two integer vectors of \code{length(Fx)}.
#' They contain the symbols and colours specific to each assembly motif
#' recorded in \code{assMotif}.
#'
#' @param nbcl an integer.
#' The integer indicates the number of clusters of components.
#' The theoretical number \code{m} of assembly motif
#' is \code{m = 2^nbcl - 1}.
#' The number of clusters of components is only used
#' as information written on the graph.
#'
#' @param main a string. The string is used as graph title.
#'
#' @param opt.mean a character equals to \code{"amean"} or \code{"gmean"}.
#' Switch to arithmetic formula if \code{opt.mean = "amean"}.
#' Switch to geometric formula if \code{opt.mean = "gmean"}.
#'
#' @param pvalue a probability,
#'   used as threshold
#'   in the variance analysis. Then \code{pvalue} must be
#'   higher than \code{0} and lower than \code{1}.
#'   Groups significantly different
#'   (at a p-value < \code{pvalue}) are then indicated by differents letters
#'   on the right of boxplots.
#'
#' @details
#' All options are default values. If all options are omitted,
#' the function plot \code{Fy} \emph{vs} \code{Fx}, in red circle,
#' with labels of x-axis \code{xlab = "Observed performances"} and
#'                y-axis \code{ylab = "Simulated performances"}. \cr
#'
#' The two dashed blue lines are \code{mean(Fy)} and \code{mean(Fx)}.
#'
#' @return Nothing. It is a procedure.
#'
#' @importFrom graphics plot points text abline lines
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot_prediction_simple <- function(Fx, Fy,
                                   assMotif = rep(1, length(Fx)),
                                   xylab    = c("Observed performance",
                                                "Modelled performance"),
                                   xylim    = range(Fx),
                                   figs     = rep(fsymbols()[1], length(Fx)),
                                   cols     = rep(fcolours()[1], length(Fx)),
                                   nbcl     = 1,
                                   main     = "",
                                   opt.mean = "amean",
                                   pvalue   = fpvalue()) {

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Check the inputs
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  if (is.character(assMotif)) {
    setMot <- unique(assMotif)
    for (ass in seq_along(setMot)) assMotif[assMotif == setMot[ass]] <- ass
  }
  assMotif <- as.integer(assMotif)

  if (length(Fx) != length(Fy)) stop("length(Fx) == length(Fy)")

  index <- which(!is.na(Fy) == TRUE)
  if (length(index) == 0) {
    stop("The vector Fy cannot be null")
  } else {
    Fx       <- Fx[index]
    Fy       <- Fy[index]
    assMotif <- assMotif[index]
  }


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Plot the figure
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  graphics::plot(x    = Fx, xlab = xylab[1], xlim = xylim,
                 y    = Fy, ylab = xylab[2], ylim = xylim,
                 main = main,
                 type = "n", tck = 0.02, las = 1)

  graphics::abline(v = mean_fct(Fx, opt.mean),
                   h = mean_fct(Fy, opt.mean),
                   lty = "dashed", col = "blue")

  graphics::lines(x = xylim, y = xylim, lty = "solid", col = "red")

  graphics::points(x = Fx, y = Fy,
                   pch = figs[assMotif], col = cols[assMotif],
                   bg  = "white", cex = 2)


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Add various useful informations
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  width <- xylim[2] - xylim[1]

  posX1 <- xylim[1] + 0.75 * width
  posX2 <- xylim[1] + 1.02 * width
  #  posX3 <- xylim[1] - 0.02 * width
  #  posX4 <- xylim[1] + 0.02 * width
  posX5 <- xylim[1] + 0.25 * width

  posY1 <- xylim[1] + 0.05 * width
  posY2 <- xylim[1]
  posY3 <- xylim[1] + 1.02 * width
  posY4 <- xylim[1] + 0.98 * width
  posY5 <- xylim[1] + 0.93 * width


  # predicting ratio
  graphics::text(x = posX1, y = posY1,
                 labels = paste0("predicting = ",
                                 sum(!is.na(Fy)), " / ", length(Fy)),
                 col = "red")

  # R2 value
  if (length(index) > 1)
    graphics::text(x = posX1, y = posY2,
                   labels = paste0("R2 = ", signif(R2mse(Fy, Fx), digits = 3)),
                   col = "red")

  # Number of clusters
  graphics::text(paste("Nb of component clusters =", nbcl, sep = " "),
                 x = posX5, y = posY4, col = "red")

  tmp <- length(unique(assMotif))
  if (tmp > 1)
    graphics::text(paste0("Nb of observed assembly motifs = ",
                          tmp, " / ", 2 ^ nbcl - 1),
                   x = posX5, y = posY5, col = "red")


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Variance analysis
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  setMot <- sort(unique(assMotif))
  if (length(setMot) > 1) {
    graphics::axis(labels = paste0("p < ", pvalue),
                   at     = max(Fx),
                   col    = "black",
                   side   = 3,
                   tick   = FALSE, las = 1, font = 1 )

    graphics::text(x = posX2, y = posY3, labels = "pred", col = "black")

    test <- test_posthoc(Fy, assMotif, pvalue)
    if (is.list(test))
      for (mot in seq_along(setMot)) {
        motif <- setMot[mot]
        index <- which(test$motif == motif)
        graphics::text(x = posX2, y = test[index, "mean"],
                       labels = as.character(test[index, "group"]),
                       col = cols[motif], font = 3)
      }
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Plot Simulated and Predicted \emph{vs} Observed performances
#'
#' @description Take three vectors
#' corresponding to assemblage performances modelled,
#' predicted by cross-validation and observed,
#' then plot modelled assemblage performances
#' \emph{versus} observed assemblage performances.
#' The error induced by the cross-validation
#' is drawn as a line linking modelled and predicted assemblage performances.
#'
#' @usage
#' plot_prediction_LOO(Fx, Fy, Fz,
#'                     assMotif = rep(1, length(Fx)),
#'                     xylab     = c("Observed performance",
#'                                   "Modelled performance"),
#'                     xylim     = range(Fx),
#'                     figs = rep(fsymbols()[1], length(Fx)),
#'                     cols = rep(fcolours()[1], length(Fx)),
#'                     nbcl      = 1,
#'                     main      = "",
#'                     opt.mean  = "amean",
#'                     pvalue    = fpvalue())
#'
#' @param Fx,Fy,Fz three numeric vectors of \code{length(Fx)}.
#' The first vector contains reference, observed assemblage performances.
#' The second vector contains assemblage performances
#' modelled by component clustering.
#' The third vector contains assemblage performances
#' predicted by cross-validation.
#'
#' @param assMotif an integer vector.
#' The vector contains the labels of each assemblages.
#' The labels are only used to plot the assemblage performances
#' that share a same assembly motif
#' with a same symbol and color.
#' If it is omitted, all points are plotted using a same symbol and color.
#' The default symbol and color is \code{"red circle"}.
#'
#' @param figs,cols two integer vectors of \code{length(Fx)}.
#' They contain the symbols and colours specific to each assembly motif
#' recorded in \code{assMotif}.
#'
#' @param nbcl an integer.
#' The integer indicates the number of clusters of components.
#' The theoretical number \code{m} of assembly motif
#' is \code{m = 2^nbcl - 1}.
#' The number of clusters of components is only used
#' as information written on the graph.
#'
#' @param main a string. The string is used as graph title.
#'
#' @param xylab a vector of two strings.
#' The strings are the labels of the axis.
#'
#' @param xylim a vector of two numerics.
#'  The numerics are extrem values of the two axis.
#'
#' @param opt.mean a character equals to \code{"amean"} or \code{"gmean"}.
#' Switch to arithmetic formula if \code{opt.mean = "amean"}.
#' Switch to geometric formula if \code{opt.mean = "gmean"}.
#'
#' @param pvalue a probability,
#'   used as threshold
#'   in the variance analysis. Then \code{pvalue} must be
#'   higher than \code{0} and lower than \code{1}.
#'   Groups significantly different
#'   (at a p-value < \code{pvalue}) are then indicated by differents letters
#'   on the right of boxplots.
#'
#' @details
#' All options are default values. If all options are omitted,
#' the function plot \code{Fy vs Fx},
#' in red circle with \code{"Simuled performances"}
#' \emph{vs} \code{"Observed performances"}
#' as axis labels.\cr
#'
#' The two dashed blue lines are \code{mean(Fy)} and \code{mean(Fx)}.
#'
#' @return Nothing. It is a procedure.
#'
#' @importFrom graphics plot points text abline lines
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot_prediction_LOO <- function(Fx, Fy, Fz,
                                assMotif = rep(1, length(Fx)),
                                xylab     = c("Observed performance",
                                              "Modelled performance"),
                                xylim     = range(Fx),
                                figs      = rep(fsymbols()[1], length(Fx)),
                                cols      = rep(fcolours()[1], length(Fx)),
                                nbcl      = 1,
                                main      = "",
                                opt.mean  = "amean",
                                pvalue    = fpvalue()) {

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Check the inputs
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  if (is.character(assMotif)) {
    setMot <- unique(assMotif)
    for (ass in seq_along(setMot)) assMotif[assMotif == setMot[ass]] <- ass
  }
  assMotif <- as.integer(assMotif)

  if ((length(Fx) != length(Fy)) | (length(Fx) != length(Fz)) |
      (length(Fy) != length(Fz)))
    stop("length(Fx) == length(Fy) == length(Fz)")

  index <- which(!is.na(Fz) == TRUE)
  if (length(index) == 0) {
    stop("The vector Fz cannot be null")
  } else {
    Fx        <- Fx[index]
    Fy        <- Fy[index]
    Fz        <- Fz[index]
    assMotif <- assMotif[index]
  }


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Plot the graph
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  graphics::plot(
    x = Fx, xlab = xylab[1], xlim = xylim,
    y = Fy, ylab = xylab[2], ylim = xylim,
    main = main,
    las = 1, type = "n", tck = 0.02)

  for (elt in seq_along(Fx))
    graphics::lines(x = c(Fx[elt], Fx[elt]), y = c(Fy[elt], Fz[elt]),
                    col = cols[assMotif[elt]], lty = "solid")

  graphics::abline(v = mean_fct(Fx, opt.mean),
                   h = mean_fct(Fy, opt.mean),
                   lty = "dashed", col = "blue")

  graphics::lines(x = xylim, y = xylim, lty = "solid", col = "red")

  graphics::points(x   = Fx, y = Fy,
                   pch = figs[assMotif],
                   col = cols[assMotif],
                   bg  = "white", cex = 2)


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Add various useful informations
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  width <- xylim[2] - xylim[1]

  posX1 <- xylim[1] + 0.75 * width
  posX2 <- xylim[1] + 1.02 * width
  posX3 <- xylim[1] - 0.02 * width
  #  posX4 <- xylim[1] + 0.02 * width
  posX5 <- xylim[1] + 0.25 * width

  posY1 <- xylim[1] + 0.05 * width
  posY2 <- xylim[1]
  posY3 <- xylim[1] + 1.02 * width
  posY4 <- xylim[1] + 0.98 * width
  posY5 <- xylim[1] + 0.93 * width


  if ( (length(index) > 1) & (length(unique(assMotif)) > 1) ) {

    # Predicting assemblages

    graphics::text(x = posX1, y = posY1,
                   labels = paste0("predicting ratio = ",
                                   sum(!is.na(Fy)), " / ", length(Fy)),
                   col = "red")

    #  Number of clusters and R2 value

    graphics::text(
      x = posX1, y = posY2,
      labels = paste0("R2 = ",  signif(R2mse(Fy, Fx), digits = 3),
                      "    E = ", signif(R2mse(Fz, Fx), digits = 3),
                      "    E/R2 = ", signif(R2mse(Fz, Fx) /
                                              R2mse(Fy, Fx), digits = 3)),
      col = "red")

    graphics::text(paste0("Nb of component clusters = ", nbcl),
                   x = posX5, y = posY4, col = "red")

    tmp <- length(unique(assMotif))
    if (tmp > 1) graphics::text(paste("Nb of observed assembly motifs =",
                                      tmp, "/", (2^nbcl - 1), sep = " "),
                                x = posX5, y = posY5, col = "red")
  }




  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Variance analysis
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  setMot <- sort(unique(assMotif))
  if (length(setMot) > 1) {

    graphics::axis(labels = rep(paste0("p < ", pvalue), 2),
                   at     = range(Fx),
                   col    = "black",
                   side   = 3,
                   tick   = FALSE, las = 1, font = 1 )

    graphics::text(x = posX3, y = posY3, labels = "cal", col = "black")
    test <- test_posthoc(Fy, assMotif, pvalue)
    if (is.list(test))
      for (mot in seq_along(setMot)) {
        motif <- setMot[mot]
        index <- which(test$motif == motif)
        graphics::text(x = posX3, y = test[index, "mean"],
                       labels = as.character(test[index, "group"]),
                       col = cols[motif], font = 3)
      }

    graphics::text(x = posX2, y = posY3, labels = "prd", col = "black")
    test <- test_posthoc(Fz, assMotif, pvalue)
    if (is.list(test))
      for (mot in seq_along(setMot)) {
        motif <- setMot[mot]
        index <- which(test$motif == motif)
        graphics::text(x = posX2, y = test[index, "mean"],
                       labels = as.character(test[index, "group"]),
                       col = cols[motif], font = 3)
      }
  }
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Plot Option.motif                                                        ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#  nbline = number of graphics::lines by page

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Plot a list of strings on several pages
#'
#' @description Take a list of strings,
#' then divide it in several units of same length a a page.
#' A page is defined by a number of lines.
#'
#' @usage
#' plot_by_page(tab, nbline = 25)
#'
#' @param tab a list of strings to plot.
#'
#' @param nbline an integer,
#' that corresponds to a number of lines to plot by page.
#'
#' @details None.
#'
#' @return Nothing. It is a procedure.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot_by_page <- function(tab, nbline = 25) {

  page <- 0
  ltab <- length(tab)
  while (nbline * page < ltab) page <- page + 1

  for (p in 1:page) {

    tmp <- tab[intersect((p - 1) * nbline + 1:nbline, seq_along(tab))]

    graphics::plot(x = seq_len(nbline), xlab = "",
                   y = seq_len(nbline), ylab = "",
                   type = "n", axes = FALSE)

    graphics::axis(side = 4, at = seq_along(tmp), rev(tmp),
                   pos = 1, tick = FALSE, las = 2, col = "black",
                   family = "mono")
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Sort assembly motifs
#'
#' @description Sort assembly motifs,
#'  by decreasing or increasing mean performances
#'
#' @usage
#' sort_motifs(fobs, assMotif, assNames,
#'             pvalue  = fpvalue(),
#'             opt.dec = TRUE   )
#'
#' @inheritParams plot_motifs_box
#'
#' @param assNames a vector of strings of \code{length(fobs)}.
#' The vector contains the names of each assemblage.
#'
#' @param pvalue a probability,
#'   used as threshold
#'   in the variance analysis. Then \code{pvalue} must be
#'   higher than \code{0} and lower than \code{1}.
#'   Groups significantly different
#'   (at a p-value < \code{pvalue}) are then indicated by differents letters
#'   on the right of boxplots.
#'
#' @param opt.dec a logical.
#' If \code{opt.dec = TRUE},
#' assembly motifs are sorted by decreasing mean performances.
#' If \code{opt.dec = FALSE},
#' assembly motifs are sorted by increasing mean performances.
#'
#' @details None.
#'
#' @return Return a table
#' containing statistical properties of assembly motifs,
#' sorted by increasing or decreasing mean performances.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

sort_motifs <- function(fobs, assMotif, assNames,
                        pvalue  = fpvalue(),
                        opt.dec = TRUE     ) {

  # opt.dec = TRUE by default because test_posthoc() is "decreasing"

  sres       <- test_posthoc(x = fobs, clusters = assMotif, pvalue = pvalue)
  sres$motif <- as.integer(sres$motif)

  indAss <- unique(assMotif)
  indNom <- unique(assNames)

  snoms  <- character(length(sres$motif))

  for (mot in seq_along(sres$motif)) {
    motif      <- sres$motif[mot]
    snoms[mot] <- indNom[which(indAss == motif)]
  }

  sres <- cbind(snoms, sres)
  colnames(sres)[1] <- c("nom")

  if (opt.dec == FALSE)
    for (j in seq_len(dim(sres)[2])) sres[ , j] <- rev(sres[ , j])

  rownames(sres) <- c(1:dim(sres)[1])

  return(sres)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Plot histogram of assemblage performances by assembly motif
#'
#' @description Plot the histogram of performances of all assemblage,
#' and the density curves of performances sorted by assembly motif
#'
#' @usage
#' plot_motifs_histo(fobs, assMotif,
#'                   nbar    = 0,
#'                   main     = "",
#'                   opt.mean = "amean")
#'
#' @param fobs numeric vector.
#' The vector contains observed assemblage performances.
#'
#' @param assMotif an integer vector of \code{length(fobs)}.
#' The vector contains the labels of each assemblages.
#' The labels are only used to plot the assemblage performances
#' that share a same assembly motif
#' with a same symbol and color.
#'
#' @param nbar an integer.
#' It is the number of wished bars of histogram.
#'
#' @param main a string. The string is used as graph title.
#'
#' @param opt.mean a character equals to \code{"amean"} or \code{"gmean"}.
#' Switch to arithmetic formula if \code{opt.mean = "amean"}.
#' Switch to geometric formula if \code{opt.mean = "gmean"}.
#'
#' @details A colour is assigned to each assembl motif.
#' The density curve of each assembly motif
#' and its mean are plotted.
#' The location of each observation is reported on the graph bottom.
#'
#' @return Nothing. It is a procedure.
#'
#' @importFrom stats density
#'
#' @importFrom graphics hist lines abline rug
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot_motifs_histo <- function(fobs, assMotif,
                              nbar     = 0,
                              main      = "",
                              opt.mean = "amean") {

  couleurs <- extend_vector(fcolours(), length(fobs))

  ymin <- floor(min(fobs))
  ymax <- 1 + floor(max(fobs))

  if (nbar == 0) nbar <- 20
  breaks <- seq(ymin, ymax, by = (ymax - ymin) / nbar)

  graphics::hist(fobs, xlab = "Observed performances", ylab = "Density",
                 main = main, prob = TRUE, las = 1, breaks = breaks)

  setClu <- unique(assMotif)
  for (clu in seq_along(setClu)) {

    clust <- setClu[clu]
    setY  <- fobs[assMotif == clust]
    if (length(setY) > 1) {

      tmp <- stats::density(setY)
      graphics::lines(x = tmp$x, y = tmp$y * length(setY) / length(fobs),
                      col = couleurs[clust])

      graphics::abline(v   = mean_fct(setY, opt.mean),
                       col = couleurs[clust])

      graphics::rug(setY, col = couleurs[clust])
    }
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Plot reference graphs for checking that plot of sorted data are right
#'
#' @description Plot two reference graphs for checking
#' that plot of sorted data are right:
#'   the non-sorted assembly motifs
#'   and assembly motifs sorted by decreasing mean observed performances.
#'
#' @usage
#' plot_motifs_infos(res, opt.sort = "performance")
#'
#' @param res the result of a functional clustering
#'  obtained using the function \code{\link{fclust}}.
#'
#' @param opt.sort a string,
#'   that specifies the way for sorting the motifs.
#'   It can be \code{"performance"} or \code{"name"},
#'   indicating a sorting by motif performances,
#'   or a sorting by motif names.
#'
#' @details This function was useful only for setting up the R-code.
#' It is now useful for users to check the resulst,
#' and be sure to the goodness of plotting ...
#' The written values are (form up to bottom):
#' names of assembly motifs, their effectifs,
#' the mean observed performance, the order
#' and the symbols (colour x symbol) systematically associated
#' with each assembly motif in all plots produced
#' by the package \code{functClust}.
#'
#' @importFrom graphics plot axis text
#'
#' @return Nothing. It is a procedure.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot_motifs_infos <- function(res, opt.sort = "performance") {

  if ( (opt.sort != "performance") & (opt.sort != "name") )
      stop("'opt.sort' == 'performance' or 'name'")

  cols <- extend_vector(fcolours(), res$nbAss)
  pchs <- extend_vector(fsymbols(), res$nbAss)

  setMot  <- unique(res$mMotifs[res$nbOpt, ])
  index   <- seq_along(setMot)
  nameMot <- unique(name_motifs(cut_ftree(res$tree.II, res$nbOpt), res$mOccur))

  ttt   <- integer(length(setMot))
  vmean <- vsd <- numeric(length(setMot))
  for (mot in seq_along(setMot)) {
    indMot     <- (res$mMotifs[res$nbOpt, ] == setMot[mot])
    ttt[mot]   <- sum(indMot)
    vmean[mot] <- amean(res$fobs[indMot])
    vsd[mot]   <- asd(res$fobs[indMot])
  }


  # References symbols and colours sorted numerically

  if (opt.sort == "performance") {

    order_num   <- sort(vmean, decreasing = FALSE, index.return = TRUE)

    graphics::plot(main = "Assembly motifs sorted by performances",
                   x = c(index, index + 1), xlim = c(0, 7), xlab = "",
                   y = c(index, index + 1), ylab = "",
                   type = "n", axes = FALSE)

    graphics::points(x = jitter(rep(2, length(index)), factor = 5),
                     xlim = c(0, 7),
                     y = index, ylim = c(0, res$nbElt),
                     pch = pchs[setMot][order_num$ix],
                     col = cols[setMot][order_num$ix],
                     cex = 3)

    graphics::axis(labels = nameMot[index][order_num$ix],
                   tick = FALSE, las = 1, side = 2,
                   at = index, pos = 1)

    graphics::text(labels = ttt[order_num$ix], y = index, x = 3 )

    graphics::text(labels = signif(vmean[order_num$ix], digits = 4),
                   y = index, x = 4 )

    graphics::text(labels = signif(vsd[order_num$ix], digits = 4),
                   y = index, x = 5 )

    graphics::text(labels = rev(index_inturn(order_num$ix)[order_num$ix]),
                   y = index, x = 6 )

    graphics::text(
      labels = c("motif", "symbol", "effectif", "mean", "sd", "order"),
      y = length(index) + 1, x = 1:6 )
  }


  # References symbols and colours sorted alphabetically

  if (opt.sort == "name") {

    order_alpha <- sort(nameMot, decreasing = TRUE, index.return = TRUE)

    graphics::plot(main = "Assembly motifs sorted by names",
                   x = c(index, index + 1), xlim = c(0, 7), xlab = "",
                   y = c(index, index + 1), ylab = "",
                   type = "n", axes = FALSE)

    graphics::points(x = jitter(rep(2, length(index)), factor = 5),
                     xlim = c(0, 7),
                     y = index, ylim = c(0, res$nbElt),
                     pch = pchs[setMot][order_alpha$ix],
                     col = cols[setMot][order_alpha$ix],
                     cex = 3)

    graphics::axis(labels = nameMot[index][order_alpha$ix],
                   tick = FALSE, las = 1, side = 2,
                   at = index, pos = 1)

    graphics::text(labels = ttt[order_alpha$ix], y = index, x = 3 )

    graphics::text(labels = signif(vmean[order_alpha$ix], digits = 4),
                   y = index, x = 4 )

    graphics::text(labels = signif(vsd[order_alpha$ix], digits = 4),
                   y = index, x = 5 )

    graphics::text(labels = rev(index_inturn(order_alpha$ix)[order_alpha$ix]),
                   y = index, x = 6 )

    graphics::text(
      labels = c("motif", "symbol", "effectif", "mean", "sd", "order"),
      y = length(index) + 1, x = 1:6 )
  }
}



##xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Write the assemblages belonging to each assembly motif
#'
#' @description Write the assemblages belonging to each assembly motif
#' on several pages.
#'
#' @usage
#' plot_motifs_content(fobs, assMotif, assNames,
#'                     opt.sort = TRUE,
#'                     sepChar  = ", ",
#'                     nbline   = 25  )
#'
#' @inheritParams plot_motifs_box
#'
#' @param assNames a vector of strings of \code{length(fobs)}.
#' The vector contains the names of each assemblage.
#'
#' @param opt.sort a logical.
#' If \code{opt.sort = TRUE}, the names of assemblages
#' are sorted alphabetically inside each motif.
#'
#' @param sepChar a string.
#' The string is a separator between each element of list of strings.
#'
#' @param nbline an integer.
#' The number of lines to plot.
#'
#' @details None.
#'
#' @return Nothing. It is a procedure.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot_motifs_content <- function(fobs, assMotif, assNames,
                                opt.sort = TRUE,
                                sepChar  = ", ",  nbline   = 25  ) {

  tmp     <- sort_motifs(fobs, assMotif, assNames)
  tmp$nom <- as.character(tmp$nom)

  # compute the number of necessary columns
  ttp <- list()
  for (i in seq_along(tmp$nom))
    ttp[[i]] <- names(fobs)[which(assNames == tmp$nom[i])]

  if (opt.sort == TRUE)
    for (i in seq_along(tmp$nom)) ttp[[i]] <- sort(ttp[[i]])

  # adjust the length of rownames
  lmax <- max(nchar(tmp$nom))
  for (i in seq_along(tmp$nom)) {
    a <- tmp$nom[i]
    while (nchar(a) < lmax) a <- paste0(" ", a)
    names(ttp[[i]])[1] <- a
  }

  blanc <- NULL
  for (k in 1:(lmax + 5)) blanc <- paste0(blanc, " ")

  # concat the matrix by row
  tab <- c("", "Motif content:", "")
  for (i in seq_along(tmp$nom)) {

    tpp    <- concat_by_line(ttp[[i]], sepChar = sepChar)
    tpp[1] <- paste0(names(ttp[[i]])[1], " = { ", tpp[1])
    if (length(tpp) > 1)
      for (j in 2:length(tpp)) tpp[j] <- paste0(blanc, tpp[j])
    tpp[length(tpp)] <- paste0(tpp[length(tpp)], " }")
#    tpp[length(tpp)] <- paste0(delstr_end(tpp[length(tpp)], 1), " }")

    tab  <- c(tab, tpp)
  }
  tab  <- c(tab, "")

  # plot the legend
  plot_by_page(tab, nbline)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Plot performances of assembly motifs
#'
#' @description Plot performances of assembly motifs as boxplots,
#' possibly horizontally or vertically,
#' by decreasing or increasing values.
#'
#' @usage
#' plot_motifs_box(fobs, assMotif, sres,
#'                 ordre    = NULL,
#'                 xylim    = range(fobs),
#'                 main     = "",
#'                 xlabel   = "",
#'                 pvalue   = fpvalue(),
#'                 opt.mean = "amean",
#'                 opt.hor  = TRUE)
#'
#' @param fobs numeric vector.
#' The vector contains assemblage performances.
#'
#' @param assMotif an integer vector of \code{length(fobs)}.
#' The vector contains the labels of each assemblages.
#' The labels are used to sorted the assemblages
#' that share a same assembly motif.
#'
#' @param sres a table generated by the function \code{\link{sort_motifs}}.
#' The table contains statistical properties of assembly motifs,
#' sorted by increasing or decreasing mean performances.
#'
#' @param ordre an integer vector of \code{length(fobs)}.
#' The vector contains labels for each assemblages,
#' that determine the order for assembly motifs.
#'
#' @param xylim a vector of two numerics.
#'  The numerics are extrem values of performance \code{y-axis}.
#'
#' @param main a string. The string is used as graph title.
#'
#' @param xlabel a string. The string is used as label for abscissa (xlab).
#'
#' @param pvalue a probability,
#'   used as threshold
#'   in the variance analysis. Then \code{pvalue} must be
#'   higher than \code{0} and lower than \code{1}.
#'   Groups significantly different
#'   (at a p-value < \code{pvalue}) are then indicated by differents letters
#'   on the right of boxplots.
#'
#' @param opt.mean a character equals to \code{"amean"} or \code{"gmean"}.
#' Switch to arithmetic formula if \code{opt.mean = "amean"}.
#' Switch to geometric formula if \code{opt.mean = "gmean"}.
#'
#' @param opt.hor a logical.
#' If \code{opt.hor = TRUE},
#' boxplot is plotted as horizontal boxes
#' alongside a vertical y-axis of performances.
#' If \code{opt.hor = FALSE},
#' boxplot is plotted as vertical boxes
#' alongside a horizontal y-axis of performances.
#'
#' @details None.
#'
#' @return Nothing. It is a procedure.
#'
#' @importFrom graphics boxplot abline points text axis title
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot_motifs_box <- function(fobs, assMotif, sres,
                            ordre    = NULL,
                            xylim    = range(fobs),
                            main     = "",
                            xlabel   = "",
                            pvalue   = fpvalue(),
                            opt.mean = "amean",
                            opt.hor  = TRUE) {

  figures  <- extend_vector(fsymbols(), length(fobs))
  couleurs <- extend_vector(fcolours(), length(fobs))

  # re-organize the motifs to appear as continuously increasing or decreasing
  ass.sort <- integer(length(fobs))
  if (is.null(ordre)) {

    for (mot in seq_along(sres$motif)) {
      motif <- sres$motif[mot]
      ass.sort[which(assMotif == motif)] <- mot
    }

  } else {

    index <- integer(length(ordre))
    for (mot in seq_along(ordre)) {
      motif      <- ordre[mot]
      index[mot] <- which(sres$motif == motif)
      ass.sort[which(assMotif == motif)] <- mot
    }
    sres <- sres[index, ]
  }


  # horizontal plotting of boxes is the most complete possibly
  if (opt.hor == TRUE) {

    graphics::boxplot(fobs ~ ass.sort,
                      names = sres$nom,
                      ylab  = "Assembly motif",
                      xlab  = xlabel,
                      ylim  = xylim,
                      horizontal = opt.hor,
                      tck = 0.02, las = 1)

    # plot general and specific means
    graphics::abline(v = mean_fct(fobs, opt.mean), col = "blue", lty = "dashed")

    graphics::points(y = seq_along(sres$motif), x = sres$mean,
                     col = couleurs[sres$motif],
                     pch = figures[sres$motif],
                     cex = 2, bg = "white"  )

    # plot additional informations
    graphics::text(y = seq_along(sres$motif), x = xylim[1],
                   pos    = 4,    # on right
                   col    = couleurs[sres$motif],
                   labels = as.character(sres$number) )

    graphics::text(y = seq_along(sres$motif), x = xylim[2],
                   pos    = 2,   # on left
                   col    = couleurs[sres$motif],
                   labels = as.character(sres$group),
                   font   = 3)

    graphics::axis(labels = c("size", paste0("p < ", pvalue)),
                   at     = xylim,
                   col    = "black",
                   side   = 3,
                   tick   = FALSE, las = 1, font = 1 )

  } else {

    # vertical plotting of boxes

    graphics::boxplot(fobs ~ ass.sort,
                      names = sres$nom,
                      xlab  = "Assembly motif",
                      ylab  = xlabel,
                      ylim  = xylim,
                      horizontal = opt.hor, las = 2)

    # plot general and specific means
    graphics::abline(h = mean_fct(fobs, opt.mean), col = "blue", lty = "dashed")

    graphics::points(x = seq_along(sres$motif), y = sres$mean,
                     col = couleurs[sres$motif],
                     pch = figures[sres$motif],
                     cex = 2, bg = "white"  )

    # plot additional informations
    graphics::text(x = 0, y = xylim, pos = 4,
                   col    = "black",
                   labels = c("size", paste0("p < ", pvalue)) )

    graphics::text(x = seq_along(sres$motif), y = xylim[1], pos = 4,
                   col    = couleurs[sres$motif],
                   labels = as.character(sres$number) )

    graphics::text(x = seq_along(sres$motif), y = xylim[2], pos = 4,
                   col    = couleurs[sres$motif],
                   labels = as.character(sres$group),
                   font   = 3)
  }

  graphics::title(main)
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Plot Option.component                                                    ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Sort assembly centroids by decreasing or increasing mean performances
#'
#' @description Plot
#'
#' @usage
#' sort_components(fobs, mOccur, cols,
#'                 pvalue  = fpvalue(),
#'                 opt.dec = TRUE  )
#'
#' @param fobs numeric vector.
#' The vector contains observed assemblage performances.
#'
#' @param mOccur a matrix of occurrence of components within the assemblages.
#' \code{dim(mOccur)[1] = length(fobs)}.
#' \code{dim(mOccur)[2]} equals to the total number of components.
#'
#' @param cols an integer vector of \code{length(fobs)}.
#' The vector contains a colour label for each assemblages.
#'
#' @param pvalue a probability,
#'   used as threshold
#'   in the variance analysis. Then \code{pvalue} must be
#'   higher than \code{0} and lower than \code{1}.
#'   Groups significantly different
#'   (at a p-value < \code{pvalue}) are then indicated by differents letters
#'   on the right of boxplots.
#'
#' @param opt.dec a logical.
#' If \code{opt.dec = TRUE},
#' assembly motifs are sorted by decreasing mean performances.
#' If \code{opt.dec = FALSE},
#' assembly motifs are sorted by increasing mean performances.
#'
#' @details None.
#'
#' @return Return a table
#' containing statistical properties of assembly centroids,
#' sorted by increasing or decreasing mean performances.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

sort_components <- function(fobs, mOccur, cols,
                            pvalue  = fpvalue(),
                            opt.dec = TRUE  ) {

  nbElt  <- dim(mOccur)[2]
  if (length(cols) < nbElt) cols <- extend_vector(cols, length(fobs))

  nbAss  <- dim(mOccur)[1]
  index  <- integer(nbAss)

  sres   <- NULL
  for (elt in seq_len(nbElt)) {

    index[ ] <- 0
    index[which(mOccur[ , elt] == 1)] <- elt
    tmp <- test_posthoc(x = fobs, clusters = index, pvalue = pvalue)

    if (tmp[1, "motif"] == elt) {
      stmp <- tmp[1, ]
      stmp$group <- ifelse(tmp[2, "group"] == "a", "b", "a")
    } else {
      stmp <- tmp[2, ]
      stmp$group <- ifelse(tmp[2, "group"] == "a", "b", "c")
    }

    sres <- rbind(sres, stmp)
  }

  sres$motif <- as.integer(sres$motif)
  nomElt <- colnames(mOccur)
  coll   <- cols
  sres   <- cbind(nomElt, coll, sres)
  colnames(sres)[1:2] <- c("nom", "color")

  tmp  <- sort(sres$mean, index.return = TRUE)
  sres <- sres[tmp$ix, ]

  if (opt.dec == TRUE) sres <- reverse_table(sres)

  rownames(sres) <- c(1:(dim(sres)[1]))

  ind.inturn <- index_inturn(sres$motif)
  sres       <- cbind(ind.inturn, sres)

  return(sres)
}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Write the components belonging to each functional group
#'
#' @description Write the components belonging to each functional group,
#' on several pages if need.
#'
#' @usage
#' plot_components_content(mOccur,
#'                         ordre   = NULL,
#'                         sepChar = ", ",
#'                         nbline  = 25 )
#'
#' @param mOccur a matrix of occurrence of components within the assemblages.
#' \code{dim(mOccur)[1] = length(fobs)}.
#' \code{dim(mOccur)[2]} equals to the total number of components.
#'
#' @param ordre a logical.
#'
#' @param sepChar a string.
#' The string is a separator between each element of list of strings.
#'
#' @param nbline an integer,
#' that corresponds to a number of lines to plot by page.
#'
#' @details None.
#'
#' @return Nothing. It is a procedure.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot_components_content <- function(mOccur,
                                    ordre  = NULL,
                                    sepChar = ", ",  nbline = 25) {

  ordre  <- rev(ordre)
  setClu <- colnames(mOccur)
  nbrow  <- length(setClu)

  # compute the number of necessary columns
  tmp    <- list()
  for (i in seq_len(nbrow))
    tmp[[i]] <- rownames(mOccur)[mOccur[ , setClu[ordre[i]]] == 1]

  # adjust the length of rownames
  lmax <- max(nchar(setClu))
  for (i in seq_len(nbrow)) {
    a <- setClu[ordre[i]]
    while (nchar(a) < lmax) a <- paste(" ", a, sep = "")
    names(tmp[[i]])[1] <- a
  }

  blanc <- NULL
  for (k in 1:(lmax + 5)) blanc <- paste(blanc, " ", sep = "")

  # concat the matrix by row
  tab <- c("Content in assemblages:", "")
  for (i in seq_len(nbrow)) {

    tpp    <- concat_by_line(tmp[[i]], sepChar = sepChar)
    tpp[1] <- paste0(names(tmp[[i]])[1], " = { ", tpp[1])
    if (length(tpp) > 1)
      for (j in 2:length(tpp)) tpp[j] <- paste0(blanc, tpp[j])
    tpp[length(tpp)] <- paste0(tpp[length(tpp)], " }")

    tab  <- c(tab, tpp)
  }

  # plot the legend
  plot_by_page(tab, nbline)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Plot mean performances of assemblages
#' that contain a given component
#'
#' @description Plot performances of assemblages
#' that contain a given component
#' as boxplots,
#' possibly horizontally or vertically,
#' by decreasing or increasing values
#' or in the same order as the tree clustering.
#'
#' @usage
#' plot_components_box(fobs, mOccur, sres,
#'                     ordre = NULL, ylim = range(fobs),
#'                     elt.wdw = seq_len(dim(mOccur)[2]),
#'                     main     = "",
#'                     pvalue   = fpvalue(),
#'                     opt.mean = "amean", opt.hor = TRUE  )
#'
#' @inheritParams fit_ftree
#'
#' @param sres a table generated by the function \code{\link{sort_motifs}}.
#' The table contains statistical properties of assembly motifs,
#' sorted by increasing or decreasing mean performances.
#'
#' @param ordre an integer vector of \code{length(fobs)}.
#' The vector contains labels for each assemblages,
#' that determine the order for assembly motifs.
#'
#' @param ylim a vector of two numerics.
#'  The numerics are extrem values of performance \code{y-axis}.
#'
#' @param elt.wdw a vector of integers.
#' It indicates the components to plot.
#'
#' @param main a string. The string is used as graph title.
#'
#' @param pvalue a probability,
#'   used as threshold
#'   in the variance analysis. Then \code{pvalue} must be
#'   higher than \code{0} and lower than \code{1}.
#'   Groups significantly different
#'   (at a p-value < \code{pvalue}) are then indicated by differents letters
#'   on the right of boxplots.
#'
#' @param opt.mean a character equals to \code{"amean"} or \code{"gmean"}.
#' Switch to arithmetic formula if \code{opt.mean = "amean"}.
#' Switch to geometric formula if \code{opt.mean = "gmean"}.
#'
#' @param opt.hor a logical.
#' If \code{opt.hor = TRUE},
#' boxplot is plotted as horizontal boxes
#' alongside a vertical y-axis of performances.
#' If \code{opt.hor = FALSE},
#' boxplot is plotted as vertical boxes
#' alongside a horizontal y-axis of performances.
#'
#' @details None.
#'
#' @return Nothing. It is a procedure.
#'
#' @importFrom graphics boxplot abline points text axis title
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot_components_box <- function(fobs, mOccur, sres,
                                ordre    = NULL,
                                ylim     = range(fobs),
                                elt.wdw  = seq_len(dim(mOccur)[2]),

                                main     = "",
                                pvalue   = fpvalue(),
                                opt.mean = "amean",
                                opt.hor  = TRUE  ) {

  nbElt  <- dim(mOccur)[2]

  # re-organize the motifs to appear as continuously increasing or decreasing
  fctElt <- motElt <- NULL
  for (elt in seq_len(nbElt)) {
    component <- sres$motif[elt]
    indElt    <- which(mOccur[ , component] == 1)
    fctElt    <- c(fctElt, fobs[indElt])
    motElt    <- c(motElt, rep(elt, length(indElt)))
  }
  elt.sort <- motElt


  # sort the sres-table according input order
  fres <- sres
  if (length(ordre) != 0) {
    index <- integer(length(ordre))
    for (elt in seq_along(ordre)) {
      component  <- ordre[elt]
      index[elt] <- which(sres$motif == component)
    }
    fres <- sres[index, ]
  }


  # horizontal plotting of boxes is the most complete possibly
  indxy <- compact_index(elt.wdw)

  if (opt.hor == TRUE) {

    graphics::boxplot(fctElt ~ elt.sort,
                      subset = (elt.sort %in% elt.wdw),
                      names  = fres$nom[elt.wdw],
                      ylab = "Set of assemblages containing a given component",
                      xlab   = "Observed performance", ylim = ylim,
                      horizontal = opt.hor,
                      tck = 0.02, las = 1)

    graphics::abline(v = mean_fct(fobs, opt.mean), col = "blue", lty = "dashed")

    graphics::points(y   = indxy,
                     x   = fres$mean[elt.wdw],
                     col = as.character(fres$color[elt.wdw]),
                     pch = 0, cex = 2, bg = "white"  )

    graphics::axis(labels = c("size", paste0("p < ", pvalue)),
                   at     = range(fobs),
                   col    = "black",
                   side   = 3,
                   tick   = FALSE, las = 1, font = 1 )

    graphics::text(labels = as.character(fres$number[elt.wdw]),
                   y      = indxy,
                   x      = min(fobs),
                   col    = as.character(fres$color[elt.wdw]),
                   pos    = 4  )       # on right

    graphics::text(labels = as.character(fres$group[elt.wdw]),
                   y      = indxy,
                   x      = max(fobs),
                   col    = as.character(fres$color[elt.wdw]),
                   font   = 3, pos = 2) # on left

  } else {

    # vertical plotting of boxes

    graphics::boxplot(fctElt ~ elt.sort,
                      subset = (elt.sort %in% elt.wdw),
                      names = fres$nom[elt.wdw],
                      xlab = "Set of assemblages containing a given component",
                      ylab = "Observed performance", ylim = ylim,
                      horizontal = opt.hor,
                      tck = 0.02, las = 2)

    graphics::abline(h = mean_fct(fobs, opt.mean), col = "blue", lty = "dashed")

    graphics::points(x   = indxy,
                     y   = fres$mean[elt.wdw],
                     col = as.character(fres$color[elt.wdw]),
                     pch = 0, cex = 2, bg = "white"  )

    graphics::text(x = 0, y = range(fobs), pos = 4,
                   col    = "black",
                   labels = c("size", paste0("p < ", pvalue)) )

    graphics::text(labels = as.character(fres$number[elt.wdw]),
                   x      = indxy,
                   y      = min(fobs),
                   col    = as.character(fres$color[elt.wdw]),
                   pos    = 4  )

    graphics::text(labels = as.character(fres$group[elt.wdw]),
                   x      = indxy,
                   y      = max(fobs),
                   col    = as.character(fres$color[elt.wdw]),
                   font   = 3, pos = 4)
  }

  graphics::title(main)
}




#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Main Plot Function                                                       ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#'
#' @title Check options
#'
#' @description Check the consistence of options for plotting fclust-results
#'
#' @param opt a list of options for a considered plot function
#'
#' @param all a list that contains all possible options
#'   for the  considered plot function
#'
#' @param main the name of plot function
#'
#' @param windowmin a vector of indices
#'
#' @param affectElt a vector of integer.
#'
#' @details None.
#'
#' @return a list of consistant options
#'
#' @importFrom grDevices colours
#'
#' @keywords internal
#'
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

check_plot_options <- function(opt, all,
                               main = "", windowmin = NULL, affectElt = NULL) {


  # delete possible option = "FALSE"

  bool <- logical(length(opt))
  for (i in seq_along(opt)) if (opt[i] != "FALSE") bool[i] <- TRUE
  index <- which(bool == FALSE)
  if (is.list(opt) & (length(index) > 0))     opt[[index]] <- NULL
  if (is.character(opt) & (length(bool) > 0)) opt <- opt[bool]


  # check the general syntax of options

  new.opt <- vector(mode = "list", length = length(opt))

  for (i in seq_along(new.opt)) {

    if ( is.null(names(opt)[i]) || (names(opt)[i] == "") ) {

      tmp               <- check_foption(opt[[i]], names(all))
      new.opt[i]        <- tmp
      names(new.opt)[i] <- ifelse(is.null(names(tmp)), tmp, names(tmp))

    } else {

      # check if names("opt") is valid
      tmp <- check_foption(substr(names(opt[i]), 1, 3), names(all))

      # check if some names("opt") == TRUE: valid the TRUE-options
      new.opt[[i]] <- ifelse(opt[[i]] == "TRUE", names(opt[i]), opt[[i]])

      names(new.opt)[i] <- tmp
    }
  }



  # get together different elements corresponding to a same option
  #   (mainly option for "cols")

  index <- which(table(names(new.opt)) != 1)

  if (length(index) > 0) {
    nitem <- which(names(new.opt) == names(index))
    tmp   <- c(unlist(new.opt[nitem]))
    names(tmp) <- NULL

    new.opt[[nitem[1]]] <- tmp
    new.opt[nitem[2:length(nitem)]] <- NULL
  }


  # check if "opt" = "all"
  if ("all" %in% opt) {
    if ( ("pvalue" %in% names(new.opt)) && (length(opt$pvalue) == 1) )
                            all$pvalue <- opt$pvalue
    new.opt <- all
  }

  opt <- new.opt



  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # check if "hor" and "ver" are given TRUE

  if ( ("hor" %in% opt) & ("ver" %in% opt) ) opt$ver <- NULL
#    stop("you must choice between 'hor' or 'ver' option")



  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  # check if "pvalue" is given with "aov"
  if ("aov" %in% opt) {
    if (length(opt$pvalue) != 1) {
      opt$pvalue <- all$pvalue
    } else {
      if (!is.numeric(opt$pvalue)) {
        stop("pvalue is numeric")
      } else {
        if ( (opt$pvalue <= 0) | (opt$pvalue >= 1) ) {
          stop("1 > pvalue > 0")
        }
      }
    }
  }


  # check if "pvalue" is a numeric value
  if (length(opt$pvalue) == 0) opt$pvalue <- all$pvalue


  # check if "window" is given with "zoom"
  if ("zoom" %in% opt) {
    if (length(opt$window) != 1) {
      stop("option 'zoom' needs an integer 'window'")
    } else {
      if (!is.numeric(opt$window)) {
        stop("'window' is integer")
      } else {
        if ( (opt$window <= 0) ) {
          stop("window > 0")
        } else {
          if (opt$window < windowmin) opt$window <- windowmin
        }
      }
    }
  }


  # check "cols" option
  if (length(opt) > 0) {
    if (length(opt$cols) <= 1) {
      opt$cols <- extend_vector(
        fcolours(), length(affectElt))[shift_affectElt(affectElt)]
    } else {

      if (is.character(opt$cols) == TRUE) {
        # characters correspond only to colours
        if (sum(opt$cols %in% grDevices::colours()) == length(opt$cols)) {
          opt$cols <- extend_vector(
            opt$cols, length(affectElt))[shift_affectElt(affectElt)]
        } else {
          # characters are user's names
          opt$cols <- char_to_int(opt$cols)
          opt$cols <- extend_vector(
            fcolours(), length(affectElt))[shift_affectElt(opt$cols)]
        }
      }
    }
  }

  return(opt)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Plot trees resulting from functional clustering
#'
#' @description The function plots
#' primary and secondary hierarchical trees of component clustering.
#'
#' @usage
#' plot_ftrees(fres, nbcl = 0, main = "Title", opt.tree = NULL )
#'
#' @inheritParams fclust_plot
#'
#' @details None.
#'
#' @return Nothing. It is a procedure.
#'
#' @seealso
#' \code{\link{plot_ftrees}} plot primary and secondary trees
#'  resulting from functional clustering \cr
#' \code{\link{plot_fperf}} plot observed, modelled and predicted performances
#'  resulting from functional clustering \cr
#' \code{\link{plot_fass}} plot performances of some given assemblages \cr
#' \code{\link{plot_fmotif}} plot as boxplot mean performances
#' of assemblages sorted by assembly motifs \cr
#' \code{\link{plot_fcomp}} plot as boxplot mean performances
#' of assemblages containing a given component \cr
#' \code{\link{fclust_plot}} plot all possible outputs
#' of functional clustering
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot_ftrees <- function(fres,
                        nbcl     = 0,
                        main     = "Title",
                        opt.tree = NULL  )  {

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  # Check inputs
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  nbElt <- dim(fres$mOccur)[2]

  nbOpt <- fres$nbOpt
  nbMax <- first_argmax(fres$tStats[ ,"R2prd"])
  if (nbcl != 0)    nbOpt <- nbcl
  if (nbcl > nbMax) nbOpt <- nbMax

  affectElt <- cut_ftree(fres$tree.II, nbOpt)
  cutline   <- amean(fres$tree.II$cor[(nbOpt - 1):nbOpt])
  ind.tree  <- sort_ftree(fres$tree.II, index.return = TRUE)$ix

  windowmin <- which(affectElt[ind.tree] == 1)[1]


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  # Check options
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  all.tree <- list("cal", "prd", "leg", cols = fstd_colour(),
                   "zoom", window = fwindow(),
                   "all")
  names(all.tree) <- list("cal", "prd", "leg", "cols",
                          "zoom", "window",
                          "all")
  opt.tree <- check_plot_options(opt.tree, all.tree, main = "opt.tree",
                                 windowmin, affectElt)


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  # Plot the primary, calibration, whole Tree
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if ("cal" %in% opt.tree) {

    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    # plot whole primary tree
    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    plot_ftree(fres$tree.I, cols = opt.tree$cols,
              main = paste(main, "opt.tree", "cal", sep = " / "))
    graphics::lines(x = c(0, nbElt + 1),
                    y = rep(cutline, 2),
                    col = "blue", lty = "dashed")
    graphics::text(x = 1, y = cutline + 0.02,
                   labels = paste0("nbOpt = ", nbOpt), col = "blue")
    graphics::lines(x = c(0, nbElt + 1),
                    y = rep(fres$tStats[nbOpt, "R2prd"], 2),
                    col = "red", lty = "dashed")
    graphics::text(x = 1, y = fres$tStats[nbOpt, "R2prd"] + 0.02,
                   labels = "E-value", col = "red")
    graphics::text(
      x = floor(nbElt/2), y = 0,
      labels = paste0("R2 = ",
                      signif(fres$tStats[nbOpt, "R2cal"], digits = 3),
                      "    E = ",
                      signif(fres$tStats[nbOpt, "R2prd"], digits = 3),
                      "    E/R2 = ",
                      signif(fres$tStats[nbOpt, "R2prd"] /
                               fres$tStats[nbOpt, "R2cal"], digits = 3)),
      col = "red")


    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    # zoom on the first part of primary tree
    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    if ("zoom" %in% opt.tree) if (opt.tree$window < nbElt) {

      wdw <- logical(nbElt)
      wdw[ind.tree][1:opt.tree$window] <- TRUE

      plot_ftree(simplify_ftree(fres$tree.I, wdw),
                 cols  = opt.tree$col[wdw],
                main = paste(main, "opt.tree", "cal",
                              paste0("zoom=", opt.tree$window), sep = " / ") )
      graphics::lines(x = c(0, length(wdw) + 1),
                      y = rep(cutline, 2),
                      col = "blue", lty = "dashed")
      graphics::text(x = 1, y = cutline + 0.02,
                     labels = paste0("nbOpt = ", nbOpt), col = "blue")
      graphics::lines(x = c(0, length(wdw) + 1),
                      y = rep(fres$tStats[nbOpt, "R2prd"], 2),
                      col = "red")
      graphics::text(x = 1, y = fres$tStats[nbOpt, "R2prd"] + 0.02,
                     labels = "E-value", col = "red")
    }
  }



  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  # Plot the secondary, validated Tree
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if ("prd" %in% opt.tree) {

    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    # plot full secondary tree
    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    plot_ftree(fres$tree.II, cols = opt.tree$col,
              main = paste(main, "opt.tree", "prd", sep = " / "))
    graphics::lines(x = c(0, nbElt + 1),
                    y = rep(fres$tStats[nbOpt, "R2prd"], 2),
                    col = "red", lty = "dashed")
    graphics::text(x = 1, y = fres$tStats[nbOpt, "R2prd"] + 0.02,
                   labels = "E-value", col = "red")
    graphics::text(
      x = floor(nbElt/2), y = 0,
      labels = paste0("R2 = ",
                      signif(fres$tStats[nbOpt, "R2cal"], digits = 3),
                      "    E = ",
                      signif(fres$tStats[nbOpt, "R2prd"], digits = 3),
                      "    E/R2 = ",
                      signif(fres$tStats[nbOpt, "R2prd"] /
                               fres$tStats[nbOpt, "R2cal"], digits = 3)),
      col = "red")



    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    # zoom on the first part of secondary, validated tree
    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    if ("zoom" %in% opt.tree) if (opt.tree$window < nbElt) {

      wdw <- logical(nbElt)
      wdw[ind.tree][1:opt.tree$window] <- TRUE

      plot_ftree(simplify_ftree(fres$tree.II, wdw),
                cols  = opt.tree$col[wdw],
                main = paste(main, "opt.tree", "prd",
                              paste0("zoom=", opt.tree$window), sep = " / ") )
      graphics::lines(x = c(0, length(wdw) + 1),
                      y = rep(fres$tStats[nbOpt, "R2prd"], 2),
                      col = "red")
      graphics::text(x = 1, y = fres$tStats[nbOpt, "R2prd"] + 0.02,
                     labels = "E-value", col = "red")
    }
  }



  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  # Plot the cluster content
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if ("leg" %in% opt.tree)
    plot_clusters_content(affectElt, opt.sort = TRUE)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Plot modelled and predicted performances
#'  resulting from functional clustering
#'
#' @description The function plots
#' observed, modelled and predicted performances
#'  resulting from functional clustering
#'
#' @usage
#' plot_fperf(fres, nbcl = 0, main = "Title", opt.perf = NULL )
#'
#' @inheritParams fclust_plot
#'
#' @details None.
#'
#' @return Nothing. It is a procedure.
#'
#' @seealso
#' \code{\link{plot_ftrees}} plot primary and secondary trees
#'  resulting from functional clustering \cr
#' \code{\link{plot_fperf}} plot observed, modelled and predicted performances
#'  resulting from functional clustering \cr
#' \code{\link{plot_fass}} plot performances of some given assemblages \cr
#' \code{\link{plot_fmotif}} plot as boxplot mean performances
#' of assemblages sorted by assembly motifs \cr
#' \code{\link{plot_fcomp}} plot as boxplot mean performances
#' of assemblages containing a given component \cr
#' \code{\link{fclust_plot}} plot all possible outputs
#' of functional clustering
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot_fperf <- function(fres,
                       nbcl      = 0,
                       main     = "Title",
                       opt.perf  = NULL )  {

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  # Check inputs
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  nbAss    <- dim(fres$mOccur)[1]

  figures  <- extend_vector(fsymbols(), nbAss)
  couleurs <- extend_vector(fcolours(), nbAss)

  nbOpt <- fres$nbOpt
  nbMax <- first_argmax(fres$tStats[ ,"R2prd"])
  if (nbcl != 0)    nbOpt <- nbcl
  if (nbcl > nbMax) nbOpt <- nbMax

  affectElt <- cut_ftree(fres$tree.II, nbOpt)
  ind.tree  <- sort_ftree(fres$tree.II, index.return = TRUE)$ix

  setXpr <- unique(names(fres$xpr))
  xpr.ON <- ifelse(length(setXpr) > 1, TRUE, FALSE)

  windowmin <- which(affectElt[ind.tree] == 1)[1]



  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  # Check options
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  all.perf <- list("stats_I", "stats_II",
                   "cal", "prd", "missing", "pub", "calprd",
                   "seq", "ass", "xpr", "multi", "aov", pvalue = fpvalue(),
                   "all")
  names(all.perf) <- list("stats_I", "stats_II",
                          "cal", "prd", "missing", "pub", "calprd",
                          "seq", "ass", "xpr", "multi", "aov", "pvalue",
                          "all")
  opt.perf <- check_plot_options(opt.perf, all.perf, main = "opt.perf",
                                 windowmin, affectElt)


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  # Plot Statistitics of calibration goodness-of-fit
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if ("stats_I" %in% opt.perf) {

    plot_stats(fres$mStats,
               main = paste(main, "perf", "stats_I",
                             fres$opt.mean, fres$opt.model, sep = " / ") )
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  # Plot Calibrations for all the number of clusters of components
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if ("cal" %in% opt.perf) {

    nbMin <- ifelse("seq" %in% opt.perf, 1, nbOpt)

    for (nbcl in nbMin:nbOpt) {

      plot_prediction_simple(Fx       = fres$fobs,
                             Fy       = fres$mCal[nbcl, ],
                             assMotif = fres$mMotifs[nbcl, ],
                             xylim    = range(fres$fobs),
                             figs     = figures,
                             cols     = couleurs,
                             nbcl     = nbcl,
                             opt.mean = fres$opt.mean,
                             pvalue   = opt.perf$pvalue,
                             main     = paste(main, "opt.perf", "cal",
                                               fres$opt.mean, fres$opt.model,
                                               paste0("nbcl=", nbcl),
                                               sep = " / ") )

      if ("ass" %in% opt.perf)
        add_ass_names(Fx        = fres$fobs,
                      Fy        = fres$tCal[nbcl, ],
                      assMotif = fres$mMotifs[nbcl, ],
                      cols      = couleurs  )
    }
  }



  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  # Plot Statistics of Tree predictions
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if ("stats_II" %in% opt.perf) {

    plot_stats(fres$tStats,
               main = paste(main, "opt.perf", "stats_II",
                             fres$opt.mean, fres$opt.model,
                             sep = " / ") )
  }




  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  # Plot Predictions for all the number of clusters of components
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if ("prd" %in% opt.perf) {

    nbMin <- ifelse("seq" %in% opt.perf, 1, nbOpt)

    for (nbcl in nbMin:nbOpt)
      if (sum(!is.na(fres$mPrd[nbcl, ])) > 0) {

        plot_prediction_LOO(Fx        = fres$fobs,
                            Fy        = fres$tCal[nbcl, ],
                            Fz        = fres$tPrd[nbcl, ],
                            assMotif  = fres$mMotifs[nbcl, ],
                            xylim     = range(fres$fobs),
                            figs      = figures,
                            cols      = couleurs,
                            nbcl      = nbcl,
                            opt.mean  = fres$opt.mean,
                            pvalue    = opt.perf$pvalue,
                            main      = paste(main, "opt.perf", "prd",
                                              fres$opt.mean, fres$opt.model,
                                              paste0("nbcl=", nbcl),
                                              sep = " / ") )

        if ("ass" %in% opt.perf)
          add_ass_names(Fx        = fres$fobs,
                        Fy        = fres$tCal[nbcl, ],
                        assMotif = fres$mMotifs[nbcl, ],
                        cols      = couleurs  )
      }
  }



  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  # Plot Tree Predictions for all numbers of clusters of components.
  #    Different colours indicate the number of clusters needs
  #     to predict the performance of assemblages.
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if ("missing" %in% opt.perf)  {

    nbMin <- ifelse("seq" %in% opt.perf, 1, nbOpt)

    for (nbcl in nbMin:nbOpt) {

      plot_prediction_LOO(Fx       = fres$fobs,
                          Fy       = fres$tCal[nbcl, ],
                          Fz       = fres$tPrd[nbcl, ],
                          assMotif = fres$tNbcl[nbcl, ],
                          xylim    = range(fres$fobs),
                          figs     = figures,
                          cols     = couleurs,
                          nbcl     = nbcl,
                          opt.mean = fres$opt.mean,
                          main     = paste(main, "opt.perf", "missing",
                                           fres$opt.mean, fres$opt.model,
                                           paste0("nbcl=", nbcl),
                                           sep = " / ") )

      index <- which(fres$tNbcl[nbcl, ] == (nbcl - 1))
      if (length(index) > 0)
        add_ass_names(Fx       = fres$fobs[index],
                      Fy       = fres$tCal[nbcl,  index],
                      assMotif = fres$tNbcl[nbcl, index],
                      cols     = couleurs)
    }
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  # Plot Tree Predictions by cross-validation versus Tree Calibrations
  #           for all numbers of clusters of components
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if ("calprd" %in% opt.perf)  {

    nbMin <- ifelse("seq" %in% opt.perf, 1, nbOpt)

    for (nbcl in nbMin:nbOpt) {

      plot_prediction_simple(Fx       = fres$tCal[nbcl, ],
                             Fy       = fres$tPrd[nbcl, ],
                             assMotif = fres$mMotifs[nbcl, ],
                             xylab    = c("Modelled performance",
                                         "Predicted performance"),
                             xylim    = range(fres$fobs),
                             figs     = figures,
                             cols     = couleurs,
                             nbcl     = nbcl,
                             opt.mean = fres$opt.mean,
                             pvalue   = opt.perf$pvalue,
                             main     = paste(main, "opt.perf", "cal vs prd",
                                              fres$opt.mean, fres$opt.model,
                                              paste0("nbcl=", nbcl),
                                              sep = " / ") )

      if ("ass" %in% opt.perf)
        add_ass_names(Fx       = fres$tCal[nbcl, ],
                      Fy       = fres$tPrd[nbcl, ],
                      assMotif = fres$mMotifs[nbcl, ],
                      cols     = couleurs  )
    }
  }




  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  # Plot Tree Predictions for all numbers of clusters of components
  #    using only a colour*symbol for all assembly motifs
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if ("pub" %in% opt.perf)  {

    nbMin <- ifelse("seq" %in% opt.perf, 1, nbOpt)

    for (nbcl in nbMin:nbOpt) {

      plot_prediction_LOO(Fx        = fres$fobs,
                          Fy        = fres$tCal[nbcl, ],
                          Fz        = fres$tPrd[nbcl, ],
                          assMotif = extend_vector(1, nbAss),
                          xylim     = range(fres$fobs),
                          figs      = extend_vector(figures[1], nbAss),
                          cols      = extend_vector(couleurs[1], nbAss),
                          nbcl      = nbcl,
                          opt.mean  = fres$opt.mean,
                          main      = paste(main, "opt.perf", "pub",
                                            fres$opt.mean, fres$opt.model,
                                            paste0("nbOpt=", nbOpt),
                                            sep = " / ") )
#
#      if ("ass" %in% opt.perf)
#        add_ass_names(Fx       = fres$fobs,
#                      Fy       = fres$tCal[nbcl, ],
#                      assMotif = extend_vector(1, nbAss),
#                      cols     = extend_vector(couleurs[2], nbAss) )
    }
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  # Plot Tree Predictions for the optimum number of clusters of components
  #    experiment by experiment
  #    using a different colour*symbol for each experiment
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (xpr.ON == TRUE) {

    nbcl <- nbOpt
    plot_prediction_LOO(Fx        = fres$fobs,
                        Fy        = fres$tCal[nbcl, ],
                        Fz        = fres$tPrd[nbcl, ],
                        assMotif = names(fres$xpr),
                        xylim     = range(fres$fobs),
                        figs      = figures,
                        cols      = couleurs,
                        nbcl      = nbcl,
                        opt.mean  = fres$opt.mean,
                        pvalue    = opt.perf$pvalue,
                        main      = paste(main, "opt.perf", "prd",
                                          fres$opt.mean, fres$opt.model,
                                          paste0("nbcl=", nbcl),
                                          "xpr=all",
                                          sep = " / ") )


    if (("xpr" %in% opt.perf) | ("multi" %in% opt.perf)) {

      for (ipr in seq_along(setXpr)) {

        index <- which(names(fres$xpr) == setXpr[ipr])
        plot_prediction_LOO(Fx        = fres$fobs[index],
                            Fy        = fres$tCal[nbcl, index],
                            Fz        = fres$tPrd[nbcl, index],
                            assMotif = rep(ipr, length(fres$xpr[index])),
                            xylim     = range(fres$fobs),
                            figs      = figures,
                            cols      = couleurs,
                            nbcl      = nbcl,
                            opt.mean  = fres$opt.mean,
                            pvalue    = opt.perf$pvalue,
                            main      = paste(main, "opt.perf", "prd",
                                              fres$opt.mean, fres$opt.model,
                                              paste0("nbcl=", nbcl),
                                              paste0("xpr=", setXpr[ipr]),
                                              sep = " / ") )

        graphics::abline(v = mean_fct(fres$fobs, fres$opt.mean),
                         lty = "dashed", col = "red")

        if ("ass" %in% opt.perf)
          add_ass_names(Fx        = fres$fobs[index],
                        Fy        = fres$tCal[nbcl, index],
                        assMotif = fres$mMotifs[nbcl, index],
                        cols      = couleurs  )
      }
    }
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Plot performances of some given assemblages
#'
#' @description The function plots performances
#'  of some given, identified assemblages.
#'
#' @usage
#' plot_fass(fres, nbcl = 0, main = "Title", opt.ass = NULL )
#'
#' @inheritParams fclust_plot
#'
#' @details None.
#'
#' @return Nothing. It is a procedure.
#'
#' @seealso
#' \code{\link{plot_ftrees}} plot primary and secondary trees
#'  resulting from functional clustering \cr
#' \code{\link{plot_fperf}} plot observed, modelled and predicted performances
#'  resulting from functional clustering \cr
#' \code{\link{plot_fass}} plot performances of some given assemblages \cr
#' \code{\link{plot_fmotif}} plot as boxplot mean performances
#' of assemblages sorted by assembly motifs \cr
#' \code{\link{plot_fcomp}} plot as boxplot mean performances
#' of assemblages containing a given component \cr
#' \code{\link{fclust_plot}} plot all possible outputs
#' of functional clustering
#'
#' @importFrom graphics abline

#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot_fass <- function(fres,
                      nbcl    = 0,
                      main   = "Title",
                      opt.ass = NULL  )  {


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  # Check inputs
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  nbAss    <- dim(fres$mOccur)[1]

  figures  <- extend_vector(fsymbols(), nbAss)
  couleurs <- extend_vector(fcolours(), nbAss)

  nbOpt <- fres$nbOpt
  nbMax <- first_argmax(fres$tStats[ ,"R2prd"])
  if (nbcl != 0)    nbOpt <- nbcl
  if (nbcl > nbMax) nbOpt <- nbMax

  affectElt <- cut_ftree(fres$tree.II, nbOpt)
  ind.tree  <- sort_ftree(fres$tree.II, index.return = TRUE)$ix

  setXpr <- unique(names(fres$xpr))
  xpr.ON  <- ifelse(length(setXpr) > 1, TRUE, FALSE)

  windowmin <- which(affectElt[ind.tree] == 1)[1]



  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  #  Check options
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  all.ass <- list(sample = 5, who = NULL)
  names(all.ass) <- list("sample", "who")
  opt.ass <- check_plot_options(opt.ass, all.ass, main = "opt.ass",
                                windowmin, affectElt)

  if (length(opt.ass$sample) > 0)
    opt.ass$who <- c(opt.ass$who,
                     sample(names(fres$fobs), size = opt.ass$sample))



  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  # if there are several experiments,
  #  plot Tree Predictions by experiment (setXpr) with the names of Assemblages
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if ( xpr.ON & (length(opt.ass$who) > 0) ) {

    for (ass in opt.ass$who) {

      index <- which(names(fres$fobs) == ass)
      plot_prediction_LOO(Fx        = fres$fobs[index],
                          Fy        = fres$tCal[nbOpt, index],
                          Fz        = fres$tPrd[nbOpt, index],
                          assMotif = fres$mMotifs[nbOpt, index],
                          xylim     = range(fres$fobs),
                          figs      = figures,
                          cols      = couleurs,
                          nbcl      = nbOpt,
                          opt.mean  = fres$opt.mean,
                          main      = paste(main,
                                            paste0("assemblage = '", ass, "'"),
                                            sep = " / ") )

      graphics::abline(v = mean_fct(fres$fobs, fres$opt.mean),
                       h = mean_fct(fres$tCal[nbOpt, ], fres$opt.mean),
                       lty = "dashed", col = "red")

      Fx        <- fres$fobs
      names(Fx) <- names(fres$xpr)
      add_ass_names(Fx        = Fx[index],
                    Fy        = fres$tCal[nbOpt, index],
                    assMotif = fres$mMotifs[nbOpt, index],
                    cols      = couleurs  )
    }
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Plot as boxplot mean performances
#' of assemblages sorted by assembly motifs
#'
#' @description The function plots, as vertical or horizontal boxplots,
#' composition and mean performance sorted by assembly motifs.
#'
#' @usage
#' plot_fmotif(fres, nbcl = 0, main = "Title", opt.motif = NULL )
#'
#' @inheritParams fclust_plot
#'
#' @details None.
#'
#' @return Nothing. It is a procedure.
#'
#' @seealso
#' \code{\link{plot_ftrees}} plot primary and secondary trees
#'  resulting from functional clustering \cr
#' \code{\link{plot_fperf}} plot observed, modelled and predicted performances
#'  resulting from functional clustering \cr
#' \code{\link{plot_fass}} plot performances of some given assemblages \cr
#' \code{\link{plot_fmotif}} plot as boxplot mean performances
#' of assemblages sorted by assembly motifs \cr
#' \code{\link{plot_fcomp}} plot as boxplot mean performances
#' of assemblages containing a given component \cr
#' \code{\link{fclust_plot}} plot all possible outputs
#' of functional clustering
#'
#' @importFrom graphics abline
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot_fmotif <- function(fres,
                        nbcl      = 0,
                        main      = "Title",
                        opt.motif = NULL    )  {

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  # Check inputs
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  nbOpt <- fres$nbOpt
  nbMax <- first_argmax(fres$tStats[ ,"R2prd"])
  if (nbcl != 0)    nbOpt <- nbcl
  if (nbcl > nbMax) nbOpt <- nbMax

  affectElt <- cut_ftree(fres$tree.II, nbOpt)
  ind.tree  <- sort_ftree(fres$tree.II, index.return = TRUE)$ix

  setXpr <- unique(names(fres$xpr))
  xpr.ON <- ifelse(length(setXpr) > 1, TRUE, FALSE)

  windowmin <- which(affectElt[ind.tree] == 1)[1]



  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  #   Check options
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  all.motif <- list("obs", "cal", "prd", "leg",
                    "hor", "ver", "seq", "xpr", "multi",
                    "aov", pvalue = fpvalue(),
                    "all")
  names(all.motif) <- list("obs", "cal", "prd", "leg",
                           "hor", "ver", "seq", "xpr", "multi",
                           "aov", "pvalue",
                           "all")
  opt.motif <- check_plot_options(opt.motif, all.motif,
                            main = "opt.motif", windowmin, affectElt)


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  # In any case, motifs are sorted by decreasing mean observed performance
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  # compute orders and statistics

  opt.fct <- NULL
  if ("obs" %in% opt.motif) opt.fct <- c(opt.fct, "obs")
  if ("cal" %in% opt.motif) opt.fct <- c(opt.fct, "cal")
  if ("prd" %in% opt.motif) opt.fct <- c(opt.fct, "prd")

  decreasing.ON <- ifelse(("hor" %in% opt.motif), FALSE, TRUE)
  horizontal.ON <- ifelse(("hor" %in% opt.motif), TRUE,  FALSE)

  ref.fobs <- list()
  for (nbcl in 1:nbOpt) {

    assMot <- fres$mMotifs[nbcl, ]
    assNam <- name_motifs(cut_ftree(fres$tree.II, nbcl), fres$mOccur)

    ref.fobs[[nbcl]] <- sort_motifs(fres$fobs, assMot, assNam,
                                    pvalue  = opt.motif$pvalue,
                                    opt.dec = decreasing.ON  )
  }


    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    #
    # Assembly motifs are plotted
    #            versus Observed, Simulated or Predicted Performances
    #
    # Loop on "obs", "cal" or "prd"
    #
    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (length(opt.fct) > 0) {

    for (indFct in opt.fct) {

      main3 <- switch(indFct,
                      obs = "Observed performance",
                      cal = "Modelled performance",
                      prd = "Predicted performance" )

      nbMin <- ifelse(("seq" %in% opt.motif), 1, nbOpt)

      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      #
      # Loop on the number of clusters of components,
      #                from 2 to nbOpt-1 clusters of components
      #
      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      for (nbcl in nbMin:nbOpt) {

        fct   <- switch(indFct,
                        obs = fres$fobs,
                        cal = fres$tCal[nbcl, ],
                        prd = fres$tPrd[nbcl, ])

        assMot <- fres$mMotifs[nbcl, ]
        assNam <- name_motifs(cut_ftree(fres$tree.II, nbcl), fres$mOccur)

        #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
        # plot all components over all experiments
        #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

        ref.nbcl <- sort_motifs(fct, assMot, assNam,
                                pvalue  = opt.motif$pvalue,
                                opt.dec = decreasing.ON  )

        main2 <- paste(main, "opt.motif", indFct,
                       paste0("nbcl=", nbcl), sep = " / ")
        if (xpr.ON == TRUE) main2 <- paste(main2, "xpr=all", sep = " / ")

        plot_motifs_box(fct, assMot,
                        sres     = ref.nbcl,
                        ordre    = as.integer(ref.fobs[[nbcl]]$motif),
                        xylim    = range(fres$fobs),
                        main     = main2,
                        xlabel   = main3,
                        pvalue   = opt.motif$pvalue,
                        opt.mean = fres$opt.mean,
                        opt.hor  = horizontal.ON  )

        if (horizontal.ON == TRUE) {
          graphics::abline(v = mean_fct(fres$fobs, fres$opt.mean),
                           lty = "dashed", col = "red")
        } else {
          graphics::abline(h = mean_fct(fres$fobs, fres$opt.mean),
                           lty = "dashed", col = "red")
        }

        if ("leg" %in% opt.motif) {
          index <- which(names(fres$xpr) == setXpr[1])
          plot_motifs_content(fct[index], assMot[index], assNam[index])
        }
      }

      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      # plot motifs for each experiment
      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      if ((xpr.ON == TRUE) &
           (("xpr" %in% opt.motif) | ("multi" %in% opt.motif)) ) {

        for (ipr in seq_along(setXpr)) {

          index <- which(names(fres$xpr) == setXpr[ipr])
          ref.xpr <- sort_motifs(fct[index],assMot[index], assNam[index],
                                 pvalue  = opt.motif$pvalue,
                                 opt.dec = decreasing.ON  )

          plot_motifs_box(fct[index], assMot[index],
                          sres     = ref.xpr,
                          ordre    = as.integer(ref.fobs[[nbcl]]$motif),
                          xylim    = range(fres$fobs),
                          main     = paste(main, "opt.motif", indFct,
                                           paste0("nbcl=", nbcl),
                                           paste0("xpr=", setXpr[ipr]),
                                           sep = " / "),
                          xlabel   = main3,
                          pvalue   = opt.motif$pvalue,
                          opt.mean = fres$opt.mean,
                          opt.hor  = horizontal.ON  )

          if (horizontal.ON == TRUE) {
            graphics::abline(v = mean_fct(fres$fobs, fres$opt.mean),
                             lty = "dashed", col = "red")
          } else {
            graphics::abline(h = mean_fct(fres$fobs, fres$opt.mean),
                             lty = "dashed", col = "red")
          }
        }
      }
    }
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Plot mean performances
#' of assemblages containing a given component
#'
#' @description The function plots, as vertical or horizontal boxplots,
#'  the mean performances
#'  of assemblages containing a given component.
#'
#' @usage
#' plot_fcomp(fres, nbcl = 0, main = "Title", opt.comp = NULL )
#'
#' @inheritParams fclust_plot
#'
#' @details None.
#'
#' @return Nothing. It is a procedure.
#'
#' @seealso
#' \code{\link{plot_ftrees}} plot primary and secondary trees
#'  resulting from functional clustering \cr
#' \code{\link{plot_fperf}} plot observed, modelled and predicted performances
#'  resulting from functional clustering \cr
#' \code{\link{plot_fass}} plot performances of some given assemblages \cr
#' \code{\link{plot_fmotif}} plot as boxplot mean performances
#' of assemblages sorted by assembly motifs \cr
#' \code{\link{plot_fcomp}} plot as boxplot mean performances
#' of assemblages containing a given component \cr
#' \code{\link{fclust_plot}} plot all possible outputs
#' of functional clustering
#'
#' @importFrom graphics abline
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot_fcomp <- function(fres,
                       nbcl      = 0,
                       main     = "Title",
                       opt.comp  = NULL  )  {

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  # Check inputs
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  nbElt <- dim(fres$mOccur)[2]

  nbOpt <- fres$nbOpt
  nbMax <- first_argmax(fres$tStats[ ,"R2prd"])
  if (nbcl != 0)    nbOpt <- nbcl
  if (nbcl > nbMax) nbOpt <- nbMax

  affectElt <- cut_ftree(fres$tree.II, nbOpt)
  ind.tree  <- sort_ftree(fres$tree.II, index.return = TRUE)$ix

  setXpr <- unique(names(fres$xpr))
  xpr.ON <- ifelse(length(setXpr) > 1, TRUE, FALSE)

  windowmin <- which(affectElt[ind.tree] == 1)[1]


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  #  Check options
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  # Plot the performances of assemblages containing a given component
  #                                     (where a given component occurs).
  # Centroids of assmblages are sorted by increasing effect on performance.
  # For comparing with assembly motifs
  #                  horizontal.ON = TRUE
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  all.component <- list("tree", "perf", "hor", "ver", "leg",
                        "xpr", "multi",
                        cols = fstd_colour(), "aov", pvalue = fpvalue(),
                        "zoom", window = fwindow(),
                        "all")
  names(all.component) <- list("tree", "perf", "hor", "ver", "leg",
                               "xpr", "multi",
                               "cols", "aov", "pvalue",
                               "zoom", "window",
                               "all")
  opt.comp <- check_plot_options(opt.comp, all.component,
                                 main = "opt.comp", windowmin, affectElt)


  # compute orders and statistics

  if ( ("tree" %in% opt.comp) | ("perf" %in% opt.comp) ) {

    decreasing.ON <- ifelse(("hor" %in% opt.comp), FALSE, TRUE)
    horizontal.ON <- ifelse(("hor" %in% opt.comp), TRUE,  FALSE)

    ref.all <- sort_components(fres$fobs, fres$mOccur,
                               cols    = opt.comp$cols,
                               pvalue  = opt.comp$pvalue,
                               opt.dec = decreasing.ON)

    if (xpr.ON == TRUE) {

      ref.xpr <- list()
      for (ipr in seq_along(setXpr)) {

        index <- which(names(fres$xpr) == setXpr[ipr])
        ref.xpr[[ipr]] <- sort_components(fres$fobs[index],
                                          fres$mOccur[index, , drop = FALSE],
                                          cols    = opt.comp$cols,
                                          pvalue  = opt.comp$pvalue,
                                          opt.dec = decreasing.ON )
      }
    }
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  # components are sorted as the clustering tree
  #                 and plotted in horizontal boxplots
  #
  # horizontal.ON <- TRUE;  decreasing.ON <- FALSE
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if ( ("tree" %in% opt.comp) & (horizontal.ON == TRUE) ) {

    indOrd  <- ref.all$ind.inturn[rev(ind.tree)]
    ordre   <- as.integer(ref.all$motif[indOrd])

    if (xpr.ON == TRUE) {

      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      # plot all components over all experiments
      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      plot_components_box(fres$fobs, fres$mOccur,
                          sres     = ref.all[indOrd, ],
                          ordre    = ordre,
                          main     = paste(main, "opt.comp", "tree",
                                           "xpr=all", sep = " / "),
                          pvalue   = opt.comp$pvalue,
                          opt.mean = fres$opt.mean,
                          opt.hor  = horizontal.ON )

      if ("leg" %in% opt.comp)
        plot_components_content(
          fres$mOccur[which(names(fres$xpr) == setXpr[1]), , drop = FALSE],
          ordre = ordre)

      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      # plot all components for each experiment
      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      if ( ("xpr" %in% opt.comp) | ("multi" %in% opt.comp) ) {

        for (ipr in seq_along(setXpr)) {

          index <- which(names(fres$xpr) == setXpr[ipr])
          plot_components_box(fres$fobs[index],
                              fres$mOccur[index, , drop = FALSE],
                              sres     = ref.xpr[[ipr]][indOrd, ],
                              ordre    = ordre,
                              main     = paste(main, "opt.comp", "tree",
                                               paste0("xpr=", setXpr[ipr]),
                                               sep = " / "),
                              pvalue   = opt.comp$pvalue,
                              opt.mean = fres$opt.mean,
                              opt.hor  = horizontal.ON )

          if (horizontal.ON == TRUE) {
            graphics::abline(v = mean_fct(fres$fobs, fres$opt.mean),
                             lty = "dashed", col = "red")
          } else {
            graphics::abline(h = mean_fct(fres$fobs, fres$opt.mean),
                             lty = "dashed", col = "red")
          }
        }
      }

      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      # zoom on the first part over all experiments
      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      if ("zoom" %in% opt.comp) if (opt.comp$window < nbElt) {

        plot_components_box(fres$fobs, fres$mOccur,
                            sres     = ref.all[indOrd, ],
                            ordre    = ordre,
                            elt.wdw  = rev(nbElt - 1:opt.comp$window + 1),
                            main     = paste(main, "opt.comp", "tree",
                                             "xpr=all",
                                             paste0("zoom=1:", opt.comp$window),
                                             sep = " / "),
                            pvalue   = opt.comp$pvalue,
                            opt.mean = fres$opt.mean,
                            opt.hor  = horizontal.ON )

        #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
        # zoom on the first part of each experiment
        #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

        if ( ("xpr" %in% opt.comp) | ("multi" %in% opt.comp) ) {

          for (ipr in seq_along(setXpr)) {

            index <- which(names(fres$xpr) == setXpr[ipr])
            plot_components_box(fres$fobs[index],
                                fres$mOccur[index, , drop = FALSE],
                                sres    = ref.xpr[[ipr]][indOrd, ],
                                ordre   = ordre,
                                elt.wdw = rev(nbElt - 1:opt.comp$window + 1),
                                main    = paste(main, "opt.comp", "tree",
                                                paste0("xpr=", setXpr[ipr]),
                                                paste0("zoom=1:",
                                                       opt.comp$window),
                                                sep = " / "),
                                pvalue   = opt.comp$pvalue,
                                opt.mean = fres$opt.mean,
                                opt.hor  = horizontal.ON )

            if (horizontal.ON == TRUE) {
              graphics::abline(v = mean_fct(fres$fobs, fres$opt.mean),
                               lty = "dashed", col = "red")
            } else {
              graphics::abline(h = mean_fct(fres$fobs, fres$opt.mean),
                               lty = "dashed", col = "red")
            }
          }
        }
      }

    } else {

      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      # plot all components
      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      plot_components_box(fres$fobs, fres$mOccur,
                          sres     = ref.all[indOrd, ],
                          ordre    = ordre,
                          main     = paste(main, "opt.comp", "tree",
                                           sep = " / "),
                          pvalue   = opt.comp$pvalue,
                          opt.mean = fres$opt.mean,
                          opt.hor  = horizontal.ON )

      if ("leg" %in% opt.comp)
        plot_components_content(fres$mOccur, ordre = ordre)

      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      # zoom on the first part of primary tree
      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      if ("zoom" %in% opt.comp) if (opt.comp$window < nbElt) {

        plot_components_box(fres$fobs, fres$mOccur,
                            sres    = ref.all[indOrd, ],
                            ordre   = ordre,
                            elt.wdw = rev(nbElt - 1:opt.comp$window + 1),
                            main    = paste(main, "opt.comp", "tree",
                                            paste0("zoom=1:", opt.comp$window),
                                            sep = " / "),
                            pvalue   = opt.comp$pvalue,
                            opt.mean = fres$opt.mean,
                            opt.hor  = horizontal.ON )
      }
    }
  }



  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  # components are sorted as the clustering tree
  #                 and plotted in vertical boxplots
  #
  # horizontal.ON <- FALSE;  decreasing.ON <- TRUE
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if ( ("tree" %in% opt.comp) & (horizontal.ON == FALSE) ) {

    indOrd  <- ref.all$ind.inturn[ind.tree]
    ordre   <- as.integer(ref.all$motif[indOrd])

    if (xpr.ON == TRUE) {

      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      # plot all components over all experiments
      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      plot_components_box(fres$fobs, fres$mOccur,
                          sres    = ref.all[indOrd, ],
                          ordre   = ordre,
                          main    = paste(main, "opt.comp", "tree",
                                          "xpr=all", sep = " / "),
                          pvalue  = opt.comp$pvalue,
                          opt.mean = fres$opt.mean,
                          opt.hor = horizontal.ON )

      if ("leg" %in% opt.comp)
        plot_components_content(
          fres$mOccur[which(names(fres$xpr) == setXpr[1]), , drop = FALSE],
          ordre = rev(ordre))

      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      # plot all components for each experiment
      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      if ( ("xpr" %in% opt.comp) | ("multi" %in% opt.comp) ) {

        for (ipr in seq_along(setXpr)) {

          index      <- which(names(fres$xpr) == setXpr[ipr])
          plot_components_box(fres$fobs[index],
                              fres$mOccur[index, , drop = FALSE],
                              sres    = ref.xpr[[ipr]][indOrd, ],
                              ordre   = ordre,
                              main    = paste(main, "opt.comp", "tree",
                                              paste0("xpr=", setXpr[ipr]),
                                              sep = " / "),
                              pvalue  = opt.comp$pvalue,
                              opt.mean = fres$opt.mean,
                              opt.hor = horizontal.ON )

          if (horizontal.ON == TRUE) {
            graphics::abline(v = mean_fct(fres$fobs, fres$opt.mean),
                             lty = "dashed", col = "red")
          } else {
            graphics::abline(h = mean_fct(fres$fobs, fres$opt.mean),
                             lty = "dashed", col = "red")
          }
        }
      }


      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      # zoom on the first part over all experiments
      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      if ("zoom" %in% opt.comp) if (opt.comp$window < nbElt) {

        plot_components_box(fres$fobs, fres$mOccur,
                            sres    = ref.all[indOrd, ],
                            ordre   = ordre,
                            elt.wdw = 1:opt.comp$window,
                            main    = paste(main, "opt.comp", "tree",
                                           "xpr=all",
                                           paste0("zoom=1:", opt.comp$window),
                                           sep = " / "),
                            pvalue  = opt.comp$pvalue,
                            opt.mean = fres$opt.mean,
                            opt.hor = horizontal.ON )

        #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
        # zoom on the first part of each experiment
        #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

        if ( ("xpr" %in% opt.comp) | ("multi" %in% opt.comp) ) {

          for (ipr in seq_along(setXpr)) {

            index   <- which(names(fres$xpr) == setXpr[ipr])
            plot_components_box(fres$fobs[index],
                                fres$mOccur[index, , drop = FALSE],
                                sres    = ref.xpr[[ipr]][indOrd, ],
                                ordre   = ordre,
                                elt.wdw = 1:opt.comp$window,
                                main    = paste(main, "opt.comp", "tree",
                                               paste0("xpr=", setXpr[ipr]),
                                               paste0("zoom=1:",
                                                      opt.comp$window),
                                               sep = " / "),
                                pvalue  = opt.comp$pvalue,
                                opt.mean = fres$opt.mean,
                                opt.hor = horizontal.ON )

            if (horizontal.ON == TRUE) {
              graphics::abline(v = mean_fct(fres$fobs, fres$opt.mean),
                               lty = "dashed", col = "red")
            } else {
              graphics::abline(h = mean_fct(fres$fobs, fres$opt.mean),
                               lty = "dashed", col = "red")
            }
          }
        }
      }

    } else {

      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      # plot all components
      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      plot_components_box(fres$fobs, fres$mOccur,
                          sres    = ref.all[indOrd, ],
                          ordre   = ordre,
                          main    = paste(main, "opt.comp", "tree",
                                          sep = " / "),
                          pvalue  = opt.comp$pvalue,
                          opt.mean = fres$opt.mean,
                          opt.hor = horizontal.ON )

      if ("leg" %in% opt.comp)
        plot_components_content(fres$mOccur, ordre = rev(ordre))

      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      # zoom on the first part of primary tree
      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      if ("zoom" %in% opt.comp) if (opt.comp$window < nbElt) {

        plot_components_box(fres$fobs, fres$mOccur,
                            sres    = ref.all[indOrd, ],
                            ordre   = ordre,
                            elt.wdw = 1:opt.comp$window,
                            main    = paste(main, "opt.comp", "tree",
                                           paste0("zoom=1:", opt.comp$window),
                                           sep = " / "),
                            pvalue  = opt.comp$pvalue,
                            opt.mean = fres$opt.mean,
                            opt.hor = horizontal.ON )
      }
    }
  }



  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  # components are sorted by increasing mean performance of components
  #               and plotted in horizontal boxplots
  #
  # horizontal.ON <- TRUE;  decreasing.ON <- FALSE
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if ( ("perf" %in% opt.comp) & (horizontal.ON == TRUE) ) {

    ordre <- as.integer(ref.all$motif)

    if (xpr.ON == TRUE) {

      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      # plot all components over all experiments
      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      plot_components_box(fres$fobs, fres$mOccur,
                          sres    = ref.all,
                          ordre   = ordre,
                          main    = paste(main, "opt.comp", "perf",
                                          "xpr=all", sep = " / "),
                          pvalue  = opt.comp$pvalue,
                          opt.mean = fres$opt.mean,
                          opt.hor = horizontal.ON)

      if ("leg" %in% opt.comp)
        plot_components_content(
          fres$mOccur[which(names(fres$xpr) == setXpr[ipr]), , drop = FALSE],
          ordre = ordre)

      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      # plot all components for each experiment
      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      if ( ("xpr" %in% opt.comp) | ("multi" %in% opt.comp) ) {

        for (ipr in seq_along(setXpr)) {

          index <- which(names(fres$xpr) == setXpr[ipr])
          plot_components_box(fres$fobs[index],
                              fres$mOccur[index, , drop = FALSE],
                              sres    = ref.xpr[[ipr]],
                              ordre   = ordre,
                              main    = paste(main, "opt.comp", "perf",
                                              setXpr[ipr], sep = " / "),
                              pvalue  = opt.comp$pvalue,
                              opt.mean = fres$opt.mean,
                              opt.hor = horizontal.ON )

          if (horizontal.ON == TRUE) {
            graphics::abline(v = mean_fct(fres$fobs, fres$opt.mean),
                             lty = "dashed", col = "red")
          } else {
            graphics::abline(h = mean_fct(fres$fobs, fres$opt.mean),
                             lty = "dashed", col = "red")
          }
        }
      }

      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      # zoom on the first part over all experiments
      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      if ("zoom" %in% opt.comp) if (opt.comp$window < nbElt) {

        plot_components_box(fres$fobs, fres$mOccur,
                            sres    = ref.all,
                            ordre   = ordre,
                            elt.wdw = rev(nbElt - 1:opt.comp$window + 1),
                            main    = paste(main, "opt.comp", "perf",
                                           "xpr=all",
                                           paste0("zoom=1:", opt.comp$window),
                                           sep = " / "),
                            pvalue  = opt.comp$pvalue,
                            opt.mean = fres$opt.mean,
                            opt.hor = horizontal.ON )

        #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
        # zoom on the first part of each experiment
        #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

        if ( ("xpr" %in% opt.comp) | ("multi" %in% opt.comp) ) {

          for (ipr in seq_along(setXpr)) {

            index <- which(names(fres$xpr) == setXpr[ipr])
            plot_components_box(fres$fobs[index],
                                fres$mOccur[index, , drop = FALSE],
                                sres    = ref.xpr[[ipr]],
                                ordre   = ordre,
                                elt.wdw = rev(nbElt - 1:opt.comp$window + 1),
                                main    = paste(main, "opt.comp", "perf",
                                                paste0("xpr=", setXpr[ipr]),
                                                paste0("zoom=1:",
                                                       opt.comp$window),
                                                sep = " / "),
                                pvalue  = opt.comp$pvalue,
                                opt.mean = fres$opt.mean,
                                opt.hor = horizontal.ON )

            if (horizontal.ON == TRUE) {
              graphics::abline(v = mean_fct(fres$fobs, fres$opt.mean),
                               lty = "dashed", col = "red")
            } else {
              graphics::abline(h = mean_fct(fres$fobs, fres$opt.mean),
                               lty = "dashed", col = "red")
            }
          }
        }
      }

    } else {

      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      # plot all components
      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      plot_components_box(fres$fobs, fres$mOccur,
                          sres    = ref.all,
                          ordre   = ordre,
                          main    = paste(main, "opt.comp", "perf",
                                          sep = " / "),
                          pvalue  = opt.comp$pvalue,
                          opt.mean = fres$opt.mean,
                          opt.hor = horizontal.ON )

      if ("leg" %in% opt.comp)
        plot_components_content(fres$mOccur, ordre = ordre)

      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      # zoom on the first part of primary tree
      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      if ("zoom" %in% opt.comp) if (opt.comp$window < nbElt) {

        plot_components_box(fres$fobs, fres$mOccur,
                            sres    = ref.all,
                            ordre   = ordre,
                            elt.wdw = rev(nbElt - 1:opt.comp$window + 1),
                            main    = paste(main, "opt.comp", "perf",
                                           paste0("zoom=1:", opt.comp$window),
                                           sep = " / "),
                            pvalue  = opt.comp$pvalue,
                            opt.mean = fres$opt.mean,
                            opt.hor = horizontal.ON )
      }
    }
  }



  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  #
  # components are sorted by decreasing mean performance of components
  #               and plotted in vertical boxplots
  #
  #  horizontal.ON <- FALSE;  decreasing.ON <- TRUE
  #
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if ( ("perf" %in% opt.comp) & (horizontal.ON == FALSE) ) {

    ordre   <- as.integer(ref.all$motif)

    if (xpr.ON == TRUE) {

      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      # plot all components over all experiments
      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      plot_components_box(fres$fobs, fres$mOccur,
                          sres    = ref.all,
                          ordre   = ordre,
                          main    = paste(main, "opt.comp", "perf",
                                          "xpr=all", sep = " / "),
                          pvalue  = opt.comp$pvalue,
                          opt.mean = fres$opt.mean,
                          opt.hor = horizontal.ON )

      if ("leg" %in% opt.comp)
        plot_components_content(
          fres$mOccur[which(names(fres$xpr) == setXpr[ipr]), , drop = FALSE],
          ordre = rev(ordre))

      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      # plot all components for each experiment
      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      if ( ("xpr" %in% opt.comp) | ("multi" %in% opt.comp) ) {

        for (ipr in seq_along(setXpr)) {

          index <- which(names(fres$xpr) == setXpr[ipr])
          plot_components_box(fres$fobs[index],
                              fres$mOccur[index, , drop = FALSE],
                              sres    = ref.xpr[[ipr]],
                              ordre   = ordre,
                              main    = paste(main, "opt.comp", "perf",
                                              setXpr[ipr], sep = " / "),
                              pvalue  = opt.comp$pvalue,
                              opt.mean = fres$opt.mean,
                              opt.hor = horizontal.ON )

          if (horizontal.ON == TRUE) {
            graphics::abline(v = mean_fct(fres$fobs, fres$opt.mean),
                             lty = "dashed", col = "red")
          } else {
            graphics::abline(h = mean_fct(fres$fobs, fres$opt.mean),
                             lty = "dashed", col = "red")
          }
        }
      }

      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      # zoom on the first part over all experiments
      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      if ("zoom" %in% opt.comp) if (opt.comp$window < nbElt) {

        plot_components_box(fres$fobs, fres$mOccur,
                            sres    = ref.all,
                            ordre   = ordre,
                            elt.wdw = 1:opt.comp$window,
                            main    = paste(main, "opt.comp", "perf",
                                            "xpr=all",
                                            paste0("zoom=1:", opt.comp$window),
                                            sep = " / "),
                            pvalue  = opt.comp$pvalue,
                            opt.mean = fres$opt.mean,
                            opt.hor = horizontal.ON )

        #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
        # zoom on the first part of each experiment
        #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

        if ( ("xpr" %in% opt.comp) | ("multi" %in% opt.comp) ) {

          for (ipr in seq_along(setXpr)) {

            plot_components_box(fres$fobs[index],
                                fres$mOccur[index, , drop = FALSE],
                                sres    = ref.xpr[[ipr]],
                                ordre   = ordre,
                                elt.wdw = 1:opt.comp$window,
                                main    = paste(main, "opt.comp", "perf",
                                                paste0("xpr=", setXpr[ipr]),
                                                paste0("zoom=1:",
                                                       opt.comp$window),
                                                sep = " / "),
                                pvalue  = opt.comp$pvalue,
                                opt.mean = fres$opt.mean,
                                opt.hor = horizontal.ON )

            if (horizontal.ON == TRUE) {
              graphics::abline(v = mean_fct(fres$fobs, fres$opt.mean),
                               lty = "dashed", col = "red")
            } else {
              graphics::abline(h = mean_fct(fres$fobs, fres$opt.mean),
                               lty = "dashed", col = "red")
            }
          }
        }
      }

    } else {

      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      # plot all components
      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      plot_components_box(fres$fobs, fres$mOccur,
                          sres    = ref.all,
                          ordre   = ordre,
                          main    = paste(main, "opt.comp", "perf",
                                          sep = " / "),
                          pvalue  = opt.comp$pvalue,
                          opt.mean = fres$opt.mean,
                          opt.hor = horizontal.ON )

      if ("leg" %in% opt.comp)
        plot_components_content(fres$mOccur, ordre = rev(ordre))

      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      # zoom on the first part of primary tree
      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      if ("zoom" %in% opt.comp) if (opt.comp$window < nbElt) {

        plot_components_box(fres$fobs, fres$mOccur,
                            sres    = ref.all,
                            ordre   = ordre,
                            elt.wdw = 1:opt.comp$window,
                            main    = paste(main, "opt.comp", "perf",
                                            paste0("zoom=1:", opt.comp$window),
                                            sep = " / "),
                            pvalue  = opt.comp$pvalue,
                            opt.mean = fres$opt.mean,
                            opt.hor = horizontal.ON )
      }
    }
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# End of file                                                              ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
