#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                         BOOT_CLUSTERING.R
#
#  Functions for evaluating
#    the significance of data and robustess of functional clustering       ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#' @include
#'
NULL



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Test robustness of clustering by bootstrapping                          ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Evaluate the robustness of a functional clustering
#' to perturbations of data
#'
#' @description Evaluate by bootstrapping
#' the robustness of a functional clustering
#' to perturbations of data. The perturbed data can be
#' the number of assemblages taken into account,
#' or the number of performances taken into account.
#'
#' @usage
#' fboot_assemblages(fres,
#'                   opt.nbMax = fres$nbOpt, opt.R2 = FALSE, opt.plot = FALSE,
#'                   nbIter = 1, seed = NULL, rm.number = 0)
#'
#' @inheritParams fboot_one_point
#'
#' @details The trees obtained by bootstrapping of performances to omit
#'  are compared to the reference tree obtained with all components
#'  using different criteria :
#'  "Czekanowski_Dice", "Folkes_Mallows", "Jaccard", "Kulczynski",
#'  "Precision", "Rand", "Recall", "Rogers_Tanimoto", "Russel_Rao",
#'  "Sokal_Sneath1" and "Sokal_Sneath2" index.
#'  For more informations, see the notice of R-package \code{clusterCrit}.
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

fboot_assemblages <- function(fres, opt.nbMax = fres$nbOpt,
                              opt.R2 = FALSE, opt.plot = FALSE,
                              nbIter = 1, seed = NULL, rm.number = 0) {

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
  # Combinatorial analysis of data
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  indCrit   <- c(1:2,4:5,8:14)
  critNames <- clusterCrit::getCriteriaNames(FALSE)[indCrit]    #external crits
  vCriteria <- numeric(length(critNames))
  names(vCriteria) <- c(critNames)

  mCrit <- array(0, dim = c(fres$nbElt, nbIter, length(critNames)))
  dimnames(mCrit) <- list(seq_len(fres$nbElt), seq_len(nbIter), critNames)

  if (opt.R2 == TRUE) {
    mStats <- array(0, dim = c(fres$nbElt, nbIter, 2))
    dimnames(mStats) <- list(seq_len(fres$nbElt),
                             seq_len(nbIter), c("R2", "E"))
  }

  mIndAss <- build_random_matrix(nbIter, seed, rm.number,
                                 rm.number.max = fres$nbAss)

  setAss <- unique(rownames(fres$mOccur))
  for (iter in seq_len(nbIter)) {

    indAss <- mIndAss[iter, ]
    index  <- which(rownames(fres$mOccur) %in% setdiff(setAss, setAss[indAss]))

    main   <- paste0("The assemblages ", list_in_quote(setAss[indAss]),
                     " are removed")
    if (getOption("verbose") == TRUE) cat(main, "\n")

    tree.I   <- fit_ftree(fres$fobs[index], fres$mOccur[index, ],
                          fres$xpr[index], fres$affectElt,
                          fres$opt.method, fres$opt.mean, fres$opt.model,
                          opt.nbMax )

    for (nb in seq_len(opt.nbMax)) {
      vCriteria[] <- as.vector(unlist(clusterCrit::extCriteria(
        as.integer(cut_ftree(tree.I, nb) ),
        as.integer(lAffectElt[[nb]]),
        crit = "all")))[indCrit]
      mCrit[nb, iter, ] <- vCriteria
    }

    # compute the secondary tree

    if (opt.R2 == TRUE) {
      res <- validate_ftree(tree.I,
                            fres$fobs[index], fres$mOccur[index, ],
                            fres$xpr[index],
                            fres$opt.method,  fres$opt.mean, fres$opt.model,
                            fres$opt.jack, fres$jack,
                            opt.nbMax  )

      mStats[ , iter, "R2"] <- res$tStats[ , "R2cal"]
      mStats[ , iter, "E"]  <- res$tStats[ , "R2prd"]
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


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Organize the results
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  rboot <- vector(mode = "list", length = length(critNames))
  names(rboot) <- critNames
  for (ind in seq_along(critNames))
    rboot[[ind]] <- as.matrix(mCrit[ , , critNames[ind]])

  if (opt.R2 == TRUE) {
    mStats[fres$nbElt, , "R2"] <- mStats[fres$nbElt - 1, , "R2"]
    mStats[fres$nbElt, , "E"]  <- mStats[fres$nbElt - 1, , "E"]
    rboot$R2 <- as.matrix(mStats[ , , "R2"])
    rboot$E  <- as.matrix(mStats[ , , "E"])
  }


  return(rboot)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Evaluate the robustness of a functional clustering
#' to perturbations of data
#'
#' @description Evaluate by bootstrapping
#' the robustness  of a functional clustering
#' to perturbations of data.
#'
#' @usage
#' fboot_performances(fres,
#'                    opt.nbMax = fres$nbOpt, opt.R2 = FALSE, opt.plot = FALSE,
#'                    nbIter = 1, seed = NULL, rm.number = 0)
#'
#' @inheritParams fboot_one_point
#'
#' @details The trees obtained by bootstrapping of performances to omit
#'  are compared to the reference tree obtained with all components
#'  using different criteria :
#'  "Czekanowski_Dice", "Folkes_Mallows", "Jaccard", "Kulczynski",
#'  "Precision", "Rand", "Recall", "Rogers_Tanimoto", "Russel_Rao",
#'  "Sokal_Sneath1" and "Sokal_Sneath2" index.
#'  For more informations, see the notice of R-package \code{clusterCrit}.
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

fboot_performances <- function(fres, opt.nbMax = fres$nbOpt,
                               opt.R2 = FALSE, opt.plot = FALSE,
                               nbIter = 1, seed = NULL, rm.number = 0) {

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
  # Combinatorial analysis of data
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  indCrit   <- c(1:2,4:5,8:14)
  critNames <- clusterCrit::getCriteriaNames(FALSE)[indCrit]    #external crits
  vCriteria <- numeric(length(critNames))
  names(vCriteria) <- c(critNames)

  mCrit <- array(0, dim = c(fres$nbElt, nbIter, length(critNames)))
  dimnames(mCrit) <- list(seq_len(fres$nbElt), seq_len(nbIter), critNames)

  if (opt.R2 == TRUE) {
    mStats <- array(0, dim = c(fres$nbElt, nbIter, 2))
    dimnames(mStats) <- list(seq_len(fres$nbElt),
                             seq_len(nbIter), c("R2", "E"))
  }

  mIndXpr <- build_random_matrix(nbIter, seed, rm.number,
                                 rm.number.max = fres$nbXpr)

  setXpr <- unique(names(fres$xpr))
  for (iter in seq_len(nbIter)) {

    indXpr <- mIndXpr[iter, ]
    index  <- which(names(fres$xpr) %in% setdiff(setXpr, setXpr[indXpr]))

    main   <- paste0("The performances ", list_in_quote(setXpr[indXpr]),
                     " are removed")
    if (getOption("verbose") == TRUE) cat(main, "\n")

    tree.I   <- fit_ftree(fres$fobs[index], fres$mOccur[index, ],
                          fres$xpr[index], fres$affectElt,
                          fres$opt.method, fres$opt.mean, fres$opt.model,
                          opt.nbMax )

    for (nb in seq_len(opt.nbMax)) {
      vCriteria[] <- as.vector(unlist(clusterCrit::extCriteria(
        as.integer(cut_ftree(tree.I, nb) ),
        as.integer(lAffectElt[[nb]]),
        crit = "all")))[indCrit]
      mCrit[nb, iter, ] <- vCriteria
    }

    # compute the secondary tree

    if (opt.R2 == TRUE) {
      res <- validate_ftree(tree.I,
                            fres$fobs[index], fres$mOccur[index, ],
                            fres$xpr[index],
                            fres$opt.method,  fres$opt.mean, fres$opt.model,
                            fres$opt.jack, fres$jack,
                            opt.nbMax  )

      mStats[seq_len(fres$nbElt - 1), iter, "R2"] <- res$tStats[ , "R2cal"]
      mStats[seq_len(fres$nbElt - 1), iter, "E"]  <- res$tStats[ , "R2prd"]
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


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Organize the results
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  rboot <- vector(mode = "list", length = length(critNames))
  names(rboot) <- critNames
  for (ind in seq_along(critNames))
    rboot[[ind]] <- as.matrix(mCrit[ , , critNames[ind]])

  if (opt.R2 == TRUE) {
    mStats[fres$nbElt, , "R2"] <- mStats[fres$nbElt - 1, , "R2"]
    mStats[fres$nbElt, , "E"]  <- mStats[fres$nbElt - 1, , "E"]
    rboot$R2 <- as.matrix(mStats[ , , "R2"])
    rboot$E  <- as.matrix(mStats[ , , "E"])
  }


  return(rboot)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Input/Output boot functions                                             ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Record a test of significance of functional clustering
#'
#' @description Record in a file the results of
#' a test of significance of functional clustering.
#'
#' @usage
#' fboot_write_one_point(fres, rboot, filename)
#'
#' @inheritParams fboot
#'
#' @param rboot a list of matrices resulting from the functions
#'   \code{fboot_assemblages}, \code{fboot_performances}
#'   or \code{fboot_one_point}.
#'
#' @details The functions
#'   \code{fboot_assemblages}, \code{ftest_performances}
#'   and \code{fboot_one_point}
#'   generate a list containing a matrix by clustering index
#'   ("Czekanowski_Dice", "Folkes_Mallows", "Jaccard", "Kulczynski",
#'   "Precision", "Rand", "Recall", "Rogers_Tanimoto", "Russel_Rao",
#'   "Sokal_Sneath1" and "Sokal_Sneath2" index).
#'   Only their dimensions change according the used functions.
#'   Consequently, a same function is used for recording and reading
#'   the results of both the test-functions.
#'
#' @importFrom utils write.table
#'
#' @return Nothing. It is a procedure.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

fboot_write_one_point <- function(fres, rboot, filename) {

  # save random.matrix

  utils::write.table(
    x = cbind(rownames(rboot$random.matrix), rboot$random.matrix),
    file = paste0(delstr_end(filename, 3), "random.matrix.csv"),
    append = FALSE, col.names = TRUE, row.names = FALSE,
    sep = ",")

  rboot$random.matrix <- NULL


  # all matrices of Clustering Index are saved in a same file

  col1 <- NULL
  for (i in seq_along(names(rboot)))
    col1 <- c(col1, rep(names(rboot)[i], fres$nbElt))

  col2 <- rep(seq_len(fres$nbElt), length(names(rboot)))

  col3 <- NULL
  for (i in seq_along(names(rboot))) col3 <- rbind(col3, rboot[[i]])

  tmp <- cbind(col1, col2, col3)
  colnames(tmp) <- c("mat", "nbClu", seq_len(dim(rboot[[1]])[2]))
  rownames(tmp) <- as.character(seq_len(length(col1)))

  utils::write.table(x = tmp,
                     file = filename,
                     append = FALSE, col.names = TRUE, row.names = FALSE,
                     sep = ",")
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Read a test of significance of functional clustering
#'
#' @description Read a file of results obtained by
#' a test of significance of functional clustering.
#'
#' @usage
#' fboot_read_one_point(filename)
#'
#' @inheritParams fboot
#'
#' @details The functions
#'   \code{fboot_assemblages}, \code{ftest_performances}
#'   and \code{fboot_one_point}
#'   generate a list containing a matrix by clustering index
#'   ("Czekanowski_Dice", "Folkes_Mallows", "Jaccard", "Kulczynski",
#'   "Precision", "Rand", "Recall", "Rogers_Tanimoto", "Russel_Rao",
#'   "Sokal_Sneath1" and "Sokal_Sneath2" index).
#'   Only their dimensions change according the used functions.
#'   Consequently, a same function is used for recording and reading
#'   the results of both the test-functions.
#'
#' @importFrom utils read.table
#'
#' @return a list of matrices
#' each containing the results for a clustering index.
#'   The list is identical to this resulting from the functions
#'   \code{fboot_assemblages}, \code{fboot_performances}
#'   or \code{fboot_one_point}.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

fboot_read_one_point <- function(filename) {


  # read random.matrix

  file <- paste0(delstr_end(filename, 3), "random.matrix.csv")
  ttt  <- utils::read.table(file, header = TRUE, sep = ",",
                            check.names = FALSE, stringsAsFactors = FALSE)
  ttt  <- as.matrix(ttt[ , c(-1), drop = FALSE])
  rownames(ttt) <- seq_len(dim(ttt)[1])
  colnames(ttt) <- seq_len(dim(ttt)[2])

  # read all matrices saved in a same file

  tmp    <- utils::read.table(file = filename, header = TRUE, sep = ",",
                              check.names = FALSE, stringsAsFactors = FALSE)

  lnames <- unique(tmp[ , 1])
  index  <- c(1:2)


  # build the resulting list

  rboot <- vector(mode = "list", length = 1 + length(lnames))
  names(rboot) <- c("random.matrix", lnames)

  rboot$random.matrix <- ttt

  for (i in seq_along(lnames)) {
    rboot[[1 + i]] <- as.matrix(tmp[tmp$mat == lnames[i], -index])
    rownames(rboot[[1 + i]]) <- seq_len(dim(rboot[[1 + i]])[1])
    colnames(rboot[[1 + i]]) <- seq_len(dim(rboot[[1 + i]])[2])
  }

  return(rboot)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Main boot functions                                                     ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Build a matrix of missing indices
#'
#' @description Build a matrix of indices of elements hat will be ommitted,
#'  randomly drawn among a whole set of elements.
#'
#' @usage
#' build_random_matrix(nbIter = 1, seed = NULL,
#'                     rm.number = 0, rm.number.max = 0)
#'
#' @inheritParams fboot_one_point
#'
#' @param rm.number.max an integer,
#'  that indicates the maximum number of elements
#'  among which the elements to remove are randomly drawn.
#'
#' @details Each line (each set of indices) of the matrix
#' is unique, different from other lines of the matrix.
#'
#' @return a matrix of indices, of dimensions \code{nbIter} * \code{rm.number}.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

build_random_matrix <- function(nbIter = 1, seed = NULL,
                                rm.number = 0, rm.number.max = 0) {

  if (!is.null(seed)) set.seed(seed)

  mIndex <- matrix(0L, nrow = nbIter, ncol = rm.number,
                   dimnames = list(seq_len(nbIter), seq_len(rm.number)))

  mIndex[1, ] <- sort.int(sample.int(n = rm.number.max, size = rm.number))
  if (nbIter > 1) for (iter in 2:nbIter) {

    all.different <- TRUE
    while (all.different) {
      tmp <- sort.int(sample.int(n = rm.number.max, size = rm.number))
      names(tmp) <- colnames(mIndex)

      i <- 1
      while (!isTRUE(all.equal(tmp, mIndex[i, ])) & (i < iter)) i <- i + 1
      if (i == iter) all.different <- FALSE
    }

    mIndex[iter, ] <- tmp
  }

  return(mIndex)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Evaluate the robustness of a functional clustering
#' to perturbations of data
#'
#' @description Evaluate by bootstrapping
#' the robustness of a functional clustering
#' to perturbations of data. The perturbed data can be
#' the number of assemblages taken into account,
#' or the number of performances taken into account.
#'
#' @usage
#' fboot_one_point(fres,
#'                 opt.var = c("assemblages", "performances"),
#'                 nbIter = 1, rm.number = 0, seed = NULL,
#'                 opt.nbMax = fres$nbOpt, opt.R2 = FALSE,
#'                 opt.plot = FALSE )
#'
#' @inheritParams fboot
#'
#' @param seed an integer, that fixes a seed for random drawing.
#'
#' @param rm.number an integer,
#' that indicates the number of elements to randomly remove.
#'
#' @details The trees obtained by bootstrapping of performances to omit
#'  are compared to the reference tree obtained with all components
#'  using different criteria :
#'  "Czekanowski_Dice", "Folkes_Mallows", "Jaccard", "Kulczynski",
#'  "Precision", "Rand", "Recall", "Rogers_Tanimoto", "Russel_Rao",
#'  "Sokal_Sneath1" and "Sokal_Sneath2" index.
#'  For more informations, see the notice of R-package \code{clusterCrit}.
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

fboot_one_point <- function(fres,
                            opt.var = c("assemblages", "performances"),
                            nbIter = 1, rm.number = 0, seed = NULL,
                            opt.nbMax = fres$nbOpt, opt.R2 = FALSE,
                            opt.plot = FALSE ) {

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  # Check options
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  opt.var <- check_foption(opt.var, c("assemblages", "performances"))

  if (nbIter < 1) stop("'nbIter' should be higher than 1")

  if (opt.var == "assemblages")
    if ( (rm.number < 1) | (rm.number > fres$nbAss - 1) ) {
      stop("'nbIter' should be > 0 and < number of Assemblages")
    } else {
      namesTest <- rownames(fres$mOccur)
    }

  if (opt.var == "performances")
    if ( (rm.number < 1) | (rm.number > fres$nbXpr - 1) ) {
      stop("'nbIter' should be > 0 and < number of Performances")
    } else {
      namesTest <- names(fres$xpr)
    }


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

  mCrit <- array(0, dim = c(fres$nbElt, nbIter, length(critNames)))
  dimnames(mCrit) <- list(seq_len(fres$nbElt), seq_len(nbIter), critNames)

  if (opt.R2 == TRUE) {
    mStats <- array(0, dim = c(fres$nbElt, nbIter, 2))
    dimnames(mStats) <- list(seq_len(fres$nbElt),
                             seq_len(nbIter), c("R2", "E"))
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  setTest <- unique(namesTest)
  nbTest  <- length(setTest)

  rm.max <- choose(nbTest, rm.number)
  if (nbIter >= rm.max) nbIter <- rm.max

  mIndex  <- build_random_matrix(nbIter, seed, rm.number,
                                 rm.number.max = nbTest)

  for (iter in seq_len(nbIter)) {

    indSet <- mIndex[iter, ]
    main   <- paste("The", opt.var, list_in_quote(setTest[indSet]),
                    "are removed", sep = " ")
    if (getOption("verbose") == TRUE) cat(main, "\n")

    index  <- which(namesTest %in% setdiff(setTest, setTest[indSet]))
    tree.I <- fit_ftree(fres$fobs[index], fres$mOccur[index, ],
                        fres$xpr[index],  fres$affectElt,
                        fres$opt.method,  fres$opt.mean,
                        fres$opt.model,
                        opt.nbMax )

    for (nb in seq_len(opt.nbMax)) {
      vCriteria[] <- as.vector(unlist(clusterCrit::extCriteria(
        as.integer(cut_ftree(tree.I, nb) ),
        as.integer(lAffectElt[[nb]]),
        crit = "all")))[indCrit]
      mCrit[nb, iter, ] <- vCriteria
    }

    # compute the secondary tree

    if (opt.R2 == TRUE) {
      res <- validate_ftree(tree.I,
                            fres$fobs[index], fres$mOccur[index, ],
                            fres$xpr[index],
                            fres$opt.method,  fres$opt.mean, fres$opt.model,
                            fres$opt.jack, fres$jack,
                            opt.nbMax )

      mStats[ , iter, "R2"] <- res$tStats[ , "R2cal"]
      mStats[ , iter, "E"]  <- res$tStats[ , "R2prd"]
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


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Organize the results
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (opt.R2 == TRUE) {

    mStats[fres$nbElt, , "R2"] <- mStats[fres$nbElt - 1, , "R2"]
    mStats[fres$nbElt, , "E"]  <- mStats[fres$nbElt - 1, , "E"]

    rboot <- vector(mode = "list", length = 3 + length(critNames))
    names(rboot)        <- c("random.matrix", "R2", "E", critNames)
    rboot$random.matrix <- as.matrix(mIndex)
    rboot$R2            <- as.matrix(mStats[ , , "R2"])
    rboot$E             <- as.matrix(mStats[ , , "E"])
    for (ind in seq_along(critNames))
      rboot[[3 + ind]] <- as.matrix(mCrit[ , , critNames[ind]])

  } else {

    rboot <- vector(mode = "list", length = 1 + length(critNames))
    names(rboot) <- c("random.matrix", critNames)
    rboot$random.matrix <- as.matrix(mIndex)
    for (ind in seq_along(critNames))
      rboot[[1 + ind]] <- as.matrix(mCrit[ , , critNames[ind]])
  }

  return(rboot)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# End of file                                                              ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

