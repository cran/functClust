#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                             FCLUST.R
#
#  Tools for formatting input and output files                             ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#' @include
#' clustering.R
#' calibrating.R
#' validating.R
#'
NULL


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Formatting functions                                                    ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Check for identical assemblages
#'
#' @description Check if several assemblages
#' have the same elemental composition,
#' then average the performances of assemblages
#' with an identical elemental composition.
#'
#' @usage
#' check_repeat(fobs, mOccur, opt.mean = "amean", opt.repeat = FALSE )
#'
#' @inheritParams fit_ftree
#'
#' @inheritParams fclust
#'
#' @details None.
#'
#' @return Return a list containing:
#' \itemize{
#' \code{$fobs}{: the matrix of averaged performances of unique assemblages.}
#' \code{$mOccur}{: the matrix of occurrence of unique assemblages.}
#' }
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

check_repeat <- function(fobs, mOccur,
                         opt.mean = "amean", opt.repeat = FALSE ) {

  if (is.vector(fobs) == TRUE) fobs <- matrix(fobs, ncol = 1)

  assMotif <- affect_motifs(seq_len(dim(mOccur)[2]), mOccur)
  setDual  <- which(table(assMotif) > 1)

  bool <- !logical(length(assMotif))
  for (i in seq_along(setDual)) {
    dual <- which(assMotif == setDual[i])

    if (getOption("verbose") == TRUE)
      cat("The assemblages",
          list_in_quote(rownames(mOccur)[dual]),
          "have an identical elemental composition.\n")

    if (opt.repeat == FALSE) {
      fobs[dual[1], ] <- apply(fobs[dual, , drop = FALSE],
                               MARGIN = 2, FUN = mean_fct, opt.mean)
#      fobs[dual[1], ] <- apply(fobs[dual, , drop = FALSE],
#                               MARGIN = 2, FUN = median)
#      fobs[dual[1], ] <- apply(fobs[dual, , drop = FALSE],
#                               MARGIN = 2, FUN = max)
      bool[dual[2:length(dual)]] <- FALSE
    }
  }

  res        <- list(fobs[bool, ], mOccur[bool, ])
  names(res) <- c("fobs", "mOccur")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Remove components from dataset
#'
#' @description Remove one or more components from a dataset.
#'
#' @usage
#' remove_components(fobs, mOccur, components)
#'
#' @inheritParams fit_ftree
#'
#' @param components a list of strings,
#'   that indicates the components to remove from the matrix of occurrence
#'   and from vector (or matrix) of assemblage performances.
#'   The components to remove must belong to colnames of mOccur.
#'
#' @details None.
#'
#' @return Return a list containing:
#' \itemize{
#' \code{$fobs}    : the matrix of averaged performances of unique assemblages.
#' \code{$mOccur} : the matrix of occurrence of unique assemblages.
#' }
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

remove_components <- function(fobs, mOccur, components) {

  bool <- components %in% colnames(mOccur)
  if (sum(bool) != length(components))
    stop("The element ", paste(components[!bool], collapse = ", "),
         " does not belong to colnames(mOccur)")

  indRow <- indCol <- NULL
  for (elt in seq_along(components)) {
    indRow <- c(indRow, which(mOccur[ , components[elt]] == 1))
    indCol <- c(indCol, which(components[elt] == colnames(mOccur)))
  }

  mOccur <- mOccur[-indRow, -indCol]
  fobs   <- ifelse(is.vector(fobs), fobs[-indRow], fobs[-indRow, ])

  res        <- list(fobs, mOccur)
  names(res) <- c("fobs", "mOccur")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Input/output functions                                                  ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Read the file containing options for a functional clustering
#'
#' @description Read the file \code{"filename.options.csv"}
#'   containing the initial options
#'   for a functional clustering
#'   and saved by using the function \code{fclust_write()}.
#'
#' @usage
#' read_foptions(filename = "")
#'
#' @inheritParams read_fstats
#'
#' @details The results are saved in
#'   the file \code{"filename.options.csv"}.
#'   If the file does not exist or is corrupted, the function is stopped.
#'
#' @importFrom utils read.table
#'
#' @seealso
#'   \code{\link{read_foptions}}{:
#'     read the file \code{"filename.options.csv"}
#'     generated by \code{\link{fclust_write}}.} \cr
#'   \code{\link{read_finputs}}{:
#'     read the file \code{"filename.inputs.csv"}
#'     generated by \code{\link{fclust_write}}.} \cr
#'   \code{\link{read_ftrees}}{:
#'     read the file \code{"filename.trees.csv"}
#'     generated by \code{\link{fclust_write}}.} \cr
#'   \code{\link{read_fmatrices}}{:
#'     read the file \code{"filename.matrices.csv"}
#'     generated by \code{\link{fclust_write}}.} \cr
#'   \code{\link{read_fstats}}{:
#'     read the file \code{"filename.stats.csv"}
#'     generated by \code{\link{fclust_write}}.} \cr
#'
#' @return The initial options
#'   for a functional clustering
#'   and recorded in the files \code{"filename.options.csv"}.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

read_foptions <- function(filename = "") {

  if (nchar(filename) == 0) stop("'filename' cannot be an empty string")

  file <- paste(filename, "options", "csv", sep = ".")
  if (!file.exists(file)) stop("file '", file, "' does not exist.")

  tmp <- utils::read.table(file = file, header = TRUE, sep = ",")

  nbElt      <- as.integer(tmp$nbElt)
  nbAss      <- as.integer(tmp$nbAss)
  nbXpr      <- as.integer(tmp$nbXpr)

  opt.method <- as.character(tmp$opt.method)
  opt.mean   <- as.character(tmp$opt.mean)
  opt.model  <- as.character(tmp$opt.model)
  opt.jack   <- as.logical(  tmp$opt.jack)

  jack       <- integer(2)
  jack[1]    <- as.integer(tmp$jack.1)
  jack[2]    <- as.integer(tmp$jack.2)

  opt.na     <- as.logical(tmp$opt.na)
  opt.repeat <- as.logical(tmp$opt.repeat)
  affectElt        <- as.integer(tmp[12:length(tmp)])
  names(affectElt) <- as.character(names(tmp)[12:length(tmp)])

  res <- list(nbElt, nbAss, nbXpr,
              opt.method, opt.mean, opt.model,
              opt.jack, jack, opt.na, opt.repeat, affectElt)
  names(res) <- c("nbElt", "nbAss", "nbXpr",
                  "opt.method", "opt.mean", "opt.model",
                  "opt.jack", "jack", "opt.na", "opt.repeat", "affectElt")

  return(res)
}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Read the file containing initial data
#'   of occurrence and performances for a functional clustering
#'
#' @description Read the file \code{"filename.inputs.csv"}
#'   containing the initial data of assemblage occurrence and performances
#'   for a functional clustering
#'   and saved by using the function \code{fclust_write()}.
#'
#' @usage
#' read_finputs(filename = "", nbElt)
#'
#' @inheritParams read_fstats
#'
#' @details The results are saved in
#'   the file \code{"filename.inputs.csv"}.
#'   If the file does not exist or is corrupted, the function is stopped.
#'
#' @importFrom utils read.table
#'
#' @seealso
#'   \code{\link{read_foptions}}{:
#'     read the file \code{"filename.options.csv"}
#'     generated by \code{\link{fclust_write}}.} \cr
#'   \code{\link{read_finputs}}{:
#'     read the file \code{"filename.inputs.csv"}
#'     generated by \code{\link{fclust_write}}.} \cr
#'   \code{\link{read_ftrees}}{:
#'     read the file \code{"filename.trees.csv"}
#'     generated by \code{\link{fclust_write}}.} \cr
#'   \code{\link{read_fmatrices}}{:
#'     read the file \code{"filename.matrices.csv"}
#'     generated by \code{\link{fclust_write}}.} \cr
#'   \code{\link{read_fstats}}{:
#'     read the file \code{"filename.stats.csv"}
#'     generated by \code{\link{fclust_write}}.} \cr
#'
#' @return The initial data
#'   of assemblage occurrence and performances
#'   for a functional clustering
#'   and recorded in the files \code{"filename.inputs.csv"}.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

read_finputs <- function(filename = "", nbElt) {

  if (nchar(filename) == 0) stop("'filename' cannot be an empty string")

  file <- paste(filename, "inputs", "csv", sep = ".")
  if (!file.exists(file)) stop("file '", file, "' does not exist.")

  tmp <- utils::read.table(file = file, header = TRUE, sep = ",")

  fobs   <- as.numeric(tmp$fobs)
  xpr    <- as.numeric(tmp$xpr)
  mOccur <- as.matrix(tmp[ , 4 + 1:nbElt])
  storage.mode(mOccur) <- "integer"

  names(fobs) <- rownames(mOccur) <- as.character(tmp$names.fobs)
  names(xpr) <- as.character(tmp$names.xpr)
  colnames(mOccur) <- as.character(colnames(tmp)[4 + 1:nbElt])

  res <- list(fobs, mOccur, xpr)
  names(res) <- c("fobs", "mOccur", "xpr")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Read the file containing the trees
#'   resulting from a functional clustering
#'
#' @description Read the file \code{"filename.trees.csv"}
#'   containing the trees I and II resulting from a functional clustering
#'   and saved by using the function \code{fclust_write()}.
#'
#' @usage
#' read_ftrees(filename = "", nbElt)
#'
#' @inheritParams read_fstats
#'
#' @details The results are saved in
#'   the file \code{"filename.trees.csv"}.
#'   If the file does not exist or is corrupted, the function is stopped.
#'
#' @importFrom utils read.table
#'
#' @seealso
#'   \code{\link{read_foptions}}{:
#'     read the file \code{"filename.options.csv"}
#'     generated by \code{\link{fclust_write}}.} \cr
#'   \code{\link{read_finputs}}{:
#'     read the file \code{"filename.inputs.csv"}
#'     generated by \code{\link{fclust_write}}.} \cr
#'   \code{\link{read_ftrees}}{:
#'     read the file \code{"filename.trees.csv"}
#'     generated by \code{\link{fclust_write}}.} \cr
#'   \code{\link{read_fmatrices}}{:
#'     read the file \code{"filename.matrices.csv"}
#'     generated by \code{\link{fclust_write}}.} \cr
#'   \code{\link{read_fstats}}{:
#'     read the file \code{"filename.stats.csv"}
#'     generated by \code{\link{fclust_write}}.} \cr
#'
#' @return The trees I and II
#'   resulting from a functional clustering
#'   and recorded in the files \code{"filename.trees.csv"}.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

read_ftrees <- function(filename = "", nbElt) {

  if (nchar(filename) == 0) stop("'filename' cannot be an empty string")

  file <- paste(filename, "trees", "csv", sep = ".")
  if (!file.exists(file)) stop("file '", file, "' does not exist.")

  tmp <- utils::read.table(file = file, header = TRUE, sep = ",")

  nbOpt <- as.integer(tmp$nbOpt)[1]

  aff <- as.matrix(tmp[ , 2 + 1:nbElt])
  storage.mode(aff) <- "integer"
  cor <- as.numeric(tmp$R2.I)
  colnames(aff) <- delstr_end(colnames(aff), 2)
  rownames(aff) <- names(cor) <- as.character(tmp$nbClu)
  tree.I        <- list(aff, cor)
  names(tree.I) <- c("aff", "cor")

  aff <- as.matrix(tmp[ , 3 + nbElt + 1:nbElt])
  storage.mode(aff) <- "integer"
  cor <- as.numeric(tmp$R2.II)
  colnames(aff) <- delstr_end(colnames(aff), 3)
  rownames(aff) <- names(cor) <- as.character(tmp$nbClu)
  tree.II        <- list(aff, cor)
  names(tree.II) <- c("aff", "cor")

  res <- list(tree.I, tree.II, nbOpt)
  names(res) <- c("tree.I", "tree.II", "nbOpt")

  return(res)
}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Read the file containing the Matrices
#'   resulting from a functional clustering
#'
#' @description Read the file \code{"filename.matrices.csv"}
#'   containing the matrices resulting from a functional clustering
#'   and saved by using the function \code{fclust_write()}.
#'
#' @usage
#' read_fmatrices(filename = "", nbElt)
#'
#' @inheritParams read_fstats
#'
#' @details The results are saved in
#'   the file \code{"filename.matrices.csv"}.
#'   If the file does not exist or is corrupted, the function is stopped.
#'
#' @importFrom utils read.table
#'
#' @seealso
#'   \code{\link{read_foptions}}{:
#'     read the file \code{"filename.options.csv"}
#'     generated by \code{\link{fclust_write}}.} \cr
#'   \code{\link{read_finputs}}{:
#'     read the file \code{"filename.inputs.csv"}
#'     generated by \code{\link{fclust_write}}.} \cr
#'   \code{\link{read_ftrees}}{:
#'     read the file \code{"filename.trees.csv"}
#'     generated by \code{\link{fclust_write}}.} \cr
#'   \code{\link{read_fmatrices}}{:
#'     read the file \code{"filename.matrices.csv"}
#'     generated by \code{\link{fclust_write}}.} \cr
#'   \code{\link{read_fstats}}{:
#'     read the file \code{"filename.stats.csv"}
#'     generated by \code{\link{fclust_write}}.} \cr
#'
#' @return The matrices
#'   resulting from a functional clustering
#'   and recorded in the files \code{"filename.matrices.csv"}.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

read_fmatrices <- function(filename = "", nbElt) {

  if (nchar(filename) == 0) stop("'filename' cannot be an empty string")

  file <- paste(filename, "matrices", "csv", sep = ".")
  if (!file.exists(file)) stop("file '", file, "' does not exist.")

  tmp <- utils::read.table(file = file, header = TRUE, sep = ",",
                           check.names = FALSE, stringsAsFactors = FALSE)

  index   <- c(1:2)

  mCal    <- as.matrix(tmp[tmp$mat == "mCal",    -index])
  mPrd    <- as.matrix(tmp[tmp$mat == "mPrd",    -index])
  mMotifs <- as.matrix(tmp[tmp$mat == "mMotifs", -index])
  storage.mode(mMotifs) <- "integer"

  tCal    <- as.matrix(tmp[tmp$mat == "tCal",    -index])
  tPrd    <- as.matrix(tmp[tmp$mat == "tPrd",    -index])
  tNbcl   <- as.matrix(tmp[tmp$mat == "tNbcl",   -index])
  storage.mode(tNbcl) <- "integer"

  dimnames(mCal) <- dimnames(mPrd) <- dimnames(mMotifs) <-
    dimnames(tCal) <- dimnames(tPrd) <- dimnames(tNbcl) <-
    list(seq_len(nbElt), colnames(tmp)[-index])

  res <- list(mCal, mPrd, mMotifs, tCal, tPrd, tNbcl)
  names(res) <- c("mCal", "mPrd", "mMotifs", "tCal", "tPrd", "tNbcl")

  return(res)
}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Read the file containing the statistics of a functional clustering
#'
#' @description Read the file \code{"filename.stats.csv"}
#'   containing the statistics of a functional clustering
#'   and saved by using the function \code{fclust_write()}.
#'
#' @usage
#' read_fstats(filename = "", nbElt)
#'
#' @param filename a string, used as radical for the 6 file names.
#'
#' @param nbElt an integer, that indicates the number of components.
#'
#' @details The results are saved in
#'   the file \code{"filename.stats.csv"}.
#'   If the file does not exist or is corrupted, the function is stopped.
#'
#' @importFrom utils read.table
#'
#' @seealso
#'   \code{\link{read_foptions}}{:
#'     read the file \code{"filename.options.csv"}
#'     generated by \code{\link{fclust_write}}.} \cr
#'   \code{\link{read_finputs}}{:
#'     read the file \code{"filename.inputs.csv"}
#'     generated by \code{\link{fclust_write}}.} \cr
#'   \code{\link{read_ftrees}}{:
#'     read the file \code{"filename.trees.csv"}
#'     generated by \code{\link{fclust_write}}.} \cr
#'   \code{\link{read_fmatrices}}{:
#'     read the file \code{"filename.matrices.csv"}
#'     generated by \code{\link{fclust_write}}.} \cr
#'   \code{\link{read_fstats}}{:
#'     read the file \code{"filename.stats.csv"}
#'     generated by \code{\link{fclust_write}}.} \cr
#'
#' @return The statistics
#'   resulting from a functional clustering
#'   and recorded in the files \code{"filename.stats.csv"}.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

read_fstats <- function(filename = "", nbElt) {

  if (nchar(filename) == 0) stop("'filename' cannot be an empty string")

  file <- paste(filename, "stats", "csv", sep = ".")
  if (!file.exists(file)) stop("file '", file, "' does not exist.")

  tmp  <- utils::read.table(file = file, header = TRUE, sep = ",",
                            check.names = FALSE, stringsAsFactors = FALSE)

  index   <- c(1:2)

  mStats <- as.matrix(tmp[tmp$mat == "mStats", -index])
  tStats <- as.matrix(tmp[tmp$mat == "tStats", -index])

  dimnames(mStats) <- dimnames(tStats) <-
    list(seq_len(nbElt), colnames(tmp)[-index])

  res <- list(mStats, tStats)
  names(res) <- c("mStats", "tStats")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Main functions                                                          ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Make the formatted result from a functional clustering
#'
#' @description Check and homogeneize the format of results
#'    from a functional clustering,
#'    then make the formatted result.
#'
#' @usage
#' make_fclust(res, opt.na, opt.repeat, affectElt)
#'
#' @param res an object resulting from a functional clustering.
#'
#' @inheritParams fclust
#'
#' @details The results are saved in
#'   the file \code{"filename.options.csv"}.
#'   If the file does not exist or is corrupted, the function is stopped.
#'
#' @return \code{fres}, which is an object
#'     resulting from a functional clustering.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

make_fclust <- function(res, opt.na, opt.repeat, affectElt) {

  # Trees

  res$tree.I$cor  <- signif(res$tree.I$cor,  digits = fnbdigits())
  res$tree.II$cor <- signif(res$tree.II$cor, digits = fnbdigits())

  # Matrices

  res$mCal <- signif(res$mCal, digits = fnbdigits())
  res$mPrd <- signif(res$mPrd, digits = fnbdigits())
  storage.mode(res$mMotifs) <- "integer"

  res$tCal <- signif(res$tCal, digits = fnbdigits())
  res$tPrd <- signif(res$tPrd, digits = fnbdigits())
  storage.mode(res$tNbcl) <- "integer"

  # Stats

  res$mStats <- signif(res$mStats, digits = fnbdigits())
  res$tStats <- signif(res$tStats, digits = fnbdigits())


  # Add opt.na, opt.repeat, affectElt

  res$opt.na     <- opt.na
  res$opt.repeat <- opt.repeat

  if (res$opt.method == "apriori") {
    res$affectElt <- affectElt
  } else {
    res$affectElt <- cut_ftree(res$tree.I, res$nbOpt)
  }


  # Re-organize the list

  fres <- res

  begin <- which(names(res) == "jack") + 1
  end   <- which(names(res) == "opt.na")

  names(fres)[begin:(begin + 2)] <- names(res)[end:(end + 2)]
  fres[begin:(begin + 2)]        <- res[end:(end + 2)]

  names(fres)[(begin + 3):(end + 2)] <- names(res)[begin:(end - 1)]
  fres[(begin + 3):(end + 2)]        <- res[begin:(end - 1)]


  # Add the class "fclust"

#  attr(fres, "class") <- "fclust"

  return(fres)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Format a raw dataset for a functional clustering
#'
#' @description Fit a raw dataset,
#' then format it as n input for the function \code{fclust}.
#'
#' @usage
#' format_fclust(dat, nbElt, weight = rep(1, dim(dat)[2] - nbElt - 1),
#'               opt.na = FALSE, opt.repeat = FALSE, opt.mean   )
#'
#' @inheritParams fclust
#'
#' @details None.
#'
#' @return Return a list containing formatted inputs, that are:\cr
#'
#'  \itemize{
#'  \item \code{fobs}: the vector of observed performances of assemblages.
#'  \code{names(fobs)} contains the names of assemblages.
#'
#'  \item \code{mOccur}: the binary matrix of occurrence of components
#'  within the assemblages.
#'  \code{dim(mOccur) = [nbAss, nbElt]}.
#'  \code{rownames(mOccur)} contains the names of assemblages.
#'  \code{colnames(mOccur)} contains the names of components.
#'
#'  \item \code{xpr}: the vector of weight of different experiments.
#'  \code{names(xpr)} contains the names of different experiments.
#'  }
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

format_fclust <- function(dat, nbElt,
                          weight = rep(1, dim(dat)[2] - nbElt - 1),
                          opt.na = FALSE, opt.repeat = FALSE,
                          opt.mean = "amean"   ) {

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Check the format of inputs
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  nbAss  <- as.integer(dim(dat)[1])
  indRow <- seq_len(nbAss)

  nbXpr <- as.integer(dim(dat)[2] - nbElt - 1)

  if (nbXpr < 1)
    stop("The data.frame 'dat' should include assemblage performances.")


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Check the vector Weight
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (length(weight) != nbXpr)
    stop("Length of 'weight' should equal to performance number.")

  if ( !is.numeric(weight) || (sum(is.na(weight)) > 0) )
    stop("The vector 'weight' should be numeric.")


  # check Assemblage identity
  id    <- as.character(dat[ , 1])
  index <- which(table(id) > 1)
  if (length(index) > 0)
    stop(paste0("The assemblages indiced ",
                paste0(which(id == names(index)), collapse = ","),
                " have identical identity."))


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # check the Matrix of Occurrence
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  mOccur <- as.matrix(dat[ , 1 + (1:nbElt)],
                      ncol = nbElt, byrow = TRUE)
  storage.mode(mOccur) <- "integer"
  rownames(mOccur) <- id

  # binary values
  index <- table(mOccur)
  if ( (length(index) != 2) ||
       (sum(!(c("0", "1") %in% names(index))) != 0) )
    stop("The matrix of occurrence 'mOccur' should be binary.")

  # empty lines
  index <- which(apply(mOccur, MARGIN = 1, FUN = sum, na.rm = TRUE) == 0)
  if (length(index) > 0) {
    if (opt.na == TRUE) indRow <- indRow[-index]
    else stop(paste0("The assemblages ",
                     list_in_quote(rownames(mOccur)[index]),
                     " contain no components."))
  }

  # empty columns
  index <- which(apply(mOccur, MARGIN = 2, FUN = sum, na.rm = TRUE) == 0)
  if (length(index) > 0)
    stop(paste0("In occurrence matrix, the component ",
                concat_by_line(paste0("'", names(index), "'")),
                " occurs in no assemblage."))

  # 'NA' in lines
  index <- which(is.na(mOccur), arr.ind = TRUE)
  if (length(index) > 0) {
    if (opt.na == TRUE) indRow <- indRow[-index[ , "row"]]
    else stop(paste0("The occurrence matrix of assemblages ",
                     list_in_quote(rownames(mOccur)[index]),
                     " are not completely filled."))
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Check the matrix of Performances
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  fobs <- as.matrix(dat[ , (2 + nbElt):dim(dat)[2]],
                    ncol = nbXpr, mode = "numeric")
  rownames(mOccur) <- id

  # fobs is not numeric
  if (!is.numeric(fobs)) stop("The performances should be numeric")


  # fobs does not vary
  index <- NULL
  for (xpr in seq_len(nbXpr))
    if (length(unique(fobs[ , xpr])) == 1) index <- c(index, xpr)

  if (!is.null(index))
    stop(paste("The performances", list_in_quote(colnames(fobs)[index]),
                 "do not vary with assemblages", sep = " "))


  # 'NA' in fobs
  index <- which(is.na(fobs), arr.ind = TRUE)
  if (length(index) > 0) {
    if (opt.na == TRUE) indRow <- indRow[-index[ , "row"]]
    else stop(paste0("The performances of assemblages ",
                     list_in_quote(rownames(fobs)[index]),
                     " are not completely filled."))
  }


  # several assemblages have same elemental composition
  tmp    <- check_repeat(fobs[indRow, ], mOccur[indRow, ],
                         opt.mean, opt.repeat)
  mOccur <- tmp$mOccur
  fobs   <- as.matrix(tmp$fobs, ncol = nbXpr, mode = "numeric")
  rownames(fobs) <- rownames(mOccur)


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Formating of outputs
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  nbAss  <- as.integer(dim(mOccur)[1])

  MOCCUR <- matrix(0, nrow = nbXpr * nbAss, ncol = nbElt,
                   dimnames = list(rep(rownames(mOccur), nbXpr),
                                   colnames(mOccur)))
  storage.mode(MOCCUR) <- "integer"
  FOBS   <- XPR <- numeric(nbXpr * nbAss)

  setXpr <- as.character(colnames(dat)[(nbElt + 2):dim(dat)[2]])

  for (ipr in seq_along(setXpr)) {

    index <- (ipr - 1) * nbAss + 1:nbAss

    MOCCUR[index, ] <- mOccur

    FOBS[index]        <- fobs[ , ipr]
    names(FOBS)[index] <- rownames(mOccur)

    XPR[index]        <- rep(weight[ipr], nbAss)
    names(XPR)[index] <- rep(setXpr[ipr], nbAss)
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Outputs
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  res         <- list(FOBS, MOCCUR, XPR)
  names(res ) <- c("fobs", "mOccur", "xpr")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# End of file                                                              ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
