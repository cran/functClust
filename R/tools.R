#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                                  TOOLS.R
#
# Set of Constants and Tools                                               ####
#                 used in the package functClust
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#' @include
#' constants.R
#'
NULL



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Manipulating strings                                                     ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Index of different instances of a given char in a string
#'
#' @description Take a string \code{str}, then look for
#'  the index of the different instances of the char \code{chr}.
#'
#' @usage argchr(str, chr)
#'
#' @param str a string.
#'
#' @param chr a char of \code{nchar(chr) = 1}.
#'
#' @return Return a vector of indices
#' of the different instances of the character \code{chr} in \code{str}.
#'
#' @details None.
#'
#' @seealso
#' \code{\link{argchr}} returns all the instances
#' of \code{chr} in \code{str}. \cr
#' \code{\link{first_argchr}} returns the first instance
#' of \code{chr} in \code{str}. \cr
#' \code{\link{delstr_begin}} deletes the beginning of a string. \cr
#' \code{\link{delstr_end}} deletes the end of a string. \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

argchr <- function(str, chr) {

  if (nchar(chr) > 1) stop("nchar(chr) == 1")

  v <- NULL
  for (i in seq_len(nchar(str)))
    if (substr(str, i, i) == chr) v <- c(v, i)

  return(v)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Index of a first instance of char in a string
#'
#' @description Take a string \code{str}, then look for
#'  the position index of the first instance of the char \code{chr}.
#'
#' @usage first_argchr(str, chr)
#'
#' @param str a string.
#'
#' @param chr a char of \code{length(chr) = 1}.
#'
#' @return Return the position index of the first \code{chr} in \code{str}.
#'
#' @details None.
#'
#' @seealso
#' \code{\link{argchr}} returns all the instances
#' of \code{chr} in \code{str}. \cr
#' \code{\link{first_argchr}} returns the first instance
#' of \code{chr} in \code{str}. \cr
#' \code{\link{delstr_begin}} deletes the beginning of a string. \cr
#' \code{\link{delstr_end}} deletes the end of a string. \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

first_argchr <- function(str, chr) {

  if (nchar(chr) > 1) stop("nchar(chr) == 1")

  i <- 1
  while ( (substr(str, i, i) != chr) & (i < nchar(str)) ) i <- i + 1
  if (i == nchar(str)) i <- 0

  return(i)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Delete the beginning of a string
#'
#' @description Take a string, then delete the begining of this string.
#'
#' @usage
#' delstr_begin(str, n)
#'
#' @param str a string of length \code{nchar(str)}.
#'
#' @param n an integer.
#' It is the number of characters to delete
#' at the beginning of the string \code{str}
#'
#' @details None.
#'
#' @return A string of length \code{nchar(str) - n},
#' shortened by its beginning.
#'
#' @seealso
#' \code{\link{argchr}} returns all the instances
#' of \code{chr} in \code{str}. \cr
#' \code{\link{first_argchr}} returns the first instance
#' of \code{chr} in \code{str}. \cr
#' \code{\link{delstr_begin}} deletes the beginning of a string. \cr
#' \code{\link{delstr_end}} deletes the end of a string. \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

delstr_begin <- function(str, n) {

  return(substr(str, n + 1, nchar(str)))
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Delete the end of a string
#'
#' @description Take a string, then delete the end of this string.
#'
#' @usage delstr_end(str, n)
#'
#' @param str a string of length \code{nchar(str)}.
#'
#' @param n an integer.
#' It is the number of characters to delete
#' at the end of the string \code{str}
#'
#' @details None.
#'
#' @return A string of length \code{nchar(str) - n},
#' shortened by its end.
#'
#' @seealso
#' \code{\link{argchr}} returns all the instances
#' of \code{chr} in \code{str}. \cr
#' \code{\link{first_argchr}} returns the first instance
#' of \code{chr} in \code{str}. \cr
#' \code{\link{delstr_begin}} deletes the beginning of a string. \cr
#' \code{\link{delstr_end}} deletes the end of a string. \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

delstr_end <- function(str, n) {

  return(substr(str, 1, nchar(str) - n))
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Concat a list of strings
#'
#' @description Take a list of strings,
#' quotes each string,
#' then concat all quoted strings in a single long string
#' by separating each one by the given separator \code{sepChar},
#' then segment this long string into strings of a given length \code{nChar}.
#'
#' @usage
#' list_in_quote(str, sepChar = ",")
#'
#' @param str a list or vector of strings or numeric.
#'
#' @param sepChar a string.
#' The string is a separator between each element of list of strings.
#'
#' @details None.
#'
#' @return Return a vector of strings between quoted and separated by commas.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

list_in_quote <- function(str, sepChar = ",") {

  sep <- paste0("'", sepChar, "'")
  return( paste0("'", paste(str, collapse = sep), "'") )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Concat a list of strings segmented by line
#'
#' @description Take a list of strings,
#' concat all the strings in a single long string
#' by separating each one by the given separator \code{sepChar},
#' then segment this long string into strings of a given length \code{nChar}.
#'
#' @usage
#' concat_by_line(str, sepChar = ",", nbChar = 35)
#'
#' @param str a list of strings.
#'
#' @param sepChar a string.
#' The string is a separator between each element of list of strings.
#'
#' @param nbChar an integer.
#' It is the length of each line of strings.
#'
#' @details None.
#'
#' @return Return a vector of strings,
#' the length of each string being lower than \code{nbChar}.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

concat_by_line <- function(str, sepChar = ",", nbChar = 35) {

  res  <- NULL

  i <- 1
  while (i <= length(str)) {

    if (nchar(str[i]) < nbChar) {

      lin <- ""
      while ((nchar(lin) + nchar(str[i]) < nbChar) & (i <= length(str)) ) {
        lin <- paste(lin, str[i], sep = sepChar)
        i   <- i + 1
      }

      if (i <= length(str)) lin <- paste0(lin, sepChar)

      res <- c(res, delstr_begin(lin, nchar(sepChar)))

    } else {

      lin <- str[i]

      if (i <= length(str)) lin <- paste0(lin, sepChar)

      res <- c(res, lin)
      i   <- i + 1
    }
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Check an option
#'
#' @description Take a string which is an option,
#'  check if the string is among a given set of options,
#'  then stop the program if option is wrong,
#'  or return the right and complete option.
#'
#' @usage
#' check_foption(opt, all)
#'
#' @param opt a string, that is a possible otion.
#'
#' @param all a vector of strings, that iS a set of possible options.
#'
#' @details None.
#'
#' @return Stop the program, or return a string which is a right option.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

check_foption <- function(opt, all) {

  if (is.na(pmatch(opt, all)))
        stop("'", opt, "' should be ", list_in_quote(all))

  return( match.arg(opt, all) )
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Manipulating index, vector or matrix                                     ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Reverse an indexation
#'
#' @description Take a vector of indices and return the vector of indices
#' of \code{length(index)} that allows to find in return the initial order.
#'
#' @usage index_inturn(index)
#'
#' @param index a vector of integers.
#'
#' @return Return a vector of integers that allows
#' to find in return the initial order.
#'
#' @details None.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

index_inturn <- function(index) {

  res <- index
  for (i in seq_along(index)) res[index[i]] <- i

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Extends a vector
#'
#' @description Take the vector \code{x},
#' then extend it until to \code{length(x) = long}.
#'
#' @usage extend_vector(x, long)
#'
#' @param x any vector.
#'
#' @param long an integer, that corresponds
#' to the whished \code{length(x) = long}.
#'
#' @return Return a vector of \code{length(x) = long}
#' built by recycling \code{x} as many times as necessary,
#' then truncated at \code{length(x) = long}.
#'
#' @details None.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

extend_vector <- function(x, long) {

  lambda <- ceiling(long / length(x))
  res    <- integer(lambda * length(x))
  res[ ] <- x

  return( res[seq_len(long)] )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Extends a vector of letters
#'
#' @description Take the vector of lowercase and uppercase letters,
#' then extend it until to the length of vector equal to \code{long}.
#'
#' @usage extend_letters(long)
#'
#' @param long an integer, that corresponds
#' to the whished \code{length(vector) = long}.
#'
#' @return Return a vector of \code{length(vector) = long}
#' built by concatening lowercase and uppercase letters with an integer,
#' until \code{length(vector) = long}.
#'
#' @details None.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

extend_letters <- function(long) {

  str <- fletters()

  lambda <- ceiling(long / length(str))
  res    <- character(lambda * length(str))

  for (i in seq_len(lambda))
    res[(i - 1) * length(str) + 1:length(str) ] <- paste0(str, i)

  return( res[seq_len(long)] )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Reverse a table along its first dimension
#'
#' @description Take a table, then reverse it along is first dimension,
#' that is by moving the first rows as the last ones,
#' and the last rows as the first ones.
#' The function is useful for re-ordering
#' an decreasing table into an increasing table, or the inverse.
#'
#' @usage reverse_table(mat)
#'
#' @param mat a table or a matrix.
#'
#' @details None.
#'
#' @return Return a reversed table or matrix.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

reverse_table <- function(mat) {

  for (j in seq_len(dim(mat)[2])) mat[ , j] <- rev(mat[ , j])

  return(mat)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Test if a vector is binary
#'
#' @description Take a vector or a matrix,
#'   then check if all elements are equal to "0 or "1".
#'
#' @usage is_binary(mat)
#'
#' @param mat a binary/numeric matrix.
#'
#' @return Return \code{TRUE} is the matrix is a binary matrix.
#'
#' @details None.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

is_binary <- function(mat) {

  res <- table(mat)
  return(
    ifelse((length(res) == 2) & (sum(!(names(res) %in% c("0", "1"))) == 0),
           TRUE, FALSE )    )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Convert a binary into logical matrix
#'
#' @description Take the binary matrix \code{mat} and convert it
#'   in a logical matrix.
#'
#' @usage binary_to_logical(mat)
#'
#' @param mat a binary/numeric matrix.
#'
#' @return Return a logical matrix of \code{dim(mat)}.
#'
#' @details None.
#'
#' @seealso
#' \code{\link{logical_to_binary}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

binary_to_logical <- function(mat) {

  return( matrix(as.logical(mat), nrow = dim(mat)[1], ncol = dim(mat)[2],
                 dimnames = list(rownames(mat), colnames(mat)),
                 byrow = FALSE) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Convert a logical into binary matrix
#'
#' @description Take the logical matrix \code{mat} and convert it
#'     in a binary matrix.
#'
#' @usage logical_to_binary(mat)
#'
#' @param mat a logical matrix.
#'
#' @return Return a binary matrix of \code{dim(mat)}.
#'
#' @details None.
#'
#' @seealso
#' \code{\link{binary_to_logical}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

logical_to_binary <- function(mat) {

  tmp <- matrix(0L, nrow = dim(mat)[1], ncol = dim(mat)[2],
                dimnames = list(rownames(mat), colnames(mat)))
  tmp[mat == TRUE] <- 1L

  return(tmp)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Convert a character into integer vector
#'
#' @description Take the charactervector \code{v} and convert it
#'   in an integer vector.
#'
#' @usage char_to_int(v)
#'
#' @param v a character vector.
#'
#' @return Return a integer vector.
#'
#' @details None.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

char_to_int <- function(v) {

  res <- integer(length(v))
  set <- unique(v)

  for (i in seq_along(set)) res[v == set[i]] <- i

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# End of file                                                              ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
