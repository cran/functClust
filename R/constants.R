#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                                  CONSTANTS.R
#
# Set of constants used in the package functClust                          ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Value of epsilon used in \pkg{functClust}
#'
#' @description Return the constant value of epsilon.
#'
#' @usage fepsilon()
#'
#' @return Return the constant value of epsilon
#'
#' @details Epsilon is a constant, equals to \code{1e-12}.
#'
#' @seealso
#' \code{\link{fepsilon}} returns a value for epsilon. \cr
#' \code{\link{fletters}} returns a vector of lowercase
#'   then uppercase letters. \cr
#' \code{\link{fcolours}} returns a vector of colours. \cr
#' \code{\link{fsymbols}} returns a vector of fsymbols.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

fepsilon <- function() {

  return( 1e-12 )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Standard p-value used in the package functClust
#'
#' @description Return the p-value used in the package functClust
#'
#' @usage fpvalue()
#'
#' @return Return the p-value used in the package functClust.
#'
#' @details None.
#'
#' @seealso
#' \code{\link{fepsilon}}: return a value for epsilon. \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

fpvalue <- function() {

  return( 0.001 )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Number of digits kept when writing in a file
#'
#' @description Return the number of digits kept when writing in a file.
#'
#' @usage fnbdigits()
#'
#' @return Return the number of digits kept when writing in a file.
#'
#' @details None.
#'
#' @seealso
#' \code{\link{fepsilon}}: return a value for epsilon. \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

fnbdigits <- function() {

  return( 7 )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Series of symbols used in \pkg{functClust}
#'
#' @description Return the symbols used in the package \pkg{functClust}.
#'
#' @usage fsymbols()
#'
#' @return Return a vector of symbols
#'  used in the package \pkg{functClust}.
#'
#' @details The vector of symbols is \code{c(21, 22, 24, 23, 25)}.
#' It corresponds to: \code{21 = "circle"}, \code{22 = "square"},
#' \code{23 = "diamond"}, \code{24 = "triangle point-up"},
#' \code{25 = "triangle point-down"}. \cr
#'
#'  The use of functions \code{fsymbols() * fcolours()} generates
#'  a vector of \code{5 * 6 = 30} different symbols.
#'
#' @seealso
#' \code{\link{fepsilon}} returns a value for epsilon. \cr
#' \code{\link{fletters}} returns a vector of lowercase
#' then uppercase letters. \cr
#' \code{\link{fcolours}} returns a vector of colours. \cr
#' \code{\link{fsymbols}} returns a vector of symbols.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

fsymbols <- function() {

  return( c(21, 22, 24, 23, 25) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Series of colours used in \pkg{functClust}
#'
#' @description Return the colours used in the package \pkg{functClust}.
#'
#' @usage fcolours()
#'
#' @return Return the vector of the colours
#'  used in the package \pkg{functClust}.
#'
#' @details The vector of colours is
#' \code{c("red3", "blue2", "orange2", "turquoise3", "magenta", "green4")}.\cr
#'
#'  The use of functions \code{fsymbols() * fcolours()} generates
#'  a vector of \code{5 * 6 = 30} different symbols,
#'  before to be automatically recycled by R.
#'
#' @seealso
#' \code{\link{fepsilon}} returns a value for epsilon. \cr
#' \code{\link{fletters}} returns a vector of lowercase
#' then uppercase letters. \cr
#' \code{\link{fcolours}} returns a vector of colours. \cr
#' \code{\link{fsymbols}} returns a vector of symbols.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

fcolours <- function() {

  return( c("red3", "blue2", "orange2", "turquoise3", "magenta", "green4") )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Series of letters used in \pkg{functClust}
#'
#' @description Concatenate lowercase and uppercase letters
#' used in the package \pkg{functClust}.
#'
#' @usage fletters()
#'
#' @return Return the vector of letters
#'  used in the package \pkg{functClust}.
#'
#' @details The vector of letters is \code{c(letters, LETTERS)},
#' which corresponds to lowercase letters, followed by uppercase letters.
#'
#' @seealso
#' \code{\link{fepsilon}} returns a value for epsilon. \cr
#' \code{\link{fletters}} returns a vector of lowercase
#' then uppercase letters. \cr
#' \code{\link{fcolours}} returns a vector of colours. \cr
#' \code{\link{fsymbols}} returns a vector of fsymbols.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

fletters <- function() {

  return( c(letters, LETTERS) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Standard number of components to plot in the package functClust
#'
#' @description Return the standard number of components
#' to plot in the package functClust
#'
#' @usage fwindow()
#'
#' @return Return the standard number of components
#' to plot in the package functClust.
#'
#' @details None.
#'
#' @seealso
#' \code{\link{fepsilon}}: return a value for epsilon. \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

fwindow <- function() {

  return( 25 )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Standard colour used in the package functClust
#'
#' @description Return the standard colour used in the package functClust
#'
#' @usage fstd_colour()
#'
#' @return Return the standard colour used in the package functClust.
#'
#' @details None.
#'
#' @seealso
#' \code{\link{fepsilon}}: return a value for epsilon. \cr
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

fstd_colour <- function() {

  return( "black" )
}





#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# End of file                                                              ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
