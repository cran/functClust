#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                                  STATS.R
#
# Set of statistical functions                                             ####
#                           used in the package functClust
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#' @include
#' constants.R
#'



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Basic statistics                                                         ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic mean
#'
#' @description Take a numeric vector and return the arithmetic mean.
#'
#' @usage amean(x)
#'
#' @param x a numeric vector.
#'
#' @return Return the arithmetic mean of a numeric vector.
#'
#' @details The function \code{amean()} works as the standard function
#'  \code{stats::mean()}, but with option \code{na.rm = TRUE} by default.
#'
#' @seealso
#' \code{\link{amean}}  returns an arithmetic mean. \cr
#' \code{\link{gmean}}  returns a geometric mean. \cr
#' \code{\link{wamean}} returns a weighted arithmetic mean. \cr
#' \code{\link{wgmean}} returns a weighted geometric mean.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

amean <- function(x) {

  x <- x[!is.na(x)]

  return( ifelse(length(x), sum(x) / length(x), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic standard deviation
#'
#' @description Take a numeric vector
#' and return its arithmetic standard deviation.
#'
#' @usage asd(x)
#'
#' @param x a numeric vector.
#'
#' @return Return the arithmetic standard deviation of a numeric vector,
#'  as the standard function \code{stats::sd()}.
#'
#' @details The function \code{asd()} works as the standard function
#' \code{stats::sd()}, but with option \code{na.rm = TRUE} by default.
#'
#' @seealso
#' \code{\link{wasd}} returns arithmetic standard deviation. \cr
#' \code{\link{wasd}} returns weighted arithmetic standard deviation. \cr
#' \code{\link{gsd}}  returns geometric standard deviation. \cr
#' \code{\link{wgsd}} returns weighted geometric standard deviation.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

asd <- function(x) {

  x <- x[!is.na(x)]

  return( ifelse(length(x), sqrt( sum((x - amean(x)) ^ 2) / length(x) ), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Weighted arithmetical mean
#'
#' @description Take a numeric vector and
#' return the weighted arithmetic mean.
#'
#' @usage wamean(x, w)
#'
#' @param x a numeric vector.
#'
#' @param w a vector of numeric weights, of \code{length(x)}.
#'
#' @return Return the weighted arithmetic mean of a numeric vector.
#'
#' @details None.
#'
#' @seealso
#' \code{\link{amean}}  returns an arithmetic mean. \cr
#' \code{\link{gmean}}  returns a geometric mean. \cr
#' \code{\link{wamean}} returns a weighted arithmetic mean. \cr
#' \code{\link{wgmean}} returns a weighted geometric mean.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

wamean <- function(x, w) {

  index <- !(is.na(x) | is.na(w))
  x     <- x[index]
  w     <- w[index]

  return( ifelse(length(x), sum(w * x) / sum(w), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Weighted arithmetic standard deviation
#'
#' @description Take a numeric vector
#' and return the weighted arithmetic standard deviation.
#'
#' @usage wasd(x, w)
#'
#' @param x a numeric vector.
#'
#' @param w a vector of numeric weights, of \code{length(x)}.
#'
#' @return Return the weighted arithmetic standard deviation
#' of a numeric vector.
#'
#' @details None.
#'
#' @seealso
#' \code{\link{wasd}} returns arithmetic standard deviation. \cr
#' \code{\link{wasd}} returns weighted arithmetic standard deviation. \cr
#' \code{\link{gsd}}  returns geometric standard deviation. \cr
#' \code{\link{wgsd}} returns weighted geometric standard deviation.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

wasd <- function(x, w) {

  index <- !(is.na(x) | is.na(w))
  x     <- x[index]
  w     <- w[index]

  return( sqrt(ifelse(length(x) > 0,
                      sum( w * ((x - sum(x*w)/sum(w)) ^ 2) ) / sum(w),
                      NA)) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric mean
#'
#' @description Take a numeric vector and return its geometric mean.
#'
#' @usage gmean(x)
#'
#' @param x a numeric vector.
#'
#' @return Return the geometric mean of a numeric vector.
#'
#' @details None.
#'
#' @seealso
#' \code{\link{amean}} returns an arithmetic mean. \cr
#' \code{\link{gmean}} returns a geometric mean. \cr
#' \code{\link{wamean}} returns a weighted arithmetic mean. \cr
#' \code{\link{wgmean}} returns a weighted geometric mean.
#'

#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

gmean <- function(x) {

  x <- x[!is.na(x)]

  #  return( ifelse(length(x), prod(x) ^ (1/length(x)), NA) )
  return( ifelse(length(x), exp(amean(log(x))), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric standard deviation
#'
#' @description Take a numeric vector
#'  and return its geometric standard deviation.
#'
#' @usage gsd(x)
#'
#' @param x a numeric vector.
#'
#' @return Return the geometric standard deviation of a numeric vector.
#'
#' @details None.
#'
#' @seealso
#' \code{\link{wasd}} returns arithmetic standard deviation. \cr
#' \code{\link{wasd}} returns weighted arithmetic standard deviation. \cr
#' \code{\link{gsd}}  returns geometric standard deviation. \cr
#' \code{\link{wgsd}} returns weighted geometric standard deviation.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

gsd <- function(x) {

  x <- x[!is.na(x)]

  return( ifelse(length(x),
                 exp( sqrt( sum( log(x/gmean(x)) ^ 2 ) / length(x) ) ),
                 NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Weighted geometric mean
#'
#' @description Take a numeric vector and return its weighted geometric mean.
#'
#' @usage wgmean(x, w)
#'
#' @param x a numeric vector.
#'
#' @param w a vector of weights, of \code{length(x)}.
#'
#' @return Return the weighted geometric standard deviation
#' of a numeric vector.
#'
#' @details None.
#'
#' @seealso
#' \code{\link{amean}} returns an arithmetic mean. \cr
#' \code{\link{gmean}} returns a geometric mean. \cr
#' \code{\link{wamean}} returns a weighted arithmetic mean. \cr
#' \code{\link{wgmean}} returns a weighted geometric mean.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

wgmean <- function(x, w) {

  index <- !(is.na(x) | is.na(w))
  x     <- x[index]
  w     <- w[index]

  return( ifelse(length(x), prod(x ^ w) ^ (1 / sum(w)), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Weighted geometric standard deviation
#'
#' @description Take a numeric vector
#' and return its weighted geometric standard deviation.
#'
#' @usage wgsd(x, w)
#'
#' @param x a numeric vector.
#'
#' @param w a vector of numeric weights of \code{length(x)}.
#'
#' @return Return the weighted geometric standard deviation
#' of a numeric vector.
#'
#' @details None.
#'
#' @seealso
#' \code{\link{wasd}} returns arithmetic standard deviation. \cr
#' \code{\link{wasd}} returns weighted arithmetic standard deviation. \cr
#' \code{\link{gsd}}  returns geometric standard deviation. \cr
#' \code{\link{wgsd}} returns weighted geometric standard deviation.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

wgsd <- function(x, w) {

  index <- !(is.na(x) | is.na(w))
  x     <- x[index]
  w     <- w[index]

  return( ifelse(length(x), prod(x ^ w) ^ (1 / sum(w)), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Switch mean function.
#'
#' @description Switch for arithmetic or geometric mean according
#' to \code{opt.mean = "amean"} or \code{opt.mean = "gmean"}, respectively.
#'
#' @usage mean_fct(x, opt.mean = c("amean", "gmean"))
#'
#' @param x a numeric vector.
#'
#' @param opt.mean equals to \code{"amean"} or \code{"gmean"}
#' according to that mean value must be computed with arithmetic
#' or geometric formula, respectively. There is no default value.
#'
#' @return Return the arithmetic or geometric mean
#' of the numeric vector \code{x}
#' according to \code{opt.mean = "amean"} or
#' \code{opt.mean = "gmean"}.
#'
#' @details None.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

mean_fct <- function(x, opt.mean = c("amean", "gmean")) {

  return(switch(opt.mean,
                amean = amean(x),
                gmean = gmean(x)
  ))
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
#  Root mean square functions                                              ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#' @title Pearson' R2.
#'
#' @description Take two numeric vectors and return the "Pearson"
#' coefficient of correlation for the linear regression
#' between \code{x} and \code{y}.
#'
#' @usage R2(x, y)
#'
#' @param x,y numeric vectors of same length.
#'
#' @return Return the "Pearson" coefficient of correlation
#' for the linear regression between the vectors \code{x} and \code{y}.
#'
#' @importFrom stats var cov
#'
#' @seealso
#' \code{\link{rss}}, \code{\link{mse}}, \code{\link{rmse}},
#' \code{\link{R2mse}}, \code{\link{pmse}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

R2 <- function(x, y) {

  index <- !(is.na(x) | is.na(y))
  x     <- x[index]
  y     <- y[index]

  if (length(x) > 0) {
    z1 <- stats::var(x)
    z2 <- stats::var(y)
    z3 <- stats::cov(x,y)

    if ((z1 > 0) && (z2 > 0)) res <- (z3 ^ 2) / (z1 * z2)
  } else {
    res <- NA
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Residual Sum of Squares (\code{RSS}).
#'
#' @description Take two numeric vectors
#' and return their Residual Sum of Square (RSS).
#'
#' @usage rss(x, y)
#'
#' @param x,y two numeric vectors of same length.
#'
#' @return Return the Residual Sum of Square (RSS).
#'
#' @details None.
#'
#' @seealso \code{\link{R2}}, \code{\link{mse}}, \code{\link{rmse}},
#'          \code{\link{R2mse}}, \code{\link{pmse}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

rss <- function(x, y) {

  res <- sum((x - y) ^ 2)

  return( ifelse(res > fepsilon(), res, 0) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Mean Square Error (MSE).
#'
#' @description Take two numeric vectors
#' and return their Mean Square Error (\code{MSE}).
#'
#' @usage mse(x, y)
#'
#' @param x,y numeric vectors of same length.
#'
#' @return Return the Mean Square Error (\code{MSE}).
#'
#' @details None.
#'
#' @seealso \code{\link{R2}}, \code{\link{rss}}, \code{\link{rmse}},
#'          \code{\link{R2mse}}, \code{\link{pmse}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

mse <- function(x, y) {

  index <- !(is.na(x) | is.na(y))
  x     <- x[index]
  y     <- y[index]

  return( ifelse(length(x), rss(x, y) / length(x), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Root Mean Square Error (\code{RMSE})
#'
#' @description Take two numeric vectors
#' and return their Root Mean Square Error (\code{RMSE}).
#'
#' @usage rmse(x, y)
#'
#' @param x,y numeric vectors of same length.
#'
#' @return Return the Root Mean Square Error (\code{RMSE}).
#'
#' @details None.
#'
#' @seealso \code{\link{R2}}, \code{\link{rss}}, \code{\link{mse}},
#'          \code{\link{R2mse}}, \code{\link{pmse}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

rmse <- function(x, y) {

  return( sqrt(mse(x, y)) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Coefficient of Determination (\code{R2})
#'
#' @description Take two numeric vectors
#' and return their Coefficient of determination (\code{R2}).
#'
#' @usage R2mse(fprd, fobs)
#'
#' @param fprd,fobs two numeric vectors.
#' The first vector contains predicted or modelled values.
#' The second vector contains the reference (observed) values.
#'
#' @return Return the Coefficient of Determination (\code{R2}).
#'
#' @details Be careful, the function \code{\link{R2mse}} is not symmetrical.
#' The first argument \code{fprd} is vector of estimated, modelled or
#' predicted values. The second argument \code{fobs} is the vector
#' of reference. Both the vectors have the same length:
#'  \code{length(fprd) = length(fobs)}.
#'
#' @seealso \code{\link{R2}}, \code{\link{rss}}, \code{\link{mse}},
#'          \code{\link{rmse}}, \code{\link{pmse}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

R2mse <- function(fprd, fobs) {

  index <- !(is.na(fprd) | is.na(fobs))
  fprd   <- fprd[index]
  fobs   <- fobs[index]

  if (length(fobs) > 1) {
    tss <- sum((fobs - amean(fobs)) ^ 2)
    rss <- sum((fobs - fprd) ^ 2)
    res <- (tss - rss) / tss
  } else {
    res <- NA
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Probability associated with the Coefficient of determination (R2)
#'
#' @description Take two numeric vectors and return the
#'  Probability associated with their Coefficient of determination (R2).
#'
#' @usage pmse(fprd, fobs, nbK)
#'
#' @param fprd,fobs two numeric vectors.
#' The first vector contains predicted or modelled values.
#' The second vector contains the reference (observed) values.
#'
#' @param nbK an integer. \code{nbK > 1} : returns an
#'  error if not.
#'
#' @return Return the Probability associated with the Coefficient
#' of Ddtermination (R2).
#'
#' @details Be careful, the function \code{\link{pmse}} is not symmetrical.
#' The first argument \code{fprd} is vector of estimated, modelled
#' or predicted values, the second argument \code{fobs} is
#' the vector of reference (observed values). Both the vectors have the
#' same length: \code{length(fprd) = length(fobs)}.
#'
#' @importFrom stats pf
#'
#' @seealso \code{\link{R2}}, \code{\link{rss}}, \code{\link{mse}},
#'          \code{\link{rmse}}, \code{\link{R2mse}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

pmse <- function(fprd, fobs, nbK) {

  if (nbK <= 1) stop("nbK should be > 1")

  index <- !(is.na(fprd) | is.na(fobs))
  fprd   <- fprd[index]
  fobs   <- fobs[index]

  res <- 1
  if (length(fobs) > 1) {

    tss   <- sum( (fobs - amean(fobs)) ^ 2 )
    rss   <- sum( (fobs - fprd) ^ 2 )
    ess   <- tss - rss

    nbfdt <- length(fobs)
    nbfde <- nbK - 1
    nbfdr <- nbfdt - nbfde

    if ( (nbfdr > 0) && (nbfde > 0) ) {

      Fratio <- (ess / nbfde) / (rss / nbfdr)
      res    <- 1 - stats::pf(Fratio, nbfde, nbfdr)
    }
  }

  return(res)
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
#  AIC and related functions                                               ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title AIC of two numeric vectors
#'
#' @description Take two numeric vectors and return the Akaike Information
#'  Criterion (\code{AIC}).
#'  \code{AIC} is computed from the Residuals Sum of Square
#'  \code{\link{rss}} and
#'  the number of used parameters.
#'
#' @usage AIC_(x, y, nbK)
#'
#' @param x,y two numeric vectors of same length.
#'
#' @param nbK an integer. \code{nbK} is the number of used parameters.
#'  \code{nbK > 1} : returns \code{NA} if not.
#'
#' @return Return the Akaike Information Criterion (\code{AIC}).
#'
#' @details None.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

AIC_ <- function(x, y, nbK) {

  n   <- length(y)
  RSS <- rss(x, y)

  if (RSS > fepsilon()) {

    res <- n * log(rss(x, y) / n) + 2 * nbK
    if (res == -Inf) res <- NA

  } else {
    res <- NA
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title AICc of two numeric vectors
#'
#' @description Take two numeric vectors and return the Akaike Information
#'  Criterion (\code{AICc}) corrected for the second order.
#'  \code{AICc} is computed from the Residuals Sum of Square
#'  (\code{\link{rss}}) and
#'  the number of used parameters.
#'
#' @usage AICc(x, y, nbK)
#'
#' @param x,y two numeric vectors of same length.
#'
#' @param nbK an integer. \code{nbK} is the number of used parameters.
#'  \code{nbK > 1} : returns \code{NA} if not.
#'
#' @return Return the Akaike Information Criterion ((\code{AICc})),
#' corrected for the second order bias.
#'
#' @details The correction for the second order bias
#' improves the Akaike Information Criterion
#' when (\code{length(x)}) is small.
#'
#' @note The correction for the second order bias
#' improves the Akaike Information Criterion
#' when \code{length(x)} is small.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

AICc <- function(x, y, nbK) {

  n   <- length(y)
  RSS <- rss(x, y)

  if ( (RSS > fepsilon()) & (n + 1 > nbK) ) {

    res <- n * log(RSS / n) + 2 * nbK * (n + 2) / (n + 1 - nbK)
    if (res == -Inf) res <- NA

  } else {
    res <- NA
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Index of the minimum values of a vector
#'
#' @description Take a numeric vector and return the index of
#' its minimum values.
#'
#' @usage argmin(x)
#'
#' @param x a numeric vector.
#'
#' @return Return the index of the minimum values of a numeric vector.
#'
#' @details The function \code{argmin()} works as the standard
#' function \code{which.min()}.
#'
#' @seealso \code{\link{first_argmin}}, \code{\link{argmax}},
#'          \code{\link{first_argmax}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

argmin <- function(x) {

  x <- x[!is.na(x)]

  return( which(x == min(x, na.rm = TRUE)) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Index of the first minimum value of a vector
#'
#' @description Take a numeric vector and return the index of the first
#' (the lowest index) minimum values.
#'
#' @usage first_argmin(x)
#'
#' @param x a numeric vector.
#'
#' @return Return the index of the first (the lowest index)
#' minimum values of a numeric vector.
#'
#' @details The function \code{first_argmin()} works as the
#' standard function \code{min(which.min())}.
#'
#' @seealso \code{\link{argmin}}, \code{\link{argmax}},
#'          \code{\link{first_argmax}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

first_argmin <- function(x) {

  x <- x[!is.na(x)]

  arg <- 1
  while ((x[arg + 1] < x[arg]) & (arg < length(x))) arg <- arg + 1

  return(arg)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Index of the maximum values of a vector
#'
#' @description Take a numeric vector and return the index of its
#'  maximum values.
#'
#' @usage argmax(x)
#'
#' @param x a numeric vector.
#'
#' @return Return the index of the maximum values of a numeric vector.
#'
#' @details The function \code{argmax()} works as the standard function
#'  \code{which.max()}.
#'
#' @seealso \code{\link{argmin}}, \code{\link{first_argmin}},
#'          \code{\link{first_argmax}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

argmax <- function(x) {

  x <- x[!is.na(x)]

  return( which(x == max(x, na.rm = TRUE)) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Index of the first maximum value of a vector
#'
#' @description Take a numeric vector and return the index of
#' the first (the lowest index) maximum values.
#'
#' @usage first_argmax(x)
#'
#' @param x a numeric vector.
#'
#' @return Return the index of the first (the lowest index)
#'  maximum values of a numeric vector.
#'
#' @details The function \code{first_argmax()} works as the
#' standard function \code{min(which.max())}.
#'
#' @seealso \code{\link{argmin}}, \code{\link{first_argmin}},
#'          \code{\link{argmax}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

first_argmax <- function(x) {

  x <- x[!is.na(x)]

  arg <- 1
  while ((x[arg + 1] > x[arg]) & (arg < length(x))) arg <- arg + 1

  return(arg)
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Statistical functions of high level                                      ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Test posthoc of variance analysis
#'
#' @description Proceed to a variance analysis of the vector of data
#' \code{x} associated with the vector of factors \code{clusters}. It returns
#' a \code{data.frame} containing the name of the groups
#' \code{motif}, the group size \code{number}, the group mean
#' \code{mean} and standard deviation \code{sd}, and the two-to-two
#'  differences expressed in the form of letters \code{group}.
#'
#' @usage test_posthoc(x, clusters, pvalue = fpvalue())
#'
#' @param x a numeric vector.
#'
#' @param clusters a vector of factors of \code{length(x)}.
#'
#' @param pvalue a marginal p-value.
#'
#' @return Return a \code{data.frame} containing
#'  the name of the groups \code{motif},
#'   the group size \code{number},
#'  the group mean \code{mean} and standard deviation \code{sd},
#'  and the two-to-two differences expressed
#'  in the form of letters \code{group}.
#'
#' @details The function \code{test_posthoc()} uses Tukey method.
#' Different groups are sorted by decreasing means.
#' Letter rank increases with decreasing means.
#'
#' @importFrom stats aov TukeyHSD
#'
#' @importFrom multcompView multcompLetters4
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

test_posthoc <- function(x, clusters, pvalue = fpvalue()) {

  # check the inputs
  index <- is.na(x)
  if (sum(index) > 0) {
    clusters <- clusters[!index]
    x        <- x[!index]
  }

  numbers <- table(clusters)
  levels  <- sort(unique(clusters))
  nblevel <- length(levels)

  index   <- which(numbers > 1)
  if (length(index) > 1) {

    # compute the mean values of groups
    means <- sds <- numeric(nblevel)
    for (i in seq_len(nblevel)) means[i] <- amean(x[clusters == levels[i]])
    for (i in seq_len(nblevel)) sds[i]   <- asd(x[clusters == levels[i]])

    # analyse the variances
    if (!is.factor(clusters)) clusters <- as.factor(clusters)

    model <- stats::aov(x ~ clusters)
    test  <- stats::TukeyHSD(x = model, "clusters", conf.level = (1 - pvalue))
    vlet  <- multcompView::multcompLetters4(object = model, comp = test,
                                            threshold = pvalue,
                                            reversed = FALSE)

    sindex <- sort.list(means, decreasing = TRUE)
    res    <- data.frame(numbers[sindex], means[sindex], sds[sindex],
                         (vlet$clusters)$Letters)
    # Letters are always sorted out the function multcompLetters4

  } else {

    res   <- data.frame(levels, length(x), amean(x), asd(x), "a")
  }

  colnames(res) <- c("motif", "number", "mean", "sd", "group")
  res$motif     <- as.character(res$motif)
  rownames(res) <- NULL

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Test for the dependence of two Pearson' R2
#'
#' @description Test for the dependence of two R2.
#'
#' @usage test_dependent_R2(v1, v2, v3, n = length(v1))
#'
#' @param v1,v2,v3 three numeric vectors of \code{n = length(v1)},
#' or three coefficients of determination. If inputs are three coefficients
#'  of determination, vector length \code{n} must be specified.
#'
#' @param n the length of vectors. \code{n} must be specified
#' if the inputs are three coefficients of determination.
#'
#' @return Return a p-value.
#'
#' @details Be careful, the three vectors \code{v1, v2, v3}
#' are not symmetrical.
#' The function compares R2 obtained by the models \code{v1 ~ v3} and
#'  \code{v2 ~ v3}.
#'
#' @importFrom stats cor pt
#'
#' @seealso \code{\link{test_dependent_R2mse}},
#'          \code{\link{pvalue_dependent_R2mse}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

test_dependent_R2 <- function(v1, v2, v3, n = length(v1)) {

  if ((length(v1) > 1) || (length(v2) > 1) || (length(v3) > 1)) {

    if ((length(v1) != n) || (length(v2) != n) || (length(v3) != n))
      stop("length(v1) == length(v2) == length(v3) == n")

    r13 <- abs(stats::cor(v1, v3, method = "pearson"))
    r23 <- abs(stats::cor(v2, v3, method = "pearson"))
    r12 <- abs(stats::cor(v1, v2, method = "pearson"))

  } else {

    r13 <- sqrt(v1)
    r23 <- sqrt(v2)
    r12 <- sqrt(v3)
  }

  R      <- diag(1, 3)
  R[2, 1] <- R[1, 2] <- r13
  R[3, 1] <- R[1, 3] <- r23
  R[2, 3] <- R[3, 2] <- r12

  if ((r13 == r23) || (r12 == 1)) {
    p.value <- 1
  } else {
    tobs <- abs(r13 - r23) *
      sqrt((1 + r12) /
             (2 * abs(det(R)) / (n - 3) + (r13 + r23) ^ 2 * (1 - r12) ^ 3 /
                (4 * (n - 1))))
    p.value <- 2 * (1 - stats::pt(tobs, (n - 3)))
  }

  return(p.value)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Test for the dependence of two Coefficients of Determination R2
#'
#' @description Test for the dependence of two R2.
#'
#' @usage test_dependent_R2mse(fprd2, fprd1, fobs)
#'
#' @param fprd2,fprd1,fobs three numeric vectors of same length.
#'
#' @return Return a p-value.
#'
#' @details Be careful, the three vectors are not symmetrical.
#' The function compare R2 obtained by the models \code{fprd1 ~ fobs}
#' and \code{fprd2 ~ fobs}.
#'
#' @seealso \code{\link{test_dependent_R2}},
#'          \code{\link{pvalue_dependent_R2mse}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

test_dependent_R2mse <- function(fprd2, fprd1, fobs) {

  R1 <- R2mse(fprd2, fobs)
  R2 <- R2mse(fprd1, fobs)
  R3 <- R2mse(fprd2, fprd1)

  test1 <- ((R1 < 0) || (R2 < 0) || (R3 < 0))
  test2 <- (is.na(R1) || is.na(R2) || is.na(R3))
  test3 <- ( !((test1 || test2) == TRUE) )

  return( ifelse(test3, test_dependent_R2(R1, R2, R3, length(fobs)), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Test for the dependence of two R2
#'
#' @description Test for the dependence of two R2.
#'
#' @usage pvalue_dependent_R2mse(mprd, fobs)
#'
#' @param mprd a matrix of numeric vectors.
#'
#' @param fobs the reference vector of length equals to \code{dim(mprd)[2]}.
#'
#' @return Return a vector of p-value of length equals
#'  to \code{dim(mprd)[1]}.
#'
#' @details The function compare R2 obtained by the models
#'  \code{mprd[i, ] ~ fobs}.
#'
#' @seealso \code{\link{test_dependent_R2}},
#'          \code{\link{test_dependent_R2mse}}
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

pvalue_dependent_R2mse <- function(mprd, fobs) {

  nbline <- dim(mprd)[1]
  res    <- numeric(nbline)

  for (lin in 2:nbline) {

    index <- !is.na(mprd[lin, ])
    if (sum(index) > 0) {

      fprd  <- mprd[lin,     index]
      pprd <- mprd[lin - 1, index]
      pobs <- fobs[index]

    } else {

      fprd  <- mprd[lin,     ]
      pprd <- mprd[lin - 1, ]
      pobs <- fobs
    }

    R1 <- R2mse(fprd,  pobs)
    R2 <- R2mse(pprd, pobs)
    R3 <- R2mse(pprd, fprd)

    test1 <- ((R1 < 0) || (R2 < 0) || (R3 < 0))
    test2 <- (is.na(R1) || is.na(R2) || is.na(R3))

    res[lin] <- NA
    if (!((test1 || test2) == TRUE))
      res[lin] <- test_dependent_R2(R1, R2, R3, length(pobs))
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# End of file                                                              ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
