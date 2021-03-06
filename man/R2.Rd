% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats.R
\name{R2}
\alias{R2}
\title{Pearson' R2.}
\usage{
R2(x, y)
}
\arguments{
\item{x, y}{numeric vectors of same length.}
}
\value{
Return the "Pearson" coefficient of correlation
for the linear regression between the vectors \code{x} and \code{y}.
}
\description{
Take two numeric vectors and return the "Pearson"
coefficient of correlation for the linear regression
between \code{x} and \code{y}.
}
\seealso{
\code{\link{rss}}, \code{\link{mse}}, \code{\link{rmse}},
\code{\link{R2mse}}, \code{\link{pmse}}
}
\keyword{internal}
