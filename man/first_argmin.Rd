% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats.R
\name{first_argmin}
\alias{first_argmin}
\title{Index of the first minimum value of a vector}
\usage{
first_argmin(x)
}
\arguments{
\item{x}{a numeric vector.}
}
\value{
Return the index of the first (the lowest index)
minimum values of a numeric vector.
}
\description{
Take a numeric vector and return the index of the first
(the lowest index) minimum values.
}
\details{
The function \code{first_argmin()} works as the
standard function \code{min(which.min())}.
}
\seealso{
\code{\link{argmin}}, \code{\link{argmax}},
         \code{\link{first_argmax}}
}
\keyword{internal}
