% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distributions-point.R
\name{point}
\alias{point}
\alias{dpoint}
\alias{rpoint}
\alias{ppoint}
\alias{qpoint}
\title{Point mass distribution}
\usage{
dpoint(x, location, log = FALSE)

rpoint(n, location)

ppoint(q, location, lower.tail = TRUE, log.p = FALSE)

qpoint(p, location, lower.tail = TRUE, log.p = FALSE)
}
\arguments{
\item{x, q}{vector or matrix of quantiles.}

\item{location}{vector of locations.}

\item{log, log.p}{logical; if \code{TRUE}, probabilities
\code{p} are given as \code{log(p)}.}

\item{n}{number of observations.}

\item{lower.tail}{logical; if \code{TRUE} (default), probabilities
are \eqn{P[X \le x]}, otherwise, \eqn{P[X \ge x]}.}

\item{p}{vector of probabilities.}
}
\value{
\code{dpoint} gives the density, \code{ppoint} gives the
distribution function, \code{qpoint} gives the quantile function,
and \code{rpoint} generates random deviates.
}
\description{
Density, distribution function,
quantile function and random generation for point distribution.
}
\examples{
# draw samples from a point distribution
rpoint(10, location = 1)

}
