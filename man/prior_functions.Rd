% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/priors.R
\name{prior_functions}
\alias{prior_functions}
\alias{rng.prior}
\alias{cdf.prior}
\alias{ccdf.prior}
\alias{lpdf.prior}
\alias{pdf.prior}
\alias{quant.prior}
\alias{mcdf.prior}
\alias{mccdf.prior}
\alias{mlpdf.prior}
\alias{mpdf.prior}
\alias{mquant.prior}
\title{Elementary prior related functions}
\usage{
\method{rng}{prior}(x, n, ...)

\method{cdf}{prior}(x, q, ...)

\method{ccdf}{prior}(x, q, ...)

\method{lpdf}{prior}(x, y, ...)

\method{pdf}{prior}(x, y, ...)

\method{quant}{prior}(x, p, ...)

\method{mcdf}{prior}(x, q, ...)

\method{mccdf}{prior}(x, q, ...)

\method{mlpdf}{prior}(x, y, ...)

\method{mpdf}{prior}(x, y, ...)

\method{mquant}{prior}(x, p, ...)
}
\arguments{
\item{x}{prior distribution}

\item{n}{number of observations}

\item{...}{unused arguments}

\item{q}{vector or matrix of quantiles}

\item{y}{vector of observations}

\item{p}{vector of probabilities}
}
\value{
\code{pdf} (\code{mpdf}) and \code{lpdf} (\code{mlpdf}) give
the (marginal) density and the log of (marginal) density,
\code{cdf} (\code{mcdf}) and \code{ccdf} (\code{mccdf}) give the
(marginal) distribution and the complement of (marginal) distribution function,
\code{quant} (\code{mquant}) give the (marginal) quantile function,
and \code{rng} generates random deviates for an object of class 'prior'.
}
\description{
Density (pdf / lpdf), distribution
function (cdf / ccdf), quantile function (quant),
random generation (rng), mean, standard deviation (sd),
and marginal variants of the functions (mpdf, mlpf, mcdf,
mccdf, mquant) for prior distributions.
}
\examples{
# create a standard normal prior distribution
p1 <- prior(distribution = "normal", parameters = list(mean = 1, sd = 1))

# generate a random sample from the prior
rng(p1, 10)

# compute cumulative density function
cdf(p1, 0)

# obtain quantile
quant(p1, .5)

# compute probability density
pdf(p1, c(0, 1, 2))

}
