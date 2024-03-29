% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/JAGS-marglik.R
\name{JAGS_bridgesampling_posterior}
\alias{JAGS_bridgesampling_posterior}
\title{Prepare 'JAGS' posterior for 'bridgesampling'}
\usage{
JAGS_bridgesampling_posterior(
  posterior,
  prior_list,
  add_parameters = NULL,
  add_bounds = NULL
)
}
\arguments{
\item{posterior}{matrix of mcmc samples from the posterior
distribution}

\item{prior_list}{named list of prior distribution
(names correspond to the parameter names) of parameters not specified within the
\code{formula_list}}

\item{add_parameters}{vector of additional parameter names that should be used
in bridgesampling but were not specified in the \code{prior_list}}

\item{add_bounds}{list with two name vectors (\code{"lb"} and \code{"up"})
containing lower and upper bounds of the additional parameters that were not
specified in the \code{prior_list}}
}
\value{
\code{JAGS_bridgesampling_posterior} returns a matrix of
posterior samples with 'lb' and 'ub' attributes carrying the
lower and upper boundaries.
}
\description{
prepares posterior distribution for 'bridgesampling'
by removing unnecessary parameters and attaching lower and upper
bounds of parameters based on a list of prior distributions.
}
