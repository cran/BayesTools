% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/JAGS-fit.R
\name{JAGS_get_inits}
\alias{JAGS_get_inits}
\title{Create initial values for 'JAGS' model}
\usage{
JAGS_get_inits(prior_list, chains, seed)
}
\arguments{
\item{prior_list}{named list of prior distribution
(names correspond to the parameter names)}

\item{chains}{number of chains}

\item{seed}{seed for random number generation}
}
\value{
\code{JAGS_add_priors} returns a list of JAGS
initial values.
}
\description{
Creates initial values for priors in
a 'JAGS' model.
}
