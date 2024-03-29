% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/JAGS-fit.R
\name{JAGS_to_monitor}
\alias{JAGS_to_monitor}
\title{Create list of monitored parameters for 'JAGS' model}
\usage{
JAGS_to_monitor(prior_list)
}
\arguments{
\item{prior_list}{named list of prior distribution
(names correspond to the parameter names)}
}
\value{
\code{JAGS_to_monitor} returns a character vector of
parameter names.
}
\description{
Creates a vector of parameter names to be
monitored in a 'JAGS' model.
}
