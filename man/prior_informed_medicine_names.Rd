% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/priors-informed.R
\docType{data}
\name{prior_informed_medicine_names}
\alias{prior_informed_medicine_names}
\title{Names of medical subfields from the Cochrane database of systematic reviews}
\format{
An object of class \code{character} of length 57.
}
\usage{
prior_informed_medicine_names
}
\value{
returns a character vector with names of medical subfields from Cochrane database of
systematic reviews.
}
\description{
Contain names identifying the individual subfields from the Cochrane database
of systematic reviews. The individual elements correspond to valid \code{name} arguments
for the \code{\link[=prior_informed]{prior_informed()}} function.
}
\examples{
print(prior_informed_medicine_names)

}
\seealso{
\code{\link[=prior_informed]{prior_informed()}}
}
\keyword{datasets}
