% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{fitted.GAM}
\alias{fitted.GAM}
\title{Extract fitted values from a fable model}
\usage{
\method{fitted}{GAM}(object, ...)
}
\arguments{
\item{object}{A model for which forecasts are required.}

\item{...}{Other arguments passed to methods}
}
\value{
A vector of fitted values.
}
\description{
Extracts the fitted values.
}
\examples{
library(fabletools)
as_tsibble(USAccDeaths) \%>\%
  model(gam = GAM(log(value) ~ trend() + season())) \%>\%
  fitted()
}
