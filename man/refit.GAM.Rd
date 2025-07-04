% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{refit.GAM}
\alias{refit.GAM}
\title{Refit a `GAM`}
\usage{
\method{refit}{GAM}(object, new_data, specials = NULL, ...)
}
\arguments{
\item{object}{A model for which forecasts are required.}

\item{new_data}{A tsibble containing the time points and exogenous regressors to produce forecasts for.}

\item{specials}{(passed by \code{\link[fabletools:forecast]{fabletools::forecast.mdl_df()}}).}

\item{...}{Other arguments passed to methods}
}
\description{
Applies a fitted `GAM` to a new dataset.
}
\examples{
library(fabletools)
lung_deaths_male <- as_tsibble(mdeaths)
lung_deaths_female <- as_tsibble(fdeaths)

fit <- lung_deaths_male \%>\%
  model(GAM(value ~ trend() + season()))

report(fit)

fit \%>\%
  refit(lung_deaths_female) \%>\%
  report()
}
