
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fable.GAM

<!-- badges: start -->

<!-- badges: end -->

This package provides a tidy R interface to the prophet forecasting
procedure using [fable](https://github.com/tidyverts/fable). This
package makes use of the [gam
package](https://cran.r-project.org/package=gam) for R.

## Installation

You can install the development version of fable.GAM from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("MauroFalc1/fable.GAM")
```

## Example

``` r
library(fable.GAM)
library(fable)
library(tsibble)

deaths <- as_tsibble(USAccDeaths)
fit <- deaths %>%
  model(GAM = GAM(log(value) ~ trend() + season()))

fc <- fit %>% forecast(h = "2 years")
autoplot(fc, deaths)
```

<img src="man/figures/README-example-1.png" width="100%" />

See `vignettes/intro.Rmd` for a full introduction.
