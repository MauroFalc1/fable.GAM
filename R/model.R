#' @docType package
#' @keywords package
"_PACKAGE"

globalVariables(c("self", "origin"))

#' Training function of GAM model
#'
#' @description
#' `train_gam()` is the internal engine that powers the user‑facing [GAM()]
#' model specification. It prepares the data, translates fable specials into
#' smooth terms, and calls [gam::gam()] from the `gam` package to estimate a
#' univariate additive model.
#'
#' @usage train_gam(.data, specials, ...)
#'
#' @param .data tsibble containing the transformed response
#' @param specials list of specials defined for GAM
#' @param ... other variables to be passed to `gam::gam()`
#'
#' @returns a model object of class GAM
#' @seealso
#'  [gam::gam()],
#' <https://www.statlearning.com/>
#' @keywords internal
#' @importFrom stats predict
train_gam <- function(.data, specials, ...) {
  dtt <- self$data
  fml <- self$formula
  idx <- dtt %>% dplyr::pull(tsibble::index(dtt))
  av <- all.vars(stats::terms.formula(fml))
  mv <- tsibble::measured_vars(.data)
  if (length(mv) > 1) stop("GAM() is a univariate model.")
  y <- .data[[mv]]
  dtt <- dplyr::bind_cols(
    tibble("t_response" = y),
    dplyr::as_tibble(dtt) %>% dplyr::select(dplyr::all_of(av[-1]))
  )
  fml <- rlang::new_formula(lhs = quote(t_response), rhs = fml[[3]])
  obj <- build_gam_vars(data = dtt, fml = fml, specials = specials)
  gam_data <- obj$gam_data
  gam_formula <- obj$gam_formula
  fit <- gam::gam(gam_formula, data = gam_data, ...)
  fit$index <- idx
  structure(fit, class = c("GAM", class(fit)))
}

#' Special functions for the GAM model
#'
#' These special functions provide interfaces to more complicated functions
#' within the model formulae interface.
#'
#' @section Specials:
#'
#' \subsection{trend}{
#' The `trend` special includes common linear and smoothened trend regressors in the model.
#' It also supports piecewise linear trend via the `knots` argument.
#' \preformatted{
#' trend(knots = NULL, origin = NULL, linear=FALSE, ...)
#' }
#'
#' \tabular{ll}{
#'   `knots`    \tab A vector of times (same class as the data's time index)
#'   identifying the position of knots for a piecewise linear trend.\cr
#'   `origin`   \tab An optional starting time value for the trend. \cr
#'   `linear` \tab logical indicating whether using linear or smoothened trend.\cr
#'   `...` \tab other variables to be passed to `gam::s()`
#' }
#' }
#'
#' \subsection{season}{
#' The `season` special includes seasonal dummy variables in the model.
#' \preformatted{
#' season(period = NULL)
#' }
#'
#' \tabular{ll}{
#'   `period`   \tab The periodic nature of the seasonality. This can be either a number indicating the number of observations in each seasonal period, or text to indicate the duration of the seasonal window (for example, annual seasonality would be "1 year").
#' }
#' }
#'
#' \subsection{fourier}{
#' The `fourier` special includes seasonal fourier terms in the model. The maximum order of the fourier terms must be specified using `K`.
#' \preformatted{
#' fourier(period = NULL, K, origin = NULL)
#' }
#'
#' \tabular{ll}{
#'   `period`   \tab The periodic nature of the seasonality. This can be either a number indicating the number of observations in each seasonal period, or text to indicate the duration of the seasonal window (for example, annual seasonality would be "1 year"). \cr
#'   `K`        \tab The maximum order of the fourier terms.\cr
#'   `origin`   \tab An optional time value to act as the starting time for the fourier series.
#' }
#' }
#'
#' @format NULL
#' @keywords internal
#' @rdname specials_gam
#' @seealso
#' * [fabletools::new_specials()],
#' * <https://fabletools.tidyverts.org/articles/extension_models.html>
specials_gam <- fabletools::new_specials(
  trend = function(knots = NULL, origin = NULL, linear = FALSE, ...) {
    if (is.null(origin)) {
      if (is.null(self$origin)) {
        self$origin <- self$data[[index_var(self$data)]][[1]]
      }
      origin <- self$origin
    }
    out.dat <- fabletools:::fbl_trend(self$data, knots, origin)
    if (linear) {
      out.str <- paste0(names(out.dat), collapse = " + ")
      out.expr <- rlang::expr(!!paste0(names(out.dat), collapse = " + "))
    } else {
      vars <- names(out.dat)
      calls <- lapply(vars, function(v) rlang::expr(s(!!rlang::sym(v), !!!list(...))))
      # Build a flat addition expression instead of nested reduction
      if (length(calls) == 1) {
        out.expr <- calls[[1]]
      } else {
        out.expr <- rlang::call2("+", !!!calls)
      }
      out.str <- deparse(out.expr)
    }
    return(list(data = out.dat, expr = out.expr, str = out.str))
  },
  season = function(period = NULL) {
    out <- as_model_matrix(fabletools:::fbl_season(
      self$data,
      period
    ))
    stats::model.matrix(~., data = as.data.frame(out))[, -1, drop = FALSE]
  },
  fourier = function(period = NULL, K, origin = NULL) {
    if (is.null(origin)) {
      if (is.null(self$origin)) {
        self$origin <- self$data[[index_var(self$data)]][[1]]
      }
      origin <- self$origin
    }
    as.matrix(fabletools:::fbl_fourier(
      self$data, period, K,
      origin
    ))
  },
  xreg = fabletools::special_xreg(default_intercept = FALSE),
  .required_specials = NULL # ,
  # .xreg_specials = names(fabletools::common_xregs)
)

#' Fit a Generalized Additive Model (GAM) with time series components
#'
#' The model formula will be handled using [`gam::gam()`], and so
#' the same approach to include smoothing functions in [`gam::gam()`] applies when
#' specifying the `formula`. In addition, it is possible to
#' include [`specials_gam`] in the model formula, such as `trend()`, `season()`,
#' and `fourier()`.
#'
#' @param formula Model specification.
#' @param ... Any other orguments being passed to [`gam::gam()`]
#'
#' @section Specials:
#'
#' \subsection{trend}{
#' The `trend` special includes common linear and smoothened trend regressors in the model.
#' It also supports piecewise linear trend via the `knots` argument.
#' \preformatted{
#' trend(knots = NULL, origin = NULL, linear=FALSE, ...)
#' }
#'
#' \tabular{ll}{
#'   `knots`    \tab A vector of times (same class as the data's time index)
#'   identifying the position of knots for a piecewise linear trend.\cr
#'   `origin`   \tab An optional starting time value for the trend. \cr
#'   `linear` \tab logical indicating whether using linear or smoothened trend.\cr
#'   `...` \tab other variables to be passed to `gam::s()`
#' }
#' }
#'
#' \subsection{season}{
#' The `season` special includes seasonal dummy variables in the model.
#' \preformatted{
#' season(period = NULL)
#' }
#'
#' \tabular{ll}{
#'   `period`   \tab The periodic nature of the seasonality. This can be either a number indicating the number of observations in each seasonal period, or text to indicate the duration of the seasonal window (for example, annual seasonality would be "1 year").
#' }
#' }
#'
#' \subsection{fourier}{
#' The `fourier` special includes seasonal fourier terms in the model. The maximum order of the fourier terms must be specified using `K`.
#' \preformatted{
#' fourier(period = NULL, K, origin = NULL)
#' }
#'
#' \tabular{ll}{
#'   `period`   \tab The periodic nature of the seasonality. This can be either a number indicating the number of observations in each seasonal period, or text to indicate the duration of the seasonal window (for example, annual seasonality would be "1 year"). \cr
#'   `K`        \tab The maximum order of the fourier terms.\cr
#'   `origin`   \tab An optional time value to act as the starting time for the fourier series.
#' }
#' }
#'
#' \subsection{xreg}{
#' Exogenous regressors can be included in a GAM model without explicitly using the `xreg()` special.
#' \preformatted{
#' xreg(...)
#' }
#'
#' \tabular{ll}{
#'   `...`      \tab Bare expressions for the exogenous regressors (such as `s(x)`)
#' }
#' }
#'
#' @return A model specification.
#'
#' @seealso
#' [`gam::gam()`],
#' [An Introduction to Statistical Learning, Moving Beyond Linearity (chapter 7)](https://www.statlearning.com/)
#'
#' @examples
#' library(fabletools)
#' as_tsibble(USAccDeaths) %>%
#'   model(gam = GAM(log(value) ~ trend() + season()))
#'
#' library(tsibbledata)
#' olympic_running %>%
#'   model(GAM(Time ~ trend())) %>%
#'   interpolate(olympic_running)
#' @export
GAM <- function(formula, ...) {
  gam_model <- fabletools::new_model_class("GAM",
    train = train_gam,
    specials = specials_gam,
    origin = NULL
  )
  fabletools::new_model_definition(gam_model, !!enquo(formula), ...)
}

#' Produce forecasts from a GAM model
#'
#' Produces forecasts from a trained model as in the `fable` package.
#'
#' @inheritParams fable::forecast.TSLM
#' @param approx_normal  Should Gaussian forecasts be approximated
#'   by a Normal (default) rather than a Student‑t?
#'
#' @returns A list of forecasts.
#' @examples
#' library(fabletools)
#' USAccDeaths %>%
#'   as_tsibble() %>%
#'   model(gam = GAM(log(value) ~ trend() + season())) %>%
#'   forecast()
#' @export
forecast.GAM <- function(object,
                         new_data,
                         specials = NULL,
                         bootstrap = FALSE,
                         approx_normal = TRUE,
                         times = 5e3,
                         ...) {
  # 1 ────────────────────────────────────────────────────────────────
  # Re‑create the design matrix for the *new* data
  # ──────────────────────────────────────────────────────────────────
  gam_data <- build_gam_data(new_data, specials)

  # 2 ────────────────────────────────────────────────────────────────
  # Point forecasts & standard errors on the *response* scale
  # ──────────────────────────────────────────────────────────────────
  pr <- predict(
    object,
    newdata = gam_data,
    se.fit  = FALSE,
    type    = "response"
  )

  fc_mean <- pr
  se_fit <- gam_se(object, gam_data)

  # Residual variance (for Gaussian‑type models)
  resvar <- mean(object$residuals^2, na.rm = TRUE)

  # 3 ────────────────────────────────────────────────────────────────
  # Build forecast distribution
  # ──────────────────────────────────────────────────────────────────
  if (bootstrap) {
    # ─ Bootstrap / simulated predictive distribution ────────────────
    res <- stats::na.omit(object$residuals) - mean(object$residuals, na.rm = TRUE)

    sims <- replicate(
      times,
      fc_mean + sample(res, size = length(fc_mean), replace = TRUE),
      simplify = FALSE
    )
    return(distributional::dist_sample(sims))
  }

  # ─ Analytic distribution informed by the model family ──────────────
  fam <- object$family$family

  switch(fam,
    # ──────────────────── Gaussian & quasi families ──────────────────
    "gaussian" = {
      sd_tot <- sqrt(se_fit^2 + resvar)
      if (approx_normal) {
        distributional::dist_normal(fc_mean, sd_tot)
      } else {
        df <- if (!is.null(object$df.residual)) {
          object$df.residual
        } else {
          length(object$residuals) - length(object$coefficients)
        }
        distributional::dist_student_t(df, fc_mean, sd_tot)
      }
    },
    "quasi" = {
      sd_tot <- sqrt(se_fit^2 + resvar)
      distributional::dist_normal(fc_mean, sd_tot)
    },

    # ─────────────────────────── Poisson ─────────────────────────────
    "poisson" = {
      distributional::dist_poisson(fc_mean)
    },

    # ────────────────────────── Binomial ─────────────────────────────
    "binomial" = {
      size <- if (!is.null(object$prior.weights)) object$prior.weights else 1
      prob <- pmin(pmax(fc_mean / size, 0), 1)
      distributional::dist_binomial(size, prob)
    },

    # ─────────────────────────── Gamma ───────────────────────────────
    "Gamma" = {
      phi <- tryCatch(summary(object)$dispersion, error = function(e) 1)
      shape <- 1 / phi
      scale <- fc_mean / shape
      distributional::dist_gamma(shape, scale)
    },

    # ─────────────────────── Inverse‑Gaussian ────────────────────────
    "inverse.gaussian" = {
      phi <- tryCatch(summary(object)$dispersion, error = function(e) 1)
      lambda <- 1 / phi
      distributional::dist_inverse_gaussian(fc_mean, lambda)
    },

    # ───────────────────── Negative Binomial ─────────────────────────
    "Negative Binomial" = {
      theta <- tryCatch(object$family$theta, error = function(e) 1)
      distributional::dist_negative_binomial(mu = fc_mean, size = theta)
    },

    # ─────────────────────── Default fallback ───────────────────────
    {
      sd_tot <- sqrt(se_fit^2 + resvar)
      distributional::dist_normal(fc_mean, sd_tot)
    }
  )
}


## -------------------------------------------------------------------------
## fitted() -----------------------------------------------------------------

#' @inherit fable::fitted.ARIMA
#'
#' @examples
#' library(fabletools)
#' as_tsibble(USAccDeaths) %>%
#'   model(gam = GAM(log(value) ~ trend() + season())) %>%
#'   fitted()
#' @export
fitted.GAM <- function(object, ...) {
  object$fitted.values
}

## -------------------------------------------------------------------------
## residuals() --------------------------------------------------------------

#' @inherit fable::residuals.ARIMA
#'
#' @examples
#' library(fabletools)
#' as_tsibble(USAccDeaths) %>%
#'   model(gam = GAM(log(value) ~ trend() + season())) %>%
#'   residuals()
#' @export
residuals.GAM <- function(object,
                          type = c(
                            "innovation", "response", "deviance",
                            "pearson", "working", "partial"
                          ),
                          ...) {
  type <- match.arg(type)
  object$residuals # all flavours are identical for a GAM
}
# register an alias for the underlying "Gam" class too
#' @export
residuals.Gam <- residuals.GAM


#' Glance at a GAM model
#'
#' Returns a one-row tibble summarizing the model's fit.
#'
#' @param x A `GAM` model object.
#' @param ... Additional arguments (not currently used).
#'
#' @return A one-row `tibble` with columns such as `r.squared`, `adj.r.squared`,
#'   `deviance`, `df.residual`, `log_lik`, `AIC`, `BIC`, `dispersion`, and `nobs`.
#'
#' @examples
#' library(fabletools)
#' as_tsibble(USAccDeaths) %>%
#'   model(gam = GAM(log(value) ~ trend() + season())) %>%
#'   glance()
#' @export
glance.GAM <- function(x, ...) {
  sumgam <- summary(x)

  # Extract dispersion. For some families, summary(x)$dispersion might be NULL
  # or not the primary way it's reported (e.g., Poisson, Binomial have fixed dispersion of 1).
  # gam::summary.gam stores it, and it's used for quasi families.
  dispersion_val <- sumgam$dispersion %||% NA_real_
  if (x$family$family %in% c("poisson", "binomial") && is.na(dispersion_val)) {
    dispersion_val <- 1.0 # Theoretical dispersion
  }


  tibble::tibble(
    # r.squared      = sumgam$r.sq %||% NA_real_,
    # adj.r.squared  = sumgam$adj.r.sq %||% NA_real_,
    deviance       = stats::deviance(x),
    df.residual    = stats::df.residual(x),
    log_lik        = as.numeric(stats::logLik(x)), # Ensure it's a plain number
    AIC            = stats::AIC(x),
    BIC            = stats::BIC(x),
    dispersion     = dispersion_val,
    nobs           = stats::nobs(x) %||% length(x$residuals) # Number of observations
  )
}


#' Summarise a GAM model fit (Hastie & Tibshirani `gam` package)
#'
#' Provides a console report formatted similarly to other `fable` model
#' `report()` methods, but tailored for objects fitted with the **gam** package.
#' In particular it makes use of the `anova` and `parametric.anova` tables
#' returned by `summary.gam()`.
#'
#' @param object A fitted GAM object (class `GAM` created by `train_gam()` or
#'   equivalent helper inside this package).
#' @param digits Minimum number of significant digits to print.
#' @param ...    Further arguments passed for S3 completeness (currently
#'   ignored).
#'
#' @return (Invisibly) returns `object` so the call can be piped/assigned.
#' @examples
#' library(fabletools)
#' as_tsibble(USAccDeaths) %>%
#'   model(gam = GAM(log(value) ~ trend() + season())) %>%
#'   report()
#' @export
report.GAM <- function(object, digits = max(3L, getOption("digits") - 3L), ...) {
  # Protect against the wrong class -------------------------------------------------
  if (!inherits(object, "GAM")) {
    stop("`report.GAM()` expects an object of class 'GAM'.")
  }

  sm <- summary(object)

  cli_line <- function(char = "-", n = 60) cat(strrep(char, n), "\n", sep = "")

  # Header -------------------------------------------------------------------------
  cat("GAM Model Report", "\n", sep = "")
  cli_line("=")
  if (!is.null(object$family)) {
    fam <- object$family
    cat("Family:\t", fam$family, " (", fam$link, " link)", "\n", sep = "")
  }
  if (!is.null(object$method)) {
    cat("Fitting method:\t", object$method, "\n", sep = "")
  }
  if (!is.null(sm$rank)) {
    cat("Rank (edf):\t", sm$rank, "\n", sep = "")
  }

  # Parametric component -----------------------------------------------------------
  if (!is.null(sm$parametric.anova)) {
    cat("\nParametric terms (Wald tests)\n")
    cli_line()
    stats::printCoefmat(sm$parametric.anova,
      digits = digits, signif.stars = TRUE,
      P.values = TRUE, na.print = "NA"
    )
  }

  # Smooth component ----------------------------------------------------------------
  if (!is.null(sm$anova)) {
    cat("\nSmooth terms (Approx. chi-square tests)\n")
    cli_line()
    stats::printCoefmat(sm$anova,
      digits = digits, signif.stars = TRUE,
      P.values = TRUE, has.Pvalue = TRUE, na.print = "NA"
    )
  }

  # Goodness of fit ------------------------------------------------------------------
  cat("\nModel statistics\n")
  cli_line()
  if (!is.null(sm$deviance)) {
    cat(sprintf("%-20s %.*f\n", "Deviance:", digits, sm$deviance))
  }
  if (!is.null(sm$r.sq)) {
    cat(sprintf("%-20s %.*f\n", "R-squared:", digits, sm$r.sq))
  }
  if (!is.null(sm$dispersion)) {
    cat(sprintf("%-20s %.*f\n", "Dispersion:", digits, sm$dispersion))
  }
  if (!is.null(object$aic)) {
    cat(sprintf("%-20s %.*f\n", "AIC:", digits, object$aic))
  }
  if (!is.null(sm$gcv.ubre)) {
    cat(sprintf("%-20s %.*f\n", "GCV/UBRE:", digits, sm$gcv.ubre))
  }
  cat(sprintf("%-20s %d\n", "Observations:", length(object$y)))

  invisible(object)
}


# ---------------------------------------------------------------------------
# model_sum() ----------------------------------------------------------------

#' @export
model_sum.GAM <- function(x, ...) {
  # sm  <- summary(x)
  # fam <- family(x)
  # sprintf("GAM(%s/%s, edf = %.1f)", fam$family, fam$link, sum(sm$edf))
  "GAM"
}

#' @export
format.GAM <- function(x, ...) {
  "Generalized Additive Model (GAM)"
}


# ---------------------------------------------------------------------------
#  tidy() --------------------------------------------------------------------
# ---------------------------------------------------------------------------

#' Tidy the coefficients (and smooth term tests) from a GAM model
#'
#' Creates a tibble containing the parametric coefficients with standard errors
#' and p-values, and (optionally) the approximate chi‑square tests for each
#' smooth term. This parallels the behaviour of `broom::tidy.glm()` but is
#' adapted for the **gam** package’s output.
#'
#' @param x A fitted `GAM` object.
#' @param include_smooth Should smooth-term tests be included? (default `TRUE`).
#' @param ... Further arguments for S3 completeness (ignored).
#'
#' @return A `tibble` with one row per term.
#' @examples
#' library(fabletools)
#' as_tsibble(USAccDeaths) %>%
#'   model(gam = GAM(log(value) ~ trend() + season())) %>%
#'   tidy()
#'
#' @export
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
#' @importFrom stats coef vcov pnorm
tidy.GAM <- function(x, include_smooth = TRUE, ...) {
  if (!inherits(x, "GAM")) stop("`tidy.GAM()` expects a <GAM> object.")

  # ── Parametric coefficients ----------------------------------------------
  beta <- stats::coef(x)
  se <- sqrt(diag(stats::vcov(x)))
  stat <- beta / se
  pval <- 2 * stats::pnorm(abs(stat), lower.tail = FALSE)

  tbl_param <- tibble::tibble(
    term      = names(beta),
    estimate  = unname(beta),
    std.error = se,
    statistic = stat,
    p.value   = pval,
    component = "parametric"
  )

  if (!include_smooth) {
    return(tbl_param)
  }

  # ── Smooth‑term chi‑square / F tests -------------------------------------
  sm <- summary(x)
  tbl_smooth <- NULL
  if (!is.null(sm$anova)) {
    an <- as.data.frame(sm$anova)
    an$term <- rownames(an)

    stat_col <- if ("Chi.sq" %in% names(an)) "Chi.sq" else if ("F" %in% names(an)) "F" else NULL
    p_col <- if ("Pr(>Chi)" %in% names(an)) "Pr(>Chi)" else if ("Pr(F)" %in% names(an)) "Pr(F)" else NULL

    tbl_smooth <- tibble::tibble(
      term      = an$term,
      df        = an$Df,
      statistic = if (!is.null(stat_col)) an[[stat_col]] else NA_real_,
      p.value   = if (!is.null(p_col)) an[[p_col]] else NA_real_,
      component = "smooth"
    )
  }

  dplyr::bind_rows(tbl_param, tbl_smooth)
}

#' Generate future sample paths from a <GAM> model
#'
#' @inheritParams fable::generate.ETS
#' @param times Integer. Number of simulated paths to produce. Default 1.
#'
#' @examples
#' library(fabletools)
#' as_tsibble(USAccDeaths) %>%
#'   model(gam = GAM(log(value) ~ trend() + season())) %>%
#'   generate()
#'
#' @export
generate.GAM <- function(x,
                         new_data,
                         specials = NULL,
                         bootstrap = FALSE,
                         times = 1,
                         ...) {
  # ── 1. Rebuild the design matrix for the future horizon ────────────────
  gam_data <- build_gam_data(new_data, specials)

  # ── 2. Mean response on the future horizon ─────────────────────────────
  mu <- stats::predict(x, newdata = gam_data, type = "response")

  # ── 3. Generate a random future paths and bind them together ───────────-
  .innov <- gam_draw_innov(x, gam_data, mu, bootstrap)
  dplyr::transmute(new_data, .sim = as.numeric(mu + .innov))
}


#' Refit a `GAM`
#'
#' Applies a fitted `GAM` to a new dataset.
#'
#' @inheritParams fable::refit.ARIMA
#'
#' @examples
#' library(fabletools)
#' lung_deaths_male <- as_tsibble(mdeaths)
#' lung_deaths_female <- as_tsibble(fdeaths)
#'
#' fit <- lung_deaths_male %>%
#'   model(GAM(value ~ trend() + season()))
#'
#' report(fit)
#'
#' fit %>%
#'   refit(lung_deaths_female) %>%
#'   report()
#' @export
refit.GAM <- function(object,
                      new_data,
                      specials = NULL,
                      ...) {
  # 0 ── Early exit ───────────────────────────────────────────────────────────
  if (is.null(new_data) || nrow(new_data) == 0) {
    warning("`new_data` is empty; returning original `object` unchanged.")
    return(object)
  }

  # 1 ── Extract original specification ───────────────────────────────────────
  gam_formula <- stats::formula(object) # stored by gam::gam
  # av <- all.vars(stats::terms.formula(gam_formula))
  mv <- tsibble::measured_vars(new_data)
  if (length(mv) > 1) stop("GAM() is a univariate model.")
  y <- new_data[[mv]]
  gam_data <- build_gam_data(dplyr::bind_cols(
    tibble("t_response" = y),
    dplyr::as_tibble(new_data)
  ), specials)
  gam_family <- object$family # also stored by gam::gam

  # 3 ── Build model frame for the new data ───────────────────────────────────

  # 4 ── Fit the GAM on the updated data ──────────────────────────────────────
  refit_obj <- gam::gam(
    formula = gam_formula,
    family  = gam_family,
    data    = gam_data,
    ...
  )

  # 5 ── Preserve the custom GAM class and any needed attributes ──────────────
  #      (add back 'specials' so generate()/forecast()/interpolate() keep working)
  refit_obj$index <- new_data %>% dplyr::pull(tsibble::index(new_data))
  structure(refit_obj, class = c("GAM", class(refit_obj)))
}

#' @inherit fable::interpolate.ARIMA
#'
#' @examples
#' library(fabletools)
#' library(tsibbledata)
#'
#' olympic_running %>%
#'   model(gam = GAM(Time ~ trend())) %>%
#'   interpolate(olympic_running)
#' @export
interpolate.GAM <- function(object, new_data, specials, ...) {
  # Get missing values
  y <- unclass(new_data)[[tsibble::measured_vars(new_data)]]
  miss_val <- which(is.na(y))
  gam_data <- build_gam_data(new_data, specials)

  # 2 ────────────────────────────────────────────────────────────────
  # Point forecasts & standard errors on the *response* scale
  # ──────────────────────────────────────────────────────────────────
  pr <- predict(
    object,
    newdata = gam_data,
    se.fit  = FALSE,
    type    = "response"
  )
  fits <- pr
  if (length(y) != length(fits)) {
    abort("Interpolation for GAM models is only supported for data used to estimate the model.")
  }

  # Update data
  y[miss_val] <- fits[miss_val]
  new_data[[tsibble::measured_vars(new_data)]] <- y
  new_data
}

#' Extract components from a GAM model
#'
#' @param object A fitted model of class GAM.
#' @param simplify Choose whether to sum-up all trend and seasonal effects, respectively.
#' @param ... Currently ignored.
#'
#' @return A [`fabletools::dable`] whose columns sum to `.response` *exactly*.
#' @examples
#' library(fabletools)
#' as_tsibble(USAccDeaths) %>%
#'   model(gam = GAM(log(value) ~ trend() + season())) %>%
#'   components()
#' @export
components.GAM <- function(object, simplify = TRUE, ...) {
  # ── Guards ---------------------------------------------------------------
  if (!inherits(object, "GAM")) {
    abort("`components.GAM()` expects an object of class 'GAM'.")
  }

  # ── 1. Core pieces ------------------------------------------------------
  resp <- stats::fitted(object) # response‑scale fitted values
  n <- length(resp)

  # Per‑term contributions (link scale) —— returned matrix has attr "constant"
  term_mat <- stats::predict(object, type = "terms", se.fit = FALSE)
  term_tbl <- tibble::as_tibble(term_mat, .name_repair = "minimal")

  if (simplify) {
    trend_term_tbl <- term_tbl %>%
      dplyr::select(dplyr::starts_with("s(trend")) %>%
      dplyr::mutate("trends" = rowSums(., na.rm = TRUE), .keep = "none")
    season_term_tbl <- term_tbl %>%
      dplyr::select(
        dplyr::starts_with(c(
          "season_", "year", "quarter",
          "month", "week", "day"
        )),
        dplyr::matches("^[SC]\\d+_\\d+$")
      ) %>%
      dplyr::mutate("seasonalities" = rowSums(., na.rm = TRUE), .keep = "none")
    term_tbl <- dplyr::bind_cols(trend_term_tbl, season_term_tbl, term_tbl %>%
                                   dplyr::select(
                                     -dplyr::starts_with(c(
                                       "s(trend", "season_", "year", "quarter",
                                       "month", "week", "day"
                                     )),
                                     -dplyr::matches("^[SC]\\d+_\\d+$")
                                   ))
  }

  # Make every column name syntactically valid **and** match what `all.vars()`
  # will pull out of the alias expression (→ important for `autoplot()`)
  original <- names(term_tbl)
  cleaned <- vapply(original, make.names, FUN.VALUE = character(1), unique = TRUE)
  names(term_tbl) <- cleaned

  # Intercept (scalar) recycled to n rows
  icpt <- attr(term_mat, "constant")
  intercept <- rep_len(icpt, n)

  # ── 2. Build a tsibble ---------------------------------------------------
  cmp <- tibble::tibble(
    .idx      = object$index,
    .response = resp,
    intercept = intercept
  ) |>
    dplyr::bind_cols(term_tbl) |>
    tsibble::as_tsibble(index = .idx)

  # ── 3. Wrap as dable -----------------------------------------------------
  alias_expr <- rlang::parse_expr( # .response ≡ sum(components)
    paste(c("intercept", cleaned), collapse = " + ")
  )

  fabletools::as_dable(
    cmp,
    resp = !!rlang::sym(".response"),
    method = "GAM",
    aliases = rlang::set_names(list(alias_expr), ".response")
  )
}
