# REMOVE SPECIALS FROM FORMULA
remove_specials <- function(formula, specials = NULL) {
  # Extract response
  response <- deparse(formula[[2]])

  # Extract terms object
  tt <- stats::terms(formula, specials = if (is.null(specials)) character() else specials)
  rhs_terms <- attr(tt, "term.labels")

  # Detect specials if not provided
  if (is.null(specials)) {
    specials_found <- names(attr(tt, "specials"))
  } else {
    specials_found <- specials
  }

  # Regex to detect special function calls
  pattern <- paste0("^(", paste0(specials_found, collapse = "|"), ")\\(")

  # Filter out terms matching any special
  keep_terms <- rhs_terms[!grepl(pattern, rhs_terms)]

  # Rebuild the cleaned formula
  stats::reformulate(keep_terms, response = response)
}

# EXTRACT SPECIALS FROM FORMULA
extract_specials <- function(formula, specials = NULL) {
  # Extract terms object with or without user-provided specials
  tt <- stats::terms(formula, specials = if (is.null(specials)) character() else specials)
  rhs_terms <- attr(tt, "term.labels")

  # Detect specials automatically if not given
  if (is.null(specials)) {
    specials_found <- names(attr(tt, "specials"))
  } else {
    specials_found <- specials
  }

  if (length(specials_found) == 0) return(character(0))  # nothing to extract

  # Regex to match function calls like trend(), season("x"), etc.
  pattern <- paste0("^(", paste0(specials_found, collapse = "|"), ")\\(")

  # Extract matching terms
  rhs_terms[grepl(pattern, rhs_terms)]
}

fbl_trend   <- getFromNamespace("fbl_trend",   "fabletools")
fbl_season  <- getFromNamespace("fbl_season",  "fabletools")
fbl_fourier <- getFromNamespace("fbl_fourier", "fabletools")

trendify <- function(.dat,expr=FALSE,keep_all=TRUE, knots = NULL, origin = NULL){
  if (is.null(origin)) {
    origin <- .dat[[index_var(.dat)]][[1]]
  }
  trnd <- .dat %>%
    fbl_trend(knots = knots,origin = origin)
  if (expr) {
    frml <- str2lang(paste0(names(trnd),collapse = " + "))
  }else{frml=NULL}
  if (keep_all) {
    trnd <- dplyr::bind_cols(.dat,trnd)
  }
  return(Filter(Negate(is.null), list(data=trnd,expr=frml)))
}

seasonify <- function(.dat,expr=FALSE,keep_all=TRUE, period = NULL){
  ssn <- .dat %>%
    fbl_season(period = period)
  if (expr) {
    frml <- str2lang(paste0(names(ssn),collapse = " + "))
  }else{frml=NULL}
  if (keep_all) {
    ssn <- dplyr::bind_cols(.dat,ssn)
  }
  return(Filter(Negate(is.null), list(data=ssn,expr=frml)))
}


fourierify <- function(.dat,expr=FALSE,keep_all=TRUE,period, K, origin = NULL){
  if (is.null(origin)) {
    origin <- .dat[[index_var(.dat)]][[1]]
  }
  frr <- .dat %>%
    fbl_fourier(period = period, K = K, origin = origin)
  if (expr) {
    frml <- str2lang(paste0(names(frr),collapse = " + "))
  }else{frml=NULL}
  if (keep_all) {
    frr <- dplyr::bind_cols(.dat,frr)
  }
  return(Filter(Negate(is.null), list(data=frr,expr=frml)))
}


add_specials <- function(formula, data, left=TRUE) {
  # Pick the best available set of names
  vars <- names(data) %||% colnames(data) %||% as.character(data)

  # Nothing to add? Return the original formula unchanged
  if (length(vars) == 0L)
    return(formula)

  # Build a single RHS string like "var1 + var2 + var3"
  extras <- paste0(vars, collapse = " + ")

  # Update the formula:   lhs ~ old_rhs + var1 + var2 + ...
  if (left) {
    return(stats::update(formula, paste("~",extras,"+.")))
  }else{
    return(stats::update(formula, paste("~ . +", extras)))
  }
}

as_model_matrix <- function(tbl) {
  stats::model.matrix(~., data = tbl)[, -1, drop = FALSE]
}


# GET SE FROM GAM MODEL ---------------------------------------------------


gam_se <- function(object, new_data) {
  terms_x <- stats::delete.response(stats::terms(object))
  X_new <- stats::model.matrix(terms_x, new_data)
  Vb <- tryCatch(vcov(object), error = function(e) NULL)

  if (is.null(Vb) || ncol(X_new) != ncol(Vb)) {
    sigma <- sqrt(mean(stats::residuals(object)^2, na.rm = TRUE))
    return(rep(sigma, NROW(new_data)))  # Fix: return vector
  }

  # Linear predictor SE (on link scale)
  se_link <- sqrt(rowSums((X_new %*% Vb) * X_new))

  # Apply delta method if non-identity link
  link_fun <- object$family$link
  if (link_fun != "identity") {
    eta_hat <- as.numeric(X_new %*% coef(object))  # linear predictor
    d_inv_link <- switch(link_fun,
                         "log" = exp(eta_hat),
                         "logit" = exp(eta_hat) / (1 + exp(eta_hat))^2,
                         "probit" = stats::dnorm(eta_hat),
                         "inverse" = -1 / eta_hat^2,
                         stop(paste("Link function", link_fun, "not supported"))
    )
    se_response <- se_link * abs(d_inv_link)
    return(se_response)
  }

  # Identity case
  return(se_link)
}


# Utilities --------------------------------------------------------------

build_gam_data <- function(new_data, specials) {
  # Usable within forecast, generate -------------------------------------
  dtt       <- dplyr::as_tibble(new_data)
  .xreg     <- dplyr::bind_cols(specials$xreg)
  .trend    <- specials$trend[[1]]
  .season <- dplyr::bind_cols(specials$season)
  .fourier <- dplyr::bind_cols(specials$fourier)
  if (!is.null(.trend)) {
    dtt <- dplyr::bind_cols(dtt, .trend$data)
  }
  specialstibble <- NULL
  if (NROW(.season)  > 0) specialstibble <- dplyr::bind_cols(specialstibble, .season)
  if (NROW(.fourier) > 0) specialstibble <- dplyr::bind_cols(specialstibble, .fourier)
  # if (NROW(.xreg)    > 0) specialstibble <- dplyr::bind_cols(specialstibble, .xreg)

  dplyr::bind_cols(dtt, specialstibble)
}

build_gam_data2 <- function(new_data, specials) {
  # Usable within interpolate, refit --------------------------------------
  dtt       <- dplyr::as_tibble(new_data)
  .xreg     <- dplyr::bind_cols(specials$xreg)
  .trend    <- specials$trend[[1]]
  .season <- dplyr::bind_cols(specials$season)
  .fourier <- dplyr::bind_cols(specials$fourier)
  if (!is.null(.trend)) {
    dtt <- dplyr::bind_cols(dtt, .trend$data)
  }
  specialstibble <- NULL
  if (NROW(.season)  > 0) specialstibble <- dplyr::bind_cols(specialstibble, .season)
  if (NROW(.fourier) > 0) specialstibble <- dplyr::bind_cols(specialstibble, .fourier)
  if (NROW(.xreg)    > 0) specialstibble <- dplyr::bind_cols(specialstibble, .xreg)

  dplyr::bind_cols(dtt, specialstibble)
}


build_gam_vars <- function(data,fml, specials) {
  # Mirror the logic inside train_GAM() --------------------------------
  dtt       <- dplyr::as_tibble(data)
  .xreg     <- dplyr::bind_cols(specials$xreg)
  .trend    <- specials$trend[[1]]
  .season <- dplyr::bind_cols(specials$season)
  .fourier <- dplyr::bind_cols(specials$fourier)
  if (!is.null(.trend)) {
    dtt <- dplyr::bind_cols(dtt, .trend$data)
    fml <- add_specials(fml,.trend$str)
  }
  specialstibble <- NULL
  if (NROW(.season)  > 0) specialstibble <- dplyr::bind_cols(specialstibble, .season)
  if (NROW(.fourier) > 0) specialstibble <- dplyr::bind_cols(specialstibble, .fourier)
  if (NROW(.xreg) > 0) {
    gam_data <- dplyr::bind_cols(dtt,specialstibble,.xreg)
  }else{
    gam_data <- dplyr::bind_cols(dtt,specialstibble)
  }
  .form <- remove_specials(fml,names(fabletools::common_xregs))
  gam_formula <- add_specials(.form,specialstibble)
  return(list(gam_data=gam_data,gam_formula=gam_formula))
}


#' @importFrom statmod rinvgauss
gam_draw_innov <- function(x,gam_data,mu=NULL,bootstrap=FALSE) {
  fam <- x$family$family
  if(is.null(mu)) mu <- stats::predict(x, newdata = gam_data, type = "response")
  se_fit  <- gam_se(x, gam_data)
  res <- x$residuals
  resvar <- mean(x$residuals^2, na.rm = TRUE)
  # res <- stats::residuals(x, type = "response")
  # resvar <- mean(res^2, na.rm = TRUE)

  if (bootstrap) {
    sample(stats::na.omit(res) - mean(res, na.rm = TRUE),
           size = length(mu), replace = TRUE)
  } else {
    switch(fam,
           "gaussian" =  stats::rnorm(length(mu), 0, sqrt(se_fit^2 + resvar)),
           "quasi"    =  stats::rnorm(length(mu), 0, sqrt(se_fit^2 + resvar)),
           "poisson"  =  stats::rpois(length(mu), lambda = pmax(mu, 0)) - mu,
           "binomial" = {
             size <- if (!is.null(x$prior.weights)) x$prior.weights else 1
             prob <- pmin(pmax(mu / size, 0), 1)
             stats::rbinom(length(mu), size, prob) - mu
           },
           "Gamma" = {
             phi   <- tryCatch(summary(x)$dispersion, error = function(e) 1)
             shape <- 1 / phi
             scale <- mu / shape
             stats::rgamma(length(mu), shape, scale) - mu
           },
           "inverse.gaussian" = {
             if (!requireNamespace("statmod", quietly = TRUE))
               stop("`statmod` is needed for inverse-Gaussian draws.")
             phi <- tryCatch(summary(x)$dispersion, error = function(e) 1)
             lambda <- 1 / phi
             statmod::rinvgauss(length(mu), mean = mu, shape = lambda) - mu
           },
           "Negative Binomial" = {
             size <- tryCatch(x$family$getTheta(), error = function(e) 1)
             stats::rnbinom(length(mu), mu = mu, size = size) - mu
           },
           # Fallback
           stats::rnorm(length(mu), 0, sqrt(se_fit^2 + resvar))
    )
  }
}


#' Extract exact frequencies for common seasonal periods
#'
#' This helper is derived from [`fabletools::common_periods()`]
#' but without the approximation of weekly data.
#'
#'
#' @param x A tsibble.
#'
#' @return A named vector of frequencies appropriate for the provided data.
#' @export
exact_periods <- function(x){
  x <- tsibble::interval(x)
  if (inherits(x, "vctrs_vctr")) {
    x <- vctrs::vec_data(x)
  }
  freq_sec <- c(year = 31557600, week = 604800, day = 86400,
                hour = 3600, minute = 60, second = 1, millisecond = 0.001,
                microsecond = 1e-06, nanosecond = 1e-09)
  nm <- names(x)[x != 0]
  if (is_empty(x))
    return(NULL)
  switch(paste(nm, collapse = ""), unit = c(none = 1), year = c(year = 1),
         quarter = c(year = 4/x[["quarter"]]), month = c(year = 12/x[["month"]]),
         week = c(year = 365.25/7/x[["week"]]), day = c(year = 365.25,
                                                        week = 7)/x[["day"]], with(list(secs = freq_sec/sum(as.numeric(x) *
                                                                                                              freq_sec[nm])), secs[secs > 1]))
}



#' Construct a GAM formula from a tsibble
#'
#' This helper inspects the series index using [`fabletools::common_periods()`]
#' to include appropriate `trend()`, `season()` or `fourier()` specials.
#'
#' `season()` terms are used for short integer periods (\eqn{\le threshold}), while
#' longer or non integer periods are expressed with `fourier()` using
#' \eqn{K = \lfloor period/K_d \rfloor`}.
#'
#' @param data A tsibble containing a single measured variable.
#' @param response The name of the response variable.
#' @param box_cox Whether to apply a box_cox transformation to the response.
#' @param xregs A character vector containing regressor variables' names.
#' @param trend Does the model include trend?
#' @param season Does the model include seasonalities?
#' @param threshold Maximum seasonal periods that triggers `season()` instead of `fourier()`
#' @param K_d Denominator of the maximum number of trigonometric terms \eqn{K = \lfloor period/K_d \rfloor`}
#'
#' @return A model formula suitable for [`GAM()`].
#' @export
auto_gam_formula <- function(data,response,box_cox = FALSE,xregs = NULL,
                             trend = TRUE, season = TRUE,
                             threshold = 24,K_d = 2) {
  if (missing(response)) {
    response <- tsibble::measured_vars(data)[[1]]
  }
  if (length(response) != 1) {
    stop("`auto_gam_formula()` expects a tsibble with a single response variable")
  }
  # idx  <- tsibble::index(data)
  if (box_cox) {
    lambda_guerrero <- data %>% fabletools::features(!!str2lang(response),guerrero) %>% pull() %>% pmax(.,0)
    response <- paste0("box_cox(x = ",response,", lambda = ",lambda_guerrero,")")
  }
  periods <- exact_periods(data)
  ts_len <- data %>% dplyr::count(!!!tsibble::key(.)) %>% dplyr::pull() %>% min()
  if (trend) {
    rhs_terms <- c("trend()")
  }else{
    rhs_terms <- NULL
  }
  if (season) {
    for (p in periods) {
      p_num <- suppressWarnings(as.numeric(p))
      if (ts_len > 2*p_num + 1 ) {
        if (!is.na(p_num) && abs(p_num - round(p_num)) < sqrt(.Machine$double.eps) && p_num <= threshold) {
          rhs_terms <- c(rhs_terms, paste0("season(period =", p, ")"))
        } else {
          K <- ifelse(is.na(p_num), 1L, max(1L, floor(p_num/K_d)))
          rhs_terms <- c(rhs_terms, paste0("fourier(period =", p, ", K=", K, ")"))
        }
      }
    }
  }
  rhs <- paste(c(rhs_terms,xregs), collapse = " + ")
  stats::as.formula(paste(response, "~", rhs))
}
