# REMOVE SPECIALS FROM FORMULA
remove_specials <- function(formula, specials = NULL) {
  # Extract response
  response <- deparse(formula[[2]])

  # Extract terms object
  tt <- terms(formula, specials = if (is.null(specials)) character() else specials)
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
  reformulate(keep_terms, response = response)
}

# EXTRACT SPECIALS FROM FORMULA
extract_specials <- function(formula, specials = NULL) {
  # Extract terms object with or without user-provided specials
  tt <- terms(formula, specials = if (is.null(specials)) character() else specials)
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


trendify <- function(.dat,expr=FALSE,keep_all=TRUE, knots = NULL, origin = NULL){
  if (is.null(origin)) {
    origin <- .dat[[index_var(.dat)]][[1]]
  }
  trnd <- .dat %>%
    fabletools:::fbl_trend(knots = knots,origin = origin)
  if (expr) {
    frml <- str2lang(paste0(names(trnd),collapse = " + "))
  }else{frml=NULL}
  if (keep_all) {
    trnd <- bind_cols(.dat,trnd)
  }
  return(Filter(Negate(is.null), list(data=trnd,expr=frml)))
}

seasonify <- function(.dat,expr=FALSE,keep_all=TRUE, period = NULL){
  ssn <- .dat %>%
    fabletools:::fbl_season(period = period)
  if (expr) {
    frml <- str2lang(paste0(names(ssn),collapse = " + "))
  }else{frml=NULL}
  if (keep_all) {
    ssn <- bind_cols(.dat,ssn)
  }
  return(Filter(Negate(is.null), list(data=ssn,expr=frml)))
}


fourierify <- function(.dat,expr=FALSE,keep_all=TRUE,period, K, origin = NULL){
  if (is.null(origin)) {
    origin <- .dat[[index_var(.dat)]][[1]]
  }
  frr <- .dat %>%
    fabletools:::fbl_fourier(period = period, K = K, origin = origin)
  if (expr) {
    frml <- str2lang(paste0(names(frr),collapse = " + "))
  }else{frml=NULL}
  if (keep_all) {
    frr <- bind_cols(.dat,frr)
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
    return(update(formula, paste("~",extras,"+.")))
  }else{
    return(update(formula, paste("~ . +", extras)))
  }
}

as_model_matrix <- function(tbl) {
  stats::model.matrix(~., data = tbl)[, -1, drop = FALSE]
}


# GET SE FROM GAM MODEL ---------------------------------------------------


gam_se <- function(object, new_data) {
  terms_x <- delete.response(terms(object))
  X_new <- model.matrix(terms_x, new_data)
  Vb <- tryCatch(vcov(object), error = function(e) NULL)

  if (is.null(Vb) || ncol(X_new) != ncol(Vb)) {
    sigma <- sqrt(mean(residuals(object)^2, na.rm = TRUE))
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
                         "probit" = dnorm(eta_hat),
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
  # Mirror the logic inside forecast.GAM() --------------------------------
  dtt       <- dplyr::as_tibble(new_data)
  .xreg     <- specials$xreg[[1]]      |> dplyr::as_tibble()
  .trend    <- specials$trend[[1]]
  .season   <- specials$season[[1]]    |> dplyr::as_tibble()
  .fourier  <- specials$fourier[[1]]   |> dplyr::as_tibble()

  if (!is.null(.trend)) {
    dtt <- dplyr::bind_cols(dtt, .trend$data)
  }
  specialstibble <- NULL
  if (NROW(.season)  > 0) specialstibble <- dplyr::bind_cols(specialstibble, .season)
  if (NROW(.fourier) > 0) specialstibble <- dplyr::bind_cols(specialstibble, .fourier)

  dplyr::bind_cols(dtt, specialstibble)
}

build_gam_vars <- function(data,fml, specials) {
  # Mirror the logic inside train_GAM() --------------------------------
  dtt       <- dplyr::as_tibble(data)
  .xreg     <- specials$xreg[[1]]      |> dplyr::as_tibble()
  .trend    <- specials$trend[[1]]
  .season   <- specials$season[[1]]    |> dplyr::as_tibble()
  .fourier  <- specials$fourier[[1]]   |> dplyr::as_tibble()

  if (!is.null(.trend)) {
    dtt <- dplyr::bind_cols(dtt, .trend$data)
    fml <- add_specials(fml,.trend$str)
  }
  specialstibble <- NULL
  if (NROW(.season)  > 0) specialstibble <- dplyr::bind_cols(specialstibble, .season)
  if (NROW(.fourier) > 0) specialstibble <- dplyr::bind_cols(specialstibble, .fourier)

  gam_data <- dplyr::bind_cols(dtt,specialstibble)
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
