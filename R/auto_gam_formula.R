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
#' @param box_cox_lambda A numeric value for the box_cox transformation parameter
#' @param xregs A character vector containing regressor variables' names.
#' @param trend Does the model include trend?
#' @param trend_args A list of named arguments to pass to `trend()`
#' @param season Does the model include seasonalities?
#' @param season_spec Optional manual specification of seasonalities
#' @param threshold Maximum seasonal periods that triggers `season()` instead of `fourier()`
#' @param K_d Denominator of the maximum number of trigonometric terms \eqn{K = \lfloor period/K_d \rfloor`}
#'
#' @return A model formula suitable for [`GAM()`].
#' @export
auto_gam_formula <- function(data,response,box_cox = FALSE,box_cox_lambda = "auto",
                             xregs = NULL,
                             trend = TRUE,trend_args=list(),
                             season = TRUE,season_spec=NULL,
                             threshold = 24,K_d = 2) {
  if (missing(response)) {
    response <- tsibble::measured_vars(data)[[1]]
  }
  if (length(response) != 1) {
    stop("`auto_gam_formula()` expects a tsibble with a single response variable")
  }
  # idx  <- tsibble::index(data)
  if (box_cox) {
    if (box_cox_lambda=="auto") {
      lambda_guerrero <- data %>% fabletools::features(!!str2lang(response),guerrero) %>% pull() %>% pmax(.,0)
      box_cox_lambda <- lambda_guerrero
    }
    response <- paste0("box_cox(x = ",response,", lambda = ",box_cox_lambda,")")
  }
  periods <- exact_periods(data)
  ts_len <- data %>% dplyr::count(!!!tsibble::key(.)) %>% dplyr::pull() %>% min()
  if (trend) {
    rhs_terms <- paste0("trend(",
                        paste(names(trend_args), trend_args, sep = "=", collapse = ","),
                        ")")
    # rhs_terms <- c("trend()")
  }else{
    rhs_terms <- NULL
  }
  if (season) {
    if (is.null(season_spec)) {
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
    }else{
      rhs_terms <- c(rhs_terms, season_spec)
    }
  }
  rhs <- paste(c(rhs_terms,xregs), collapse = " + ")
  if (is.null(rhs)) {
    rhs <- "1"
  }
  stats::as.formula(paste(response, "~", rhs))
}
