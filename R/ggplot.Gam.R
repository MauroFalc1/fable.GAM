#' ggplot‑based term plots for a Gam object
#'
#' A faithful re‑implementation of \code{plot.Gam()} that
#' produces one \code{ggplot} object per model term (main effects
#' of one or two predictors).  All arguments mirror those of
#' \code{plot.Gam()} so existing code can simply replace the
#' function name.
#'
#' @inheritParams gam::plot.Gam
#' @return (invisibly) a named list of \code{ggplot} objects.
#'         When \code{ask = FALSE} the plots are printed
#'         sequentially.  When \code{ask = TRUE} a menu lets you
#'         choose which term(s) to draw, exactly as in
#'         \code{plot.Gam()}.
#' @author 2025‑07‑19, adapted from Hastie’s original
#'         \code{plot.Gam()} design.  Requires packages
#'         \pkg{ggplot2}, \pkg{rlang}, \pkg{reshape2}.
#' @export
ggplot.Gam <- function(x,
                       residuals = NULL,
                       rugplot   = TRUE,
                       se        = FALSE,
                       scale     = 0,
                       ask       = FALSE,
                       terms     = labels.Gam(x),
                       ...)
{
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))

  ## ------------------------------------------------------------------
  ## 1.  Pre‑compute the plotting data (identical to plot.Gam) …… :contentReference[oaicite:0]{index=0}
  ## ------------------------------------------------------------------
  pp <- x$preplot
  if (is.null(pp))
    pp <- gam::preplot.Gam(x, terms = terms, ...)

  all.resid <- stats::resid(x)
  if (!is.null(residuals)) {
    if (length(residuals) == 1L)          # TRUE/FALSE flag
      residuals <- if (isTRUE(residuals)) all.resid else NULL
  } else {
    residuals <- NULL
  }

  ## The workhorse: build a ggplot for ONE preplot.Gam term ----------
  build_plot <- function(obj, common_ylim = NULL)
  {
    ## obj is a single preplot.Gam element  …… :contentReference[oaicite:1]{index=1}
    xvar   <- obj$x
    fit    <- obj$y
    se.fit <- obj$se.y
    labelx <- paste(obj$xlab, collapse = ":")
    labely <- obj$ylab

    ## Work out partial residuals for this term if requested
    if (!is.null(residuals)) {
      idx <- seq_along(fit)
      presid <- fit + residuals[idx]
      res_df <- data.frame(x = as.vector(xvar), y = presid)
    }

    if (is.numeric(xvar)) {                       # ------------- 1D smooth
      df <- data.frame(x = xvar, fit = fit)
      p <- ggplot2::ggplot(df, ggplot2::aes(x, fit)) +
        ggplot2::geom_line()

      if (se && !is.null(se.fit)) {               # ±2×SE ribbon
        df$upper <- fit + 2 * se.fit
        df$lower <- fit - 2 * se.fit
        p <- p + ggplot2::geom_ribbon(
          ggplot2::aes(x,ymin = df$lower, ymax = df$upper),
          alpha = .2, inherit.aes = FALSE)
      }

      if (!is.null(residuals))
        p <- p + ggplot2::geom_point(data = res_df,
                                     ggplot2::aes(x, y),
                                     alpha = .4, inherit.aes = FALSE)

      if (rugplot)
        p <- p + ggplot2::geom_rug(sides = "b", alpha = .15)

    } else if (is.factor(xvar)) {                 # ------------- factor
      df <- data.frame(x = xvar, fit = fit)
      p <- ggplot2::ggplot(df,
                           ggplot2::aes(x = xvar, y = fit)) +
        ggplot2::geom_boxplot(size=0.2)

      if (se && !is.null(se.fit))
        p <- p + ggplot2::geom_errorbar(
          ggplot2::aes(ymin = fit - 2 * se.fit,
                       ymax = fit + 2 * se.fit),
          width = 0.2,colour="grey70")

      if (!is.null(residuals))
        p <- p + ggplot2::geom_point(data = res_df,
                                     ggplot2::aes(x, y),
                                     alpha = .15, inherit.aes = FALSE)

      if (rugplot)
        p <- p + ggplot2::geom_rug(sides = "b", alpha = .15,colour="grey70")

    } else if (is.matrix(xvar) && ncol(xvar) == 2L) {  # --------- 2D smooth
      grid  <- data.frame(x1 = xvar[, 1],
                          x2 = xvar[, 2],
                          z  = fit)
      p <- ggplot2::ggplot(grid,
                           ggplot2::aes(x = x1, y = x2, z = z)) +
        ggplot2::geom_contour_filled() +
        ggplot2::labs(fill = "fit")
    } else {
      stop("Don't know how to plot term of class ",
           paste(class(xvar), collapse = "/"))
    }

    if (!is.null(common_ylim))
      p <- p + ggplot2::coord_cartesian(ylim = common_ylim)

    p + ggplot2::labs(x = labelx, y = labely)
  }

  ## ------------------------------------------------------------------
  ## 2.  Determine a common y‑scale if requested
  ## ------------------------------------------------------------------
  if (scale > 0) {
    ranges <- vapply(pp,
                     function(o) diff(range(o$y, na.rm = TRUE)),
                     numeric(1))
    big    <- max(scale, max(ranges))
    ylim   <- c(-big / 2,  big / 2)
  } else ylim <- NULL

  ## ------------------------------------------------------------------
  ## 3.  Create one ggplot per term
  ## ------------------------------------------------------------------
  plots <- lapply(pp, build_plot, common_ylim = ylim)
  names(plots) <- names(pp)

  ## ------------------------------------------------------------------
  ## 4.  Interactive or sequential display (mirrors plot.Gam)
  ## ------------------------------------------------------------------
  if (!ask) {
    # for (p in plots) print(p)
    plots
  } else {
    choices <- c(paste0("plot: ", substr(names(plots), 1L, 40L)),
                 "plot all terms", "exit")
    repeat {
      pick <- utils::menu(choices,
                          title = "Select term to draw (or 0/exit)")
      if (pick == 0 || pick == length(choices)) break
      if (pick == length(choices) - 1) {       # all terms
        # lapply(plots, print)
        plots
      } else {
        # print(plots[[pick]])
        plots[[pick]]
      }
    }
  }

  invisible(plots)
}
