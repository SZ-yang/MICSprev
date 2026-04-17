# R/plot_indicator.R
# Exported plotting functions for indicator model outputs

.get_mean_map_direction <- function(indicator) {
  if (is.null(indicator) || length(indicator) == 0 || all(is.na(indicator))) {
    warning(
      "Missing indicator in `fit$indicator`; using default direction = -1.",
      call. = FALSE
    )
    return(-1)
  }

  indicator <- tolower(as.character(indicator)[1])

  if (indicator %in% "nmr") {
    return(-1)  # high = bad -> dark
  }

  if (indicator %in% c("anc", "dtp3")) {
    return(1)   # high = good -> light
  }

  warning(
    "Unknown indicator in `fit$indicator`; using default direction = -1.",
    call. = FALSE
  )
  -1
}

#' Plot prevalence mean and CV maps for an indicator
#'
#' Creates two maps (mean and CV) across selected models at a given admin level.
#' This function does not refit any models; it only plots from \code{run_indicator_models()} output.
#'
#' @param fit Output from \code{\link{run_indicator_models}}.
#' @param geo Output list from \code{\link{process_geo_mics}}.
#' @param admin Integer. One of 1 or 2.
#' @param models Character vector of models to plot.
#'   Supported values are \code{c("direct","fh","fh_nested","cluster","cluster_strat")}.
#' @param scale Numeric multiplier applied to means before plotting (default 1).
#' @param out_dir Optional directory to save PDFs. If NULL, does not save.
#' @param prefix Optional filename prefix.
#'
#' @return A list with ggplot objects: \code{list(mean_map = p_mean, cv_map = p_cv)}.
#' @export
plot_indicator_maps <- function(
    fit,
    geo,
    admin = 1,
    models = c("direct", "fh", "fh_nested", "cluster", "cluster_strat"),
    scale = 1,
    out_dir = NULL,
    prefix = NULL
) {
  admin <- as.integer(admin)
  if (!admin %in% c(1L, 2L)) {
    stop("plot_indicator_maps(): admin must be 1 or 2.", call. = FALSE)
  }

  supported_models <- c("direct", "fh", "fh_nested", "cluster", "cluster_strat")
  bad_models <- setdiff(models, supported_models)
  if (length(bad_models) > 0) {
    stop(
      "plot_indicator_maps(): unsupported model(s): ",
      paste(bad_models, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(prefix)) prefix <- "indicator"

  mapPlot_fun <- .get_summer_fun("mapPlot")

  dat_long <- .extract_models_long(fit, admin = admin, models = models)
  if (is.null(dat_long) || nrow(dat_long) == 0) {
    stop("plot_indicator_maps(): no model results available to plot for this admin/models.", call. = FALSE)
  }

  dat_long$mean <- dat_long$mean * scale

  join_keys <- unique(dat_long$area_key)
  if (length(join_keys) != 1) {
    stop(
      "plot_indicator_maps(): inconsistent area keys across models: ",
      paste(join_keys, collapse = ", "),
      call. = FALSE
    )
  }

  geo_join <- .resolve_geo_join(geo, admin = admin, area_col_in_data = join_keys)
  geo_sf <- geo_join$geo_sf
  by_geo <- geo_join$by_geo

  mean_direction <- .get_mean_map_direction(fit$indicator)
  cv_direction <- -1

  p_mean <- mapPlot_fun(
    data = dat_long,
    geo = geo_sf,
    by.data = "area",
    by.geo = by_geo,
    is.long = TRUE,
    variables = "model",
    values = "mean",
    legend.label = if (scale == 1) "Mean" else paste0("Mean (x", scale, ")"),
    direction = mean_direction,
    ncol = length(unique(dat_long$model))
  )

  p_cv <- mapPlot_fun(
    data = dat_long,
    geo = geo_sf,
    by.data = "area",
    by.geo = by_geo,
    is.long = TRUE,
    variables = "model",
    values = "cv",
    legend.label = "CV",
    direction = cv_direction,
    ncol = length(unique(dat_long$model))
  )

  .maybe_save_plot(
    p_mean, out_dir,
    filename = sprintf("%s_admin%d_mean.pdf", prefix, admin),
    width = 12, height = 7
  )
  .maybe_save_plot(
    p_cv, out_dir,
    filename = sprintf("%s_admin%d_cv.pdf", prefix, admin),
    width = 12, height = 7
  )

  list(mean_map = p_mean, cv_map = p_cv)
}

#' Plot ridge plot of posterior distributions for selected models
#'
#' Uses \code{ridgePlot()} on SUMMER model objects (FH/cluster/cluster_strat).
#' Direct estimates are not ridge-plotted.
#'
#' @param fit Output from \code{\link{run_indicator_models}}.
#' @param admin Integer. One of 1 or 2.
#' @param models Character vector subset of
#'   \code{c("fh","fh_nested","cluster","cluster_strat")}.
#' @param threshold Optional numeric. If provided, draws a vertical reference line.
#' @param scale Numeric multiplier applied to x-axis labels only (default 1).
#' @param out_dir Optional directory to save PDF. If NULL, does not save.
#' @param prefix Optional filename prefix.
#'
#' @return A ggplot object, patchwork object, or list of ggplots.
#' @export
plot_indicator_ridge <- function(
    fit,
    admin = 1,
    models = c("fh", "fh_nested", "cluster", "cluster_strat"),
    threshold = NULL,
    scale = 1,
    out_dir = NULL,
    prefix = NULL
) {
  admin <- as.integer(admin)
  if (!admin %in% c(1L, 2L)) {
    stop("plot_indicator_ridge(): admin must be 1 or 2.", call. = FALSE)
  }

  ridge_supported <- c("fh", "fh_nested", "cluster", "cluster_strat")
  bad_models <- setdiff(models, ridge_supported)
  if (length(bad_models) > 0) {
    stop(
      "plot_indicator_ridge(): unsupported ridge model(s): ",
      paste(bad_models, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(prefix)) prefix <- "indicator"

  ridgePlot_fun <- .get_summer_fun("ridgePlot")

  adm_chr <- as.character(admin)
  plots <- list()

  for (m in models) {
    obj <- fit$fits[[adm_chr]][[m]]
    if (is.null(obj)) next

    p <- ridgePlot_fun(x = obj, direction = -1) +
      ggplot2::ggtitle(m)

    if (!is.null(threshold)) {
      p <- p + ggplot2::geom_vline(
        xintercept = threshold,
        linetype = "dashed",
        alpha = 0.6,
        linewidth = 1
      )
    }

    if (!is.null(scale) && is.numeric(scale) && length(scale) == 1 && scale != 1) {
      p <- p + ggplot2::scale_x_continuous(labels = function(x) x * scale)
    }

    plots[[m]] <- p
  }

  if (length(plots) == 0) {
    stop("plot_indicator_ridge(): no ridge-capable model objects found in `fit` for this admin/models.", call. = FALSE)
  }

  out_plot <- if (length(plots) == 1) {
    plots[[1]]
  } else {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      patchwork::wrap_plots(plots, ncol = 1)
    } else {
      warning("patchwork not installed; returning a list of ggplots.", call. = FALSE)
      return(plots)
    }
  }

  .maybe_save_plot(
    out_plot, out_dir,
    filename = sprintf("%s_admin%d_ridge.pdf", prefix, admin),
    width = 14, height = 8
  )

  out_plot
}
