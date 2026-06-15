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
#' @param transform Character. Transformation applied to the plotted mean values after
#'   applying \code{scale}. Use \code{"none"} for the original scale or
#'   \code{"log1p"} for \code{log1p(mean * scale)}. The CV map is not transformed.
#'   When \code{transform = "log1p"}, the legend labels are shown on the original scale.
#' @param direction Integer, either \code{-1} or \code{1}. Controls the color scale
#'   direction for the mean map. Use \code{-1} (default) when darker = higher value
#'   (e.g. NMR), and \code{1} when darker = lower value (e.g. ANC, DTP3).
#' @param out_dir Optional directory to save PDFs. If NULL, does not save.
#' @param prefix Optional filename prefix.
#' @param map_ncol Integer. Number of columns in the faceted map layout.
#'   Default is 2, so four models will appear as a 2 by 2 layout.
#'
#' @return A list with ggplot objects: \code{list(mean_map = p_mean, cv_map = p_cv)}.
#' @export
plot_indicator_maps <- function(
    fit,
    geo,
    admin = 1,
    models = c("direct", "fh", "fh_nested", "cluster", "cluster_strat"),
    scale = 1,
    transform = c("none", "log1p"),
    direction = -1,
    out_dir = NULL,
    prefix = NULL,
    map_ncol = 2
) {
  admin <- as.integer(admin)
  transform <- match.arg(transform)

  if (!admin %in% c(1L, 2L)) {
    stop("plot_indicator_maps(): admin must be 1 or 2.", call. = FALSE)
  }

  if (!is.numeric(direction) || length(direction) != 1 || !direction %in% c(-1, 1)) {
    stop("plot_indicator_maps(): direction must be either -1 or 1.", call. = FALSE)
  }

  if (!is.numeric(map_ncol) || length(map_ncol) != 1 || is.na(map_ncol) || map_ncol < 1) {
    stop("plot_indicator_maps(): map_ncol must be a positive integer.", call. = FALSE)
  }

  direction <- as.integer(direction)
  map_ncol <- as.integer(map_ncol)

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

  mean_label <- if (scale == 1) "Mean" else paste0("Mean (x", scale, ")")

  dat_long$mean_raw <- dat_long$mean * scale
  dat_long$mean_plot <- dat_long$mean_raw

  if (transform == "log1p") {
    if (any(dat_long$mean_raw < 0, na.rm = TRUE)) {
      stop("plot_indicator_maps(): transform = 'log1p' requires non-negative mean values after scaling.", call. = FALSE)
    }

    # Plot colors on log1p scale, but keep the legend on the original scale.
    dat_long$mean_plot <- log1p(dat_long$mean_raw)
  }

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

  n_models <- length(unique(dat_long$model))
  map_ncol <- min(map_ncol, n_models)
  map_nrow <- ceiling(n_models / map_ncol)

  p_mean <- mapPlot_fun(
    data = dat_long,
    geo = geo_sf,
    by.data = "area",
    by.geo = by_geo,
    is.long = TRUE,
    variables = "model",
    values = "mean_plot",
    legend.label = mean_label,
    direction = direction,
    ncol = map_ncol
  )

  # For log1p maps, replace the legend scale so that tick labels are shown
  # on the original scale, not the log1p scale.
  if (transform == "log1p") {
    mean_raw_finite <- dat_long$mean_raw[is.finite(dat_long$mean_raw)]

    if (length(mean_raw_finite) > 0) {
      raw_breaks <- pretty(mean_raw_finite, n = 5)
      raw_breaks <- raw_breaks[
        raw_breaks >= min(mean_raw_finite, na.rm = TRUE) &
          raw_breaks <= max(mean_raw_finite, na.rm = TRUE)
      ]

      if (length(raw_breaks) > 0) {
        p_mean <- p_mean +
          ggplot2::scale_fill_viridis_c(
            name = mean_label,
            breaks = log1p(raw_breaks),
            labels = format(raw_breaks, trim = TRUE, scientific = FALSE),
            direction = direction
          )
      }
    }
  }

  p_cv <- mapPlot_fun(
    data = dat_long,
    geo = geo_sf,
    by.data = "area",
    by.geo = by_geo,
    is.long = TRUE,
    variables = "model",
    values = "cv",
    legend.label = "CV",
    direction = -1,
    ncol = map_ncol
  )

  plot_width <- if (map_ncol == 1) 8 else 10
  plot_height <- 4 * map_nrow + 1

  .maybe_save_plot(
    p_mean, out_dir,
    filename = sprintf(
      "%s_admin%d_mean%s.pdf",
      prefix,
      admin,
      if (transform == "none") "" else paste0("_", transform)
    ),
    width = plot_width,
    height = plot_height
  )

  .maybe_save_plot(
    p_cv, out_dir,
    filename = sprintf("%s_admin%d_cv.pdf", prefix, admin),
    width = plot_width,
    height = plot_height
  )

  list(mean_map = p_mean, cv_map = p_cv)
}
