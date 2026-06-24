# R/plot_helpers.R
# Internal helpers for extracting model outputs and saving plots

#' @keywords internal
.get_summer_fun <- function(fname) {
  if (!requireNamespace("SUMMER", quietly = TRUE)) {
    stop("Required package 'SUMMER' is not installed.", call. = FALSE)
  }

  # Try exported first
  out <- tryCatch(getExportedValue("SUMMER", fname), error = function(e) NULL)
  if (is.function(out)) return(out)

  # Fallback: internal namespace
  out <- tryCatch(getFromNamespace(fname, "SUMMER"), error = function(e) NULL)
  if (is.function(out)) return(out)

  # Final fallback: attached search path
  if (exists(fname, mode = "function")) {
    out <- get(fname, mode = "function")
    if (is.function(out)) return(out)
  }

  stop(
    "Could not find function '", fname, "' in package SUMMER.\n",
    "Your SUMMER version may be incompatible.",
    call. = FALSE
  )
}

#' @keywords internal
.ensure_admin2_name_full <- function(admin2_sf) {
  if (!inherits(admin2_sf, "sf")) return(admin2_sf)

  if (!("admin2.name.full" %in% names(admin2_sf))) {
    if (all(c("NAME_1", "NAME_2") %in% names(admin2_sf))) {
      admin2_sf$admin2.name.full <- paste0(admin2_sf$NAME_1, "_", admin2_sf$NAME_2)
    }
  }
  admin2_sf
}

#' @keywords internal
.pick_area_col <- function(df, admin) {
  cand <- c(
    paste0("admin", admin, ".name.full"),
    paste0("admin", admin, ".name"),
    paste0("admin", admin, ".name_full"),
    paste0("adm", admin, ".name.full"),
    paste0("adm", admin, ".name")
  )
  hit <- intersect(cand, names(df))
  if (length(hit) > 0) return(hit[1])

  hit2 <- grep(paste0("^admin", admin, "\\."), names(df), value = TRUE)
  if (length(hit2) > 0) return(hit2[1])

  stop("Could not find area id column in model output for admin=", admin, call. = FALSE)
}

#' @keywords internal
.extract_one_model_admin <- function(fit, admin, model_name) {
  adm_chr <- as.character(admin)

  if (is.null(fit$fits[[adm_chr]])) return(NULL)
  obj <- fit$fits[[adm_chr]][[model_name]]
  if (is.null(obj)) return(NULL)

  if (model_name == "direct") {
    res_name <- paste0("res.admin", admin)
    res <- obj[[res_name]]
    if (is.null(res)) return(NULL)

    area_col <- .pick_area_col(res, admin)

    if (!("direct.est" %in% names(res))) {
      stop("directEST output missing `direct.est`.", call. = FALSE)
    }
    if (!("cv" %in% names(res))) {
      if (all(c("direct.se", "direct.est") %in% names(res))) {
        res$cv <- with(res, direct.se / direct.est)
      } else {
        res$cv <- NA_real_
      }
    }

    if (all(c("direct.lower", "direct.upper") %in% names(res))) {
      res$interval_width <- res$direct.upper - res$direct.lower
    } else {
      res$interval_width <- NA_real_
    }

    out <- res[, c(area_col, "direct.est", "cv", "interval_width"), drop = FALSE]
    names(out) <- c("area", "mean", "cv", "interval_width")

    out$area_key <- area_col
    out$model <- "direct"
    return(out)
  }

  if (!model_name %in% c("fh", "cluster", "cluster_strat")) {
    stop(
      "Unsupported model_name in .extract_one_model_admin(): ",
      model_name,
      call. = FALSE
    )
  }

  res_name <- paste0("res.admin", admin)
  res <- obj[[res_name]]
  if (is.null(res)) return(NULL)

  if (model_name == "cluster_strat" && "type" %in% names(res)) {
    res <- res[res$type == "full", , drop = FALSE]
  }

  area_col <- .pick_area_col(res, admin)

  if (!("mean" %in% names(res))) {
    stop("SUMMER model output missing `mean`.", call. = FALSE)
  }
  if (!("cv" %in% names(res))) {
    if (all(c("sd", "mean") %in% names(res))) {
      res$cv <- with(res, sd / mean)
    } else {
      res$cv <- NA_real_
    }
  }

  if (all(c("lower", "upper") %in% names(res))) {
    res$interval_width <- res$upper - res$lower
  } else {
    res$interval_width <- NA_real_
  }

  out <- res[, c(area_col, "mean", "cv", "interval_width"), drop = FALSE]
  names(out) <- c("area", "mean", "cv", "interval_width")

  out$area_key <- area_col
  out$model <- model_name
  out
}

#' @keywords internal
.extract_models_long <- function(fit, admin, models) {
  pieces <- lapply(models, function(m) .extract_one_model_admin(fit, admin, m))
  pieces <- Filter(Negate(is.null), pieces)
  if (length(pieces) == 0) return(NULL)
  out <- do.call(rbind, pieces)
  out$model <- factor(out$model, levels = unique(out$model))
  out
}

#' @keywords internal
.resolve_geo_join <- function(geo, admin, area_col_in_data) {
  if (!is.list(geo)) stop("geo must be a list from process_geo_mics()", call. = FALSE)

  if (admin == 1) {
    geo_sf <- geo$admin1
    if (is.null(geo_sf)) stop("geo$admin1 is missing.", call. = FALSE)

    by_geo <- if (area_col_in_data %in% names(geo_sf)) area_col_in_data else "NAME_1"
    return(list(geo_sf = geo_sf, by_geo = by_geo))
  }

  if (admin == 2) {
    geo_sf <- geo$admin2
    if (is.null(geo_sf)) stop("geo$admin2 is missing.", call. = FALSE)

    geo_sf <- .ensure_admin2_name_full(geo_sf)

    if (area_col_in_data %in% names(geo_sf)) {
      by_geo <- area_col_in_data
    } else if ("admin2.name.full" %in% names(geo_sf)) {
      by_geo <- "admin2.name.full"
    } else {
      by_geo <- "NAME_2"
    }

    return(list(geo_sf = geo_sf, by_geo = by_geo))
  }

  stop("admin must be 1 or 2 for map plotting.", call. = FALSE)
}

#' @keywords internal
.maybe_save_plot <- function(p, out_dir, filename, width = 8, height = 6) {
  if (is.null(out_dir)) return(invisible(NULL))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = file.path(out_dir, filename),
    plot = p,
    width = width,
    height = height
  )
  invisible(NULL)
}
