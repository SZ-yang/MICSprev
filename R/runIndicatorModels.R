#' Run multiple SAE models for one processed indicator
#'
#' Convenience wrapper to run direct estimates, Fay-Herriot models,
#' and cluster-level models across one or more admin levels.
#'
#' @param data Processed indicator data. Must contain at least
#'   cluster, householdID, weight, strata, value.
#' @param geo Output list from process_geo_mics().
#' @param admin_levels Integer vector of admin levels to fit, e.g. c(0, 1, 2).
#' @param models Character vector of models to fit. Supported values are
#'   "direct", "fh", "cluster", "cluster_strat".
#'   "fh" fits a Fay-Herriot model; use \code{fix_var} and \code{nested}
#'   to choose the FH variant.
#' @param model Spatial random effect for fh/cluster models. Usually "bym2" or "iid".
#' @param ci Credible/confidence interval level.
#' @param fix_var Logical; passed to \code{surveyPrev::fhModel()} as
#'   \code{var.fix}. If TRUE, fit FH with fixed direct variances.
#' @param nested Logical; passed to \code{surveyPrev::fhModel()} as
#'   \code{nested}. If TRUE, fit the nested FH variant.
#' @param drop_lowvar_admin2 Logical; if TRUE, optionally drops problematic admin2
#'   units before fitting standard admin2 FH models. This is only applied when
#'   \code{fix_var = FALSE} and \code{nested = FALSE}.
#' @param lowvar_cutoff Numeric cutoff for very small direct variance at admin2.
#' @param verbose Logical; print progress messages.
#' @param ... Additional arguments passed to surveyPrev::clusterModel().
#'
#' @return A list with fitted objects by admin level and model type, plus
#'   any dropped admin2 units/clusters for the standard FH admin2 fit.
#' @export
run_indicator_models <- function(
    data,
    geo,
    admin_levels = c(0, 1, 2),
    models = c("direct", "fh", "cluster", "cluster_strat"),
    model = "bym2",
    ci = 0.95,
    fix_var = FALSE,
    nested = FALSE,
    drop_lowvar_admin2 = TRUE,
    lowvar_cutoff = 1e-30,
    verbose = TRUE,
    ...
) {
  # -----------------------------
  # checks
  # -----------------------------
  needed_cols <- c("cluster", "householdID", "weight", "strata", "value", "indicator")
  miss_cols <- setdiff(needed_cols, names(data))
  if (length(miss_cols) > 0) {
    stop(
      "run_indicator_models(): input `data` is missing required columns: ",
      paste(miss_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.list(geo)) {
    stop(
      "run_indicator_models(): `geo` must be the output list from process_geo_mics().",
      call. = FALSE
    )
  }

  if (is.null(geo$cluster.info)) {
    stop("run_indicator_models(): `geo$cluster.info` is missing.", call. = FALSE)
  }

  supported_models <- c("direct", "fh", "cluster", "cluster_strat")
  bad_models <- setdiff(models, supported_models)
  if (length(bad_models) > 0) {
    if ("fh_nested" %in% bad_models) {
      stop(
        "run_indicator_models(): `fh_nested` is no longer a separate model. ",
        "Use models = 'fh' with nested = TRUE and fix_var = TRUE instead.",
        call. = FALSE
      )
    }

    stop(
      "run_indicator_models(): unsupported model(s): ",
      paste(bad_models, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.logical(fix_var) || length(fix_var) != 1L || is.na(fix_var)) {
    stop("run_indicator_models(): `fix_var` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(nested) || length(nested) != 1L || is.na(nested)) {
    stop("run_indicator_models(): `nested` must be TRUE or FALSE.", call. = FALSE)
  }

  admin_levels <- sort(unique(as.integer(admin_levels)))
  if (any(is.na(admin_levels)) || any(!admin_levels %in% c(0L, 1L, 2L))) {
    stop(
      "run_indicator_models(): `admin_levels` must be drawn from 0, 1, 2.",
      call. = FALSE
    )
  }

  indicator_vals <- unique(stats::na.omit(as.character(data$indicator)))

  if (length(indicator_vals) != 1) {
    stop(
      "run_indicator_models(): `data$indicator` must contain exactly one indicator.",
      call. = FALSE
    )
  }

  # -----------------------------
  # dependency checks
  # -----------------------------
  if (!requireNamespace("surveyPrev", quietly = TRUE)) {
    stop(
      "run_indicator_models(): package `surveyPrev` is required but not installed.",
      call. = FALSE
    )
  }

  # directEST, fhModel, clusterModel are workflow-level functions from surveyPrev
  direct_fun  <- surveyPrev::directEST
  fh_fun      <- surveyPrev::fhModel
  cluster_fun <- surveyPrev::clusterModel

  # -----------------------------
  # helpers
  # -----------------------------
  .msg <- function(...) {
    if (isTRUE(verbose)) message(...)
  }

  .get_admin_info <- function(admin) {
    if (admin == 1L) {
      if (is.null(geo$admin.info1)) {
        stop("run_indicator_models(): geo$admin.info1 is missing.", call. = FALSE)
      }
      return(geo$admin.info1)
    }

    if (admin == 2L) {
      if (is.null(geo$admin.info2)) {
        stop("run_indicator_models(): geo$admin.info2 is missing.", call. = FALSE)
      }
      return(geo$admin.info2)
    }

    NULL
  }

  .drop_bad_admin2_clusters <- function(data_in, geo, cutoff) {
    direct_adm2 <- direct_fun(
      data = data_in,
      cluster.info = geo$cluster.info,
      admin = 2
    )

    if (is.null(direct_adm2$res.admin2)) {
      return(list(
        data = data_in,
        bad_admin2 = character(0),
        bad_clusters = numeric(0),
        direct_adm2 = direct_adm2
      ))
    }

    adm2_tab <- direct_adm2$res.admin2
    if (!("direct.var" %in% names(adm2_tab))) {
      return(list(
        data = data_in,
        bad_admin2 = character(0),
        bad_clusters = numeric(0),
        direct_adm2 = direct_adm2
      ))
    }

    bad_admin2 <- adm2_tab$admin2.name.full[
      !is.na(adm2_tab$direct.var) & adm2_tab$direct.var < cutoff
    ]

    bad_admin2 <- unique(as.character(stats::na.omit(bad_admin2)))

    if (length(bad_admin2) == 0) {
      return(list(
        data = data_in,
        bad_admin2 = character(0),
        bad_clusters = numeric(0),
        direct_adm2 = direct_adm2
      ))
    }

    if (is.null(geo$cluster.info$data) ||
        !all(c("cluster", "admin2.name.full") %in% names(geo$cluster.info$data))) {
      stop(
        "run_indicator_models(): geo$cluster.info$data must contain `cluster` ",
        "and `admin2.name.full` to drop problematic admin2 clusters.",
        call. = FALSE
      )
    }

    bad_clusters <- unique(
      geo$cluster.info$data$cluster[
        geo$cluster.info$data$admin2.name.full %in% bad_admin2
      ]
    )

    data_out <- data_in[!(data_in$cluster %in% bad_clusters), , drop = FALSE]

    list(
      data = data_out,
      bad_admin2 = bad_admin2,
      bad_clusters = bad_clusters,
      direct_adm2 = direct_adm2
    )
  }

  # -----------------------------
  # storage
  # -----------------------------
  out <- list(
    indicator = indicator_vals[1],
    fits = setNames(vector("list", length(admin_levels)), as.character(admin_levels)),
    dropped = list()
  )

  # -----------------------------
  # main loop
  # -----------------------------
  for (adm in admin_levels) {
    adm_chr <- as.character(adm)
    out$fits[[adm_chr]] <- list(
      direct = NULL,
      fh = NULL,
      cluster = NULL,
      cluster_strat = NULL
    )

    .msg("=== Fitting admin level ", adm, " ===")

    admin.info <- .get_admin_info(adm)

    # -------------------------
    # direct
    # -------------------------
    if ("direct" %in% models) {
      .msg("  - direct estimates")
      out$fits[[adm_chr]]$direct <- direct_fun(
        data = data,
        cluster.info = geo$cluster.info,
        admin = adm
      )
    }

    # -------------------------
    # FH (configurable; optionally drop bad admin2 before standard fit)
    # -------------------------
    if ("fh" %in% models) {
      if (adm == 0L) {
        .msg("  - skipping FH at admin 0")
      } else {
        fh_data <- data
        use_fh_drop <- adm == 2L &&
          isTRUE(drop_lowvar_admin2) &&
          !isTRUE(fix_var) &&
          !isTRUE(nested)

        if (adm == 2L && isTRUE(drop_lowvar_admin2) && !use_fh_drop) {
          .msg("  - not applying drop_lowvar_admin2 because fix_var or nested is TRUE")
        }

        if (use_fh_drop) {
          .msg("  - checking problematic admin2 direct variances")
          dropped_info <- .drop_bad_admin2_clusters(
            data_in = data,
            geo = geo,
            cutoff = lowvar_cutoff
          )
          fh_data <- dropped_info$data
          out$dropped$fh_admin2 <- list(
            bad_admin2 = dropped_info$bad_admin2,
            bad_clusters = dropped_info$bad_clusters
          )

          if (length(dropped_info$bad_clusters) > 0) {
            .msg(
              "    dropped ", length(dropped_info$bad_clusters),
              " cluster(s) from admin2 FH fit"
            )
          }
        }

        .msg(
          "  - Fay-Herriot model",
          " (var.fix = ", isTRUE(fix_var),
          ", nested = ", isTRUE(nested), ")"
        )
        out$fits[[adm_chr]]$fh <- fh_fun(
          data = fh_data,
          cluster.info = geo$cluster.info,
          admin.info = admin.info,
          admin = adm,
          model = model,
          aggregation = TRUE,
          var.fix = isTRUE(fix_var),
          nested = isTRUE(nested),
          CI = ci
        )
      }
    }

    # -------------------------
    # cluster (unstratified)
    # -------------------------
    if ("cluster" %in% models) {
      if (adm == 0L) {
        .msg("  - skipping cluster model at admin 0")
      } else {
        .msg("  - cluster-level model")
        out$fits[[adm_chr]]$cluster <- cluster_fun(
          data = data,
          cluster.info = geo$cluster.info,
          admin.info = admin.info,
          stratification = FALSE,
          model = model,
          admin = adm,
          aggregation = TRUE,
          CI = ci,
          ...
        )
      }
    }

    # -------------------------
    # cluster (stratified)
    # -------------------------
    if ("cluster_strat" %in% models) {
      if (adm == 0L) {
        .msg("  - skipping stratified cluster model at admin 0")
      } else {
        .msg("  - stratified cluster-level model")
        out$fits[[adm_chr]]$cluster_strat <- cluster_fun(
          data = data,
          cluster.info = geo$cluster.info,
          admin.info = admin.info,
          stratification = TRUE,
          model = model,
          admin = adm,
          aggregation = TRUE,
          CI = ci,
          ...
        )
      }
    }
  }

  class(out) <- c("mics_indicator_models", class(out))
  out
}
