#' Standardize key survey columns
#'
#' Harmonizes basic variables (cluster, household, strata, weight)
#' across different MICS survey files.
#'
#' @param df A data.frame from a MICS survey (bh, wm, or ch).
#' @param type Character string: one of "bh", "wm", or "ch".
#'
#' @return The modified data.frame with standardized columns:
#' \code{cluster}, \code{householdID}, \code{strata}, and \code{weight}
#' where available.
#' @export
standardize_columns <- function(df, type) {
  # ---- Cluster ----
  cand_cluster <- c("HH1", "cluster", "clusterno")
  hit <- intersect(cand_cluster, names(df))
  if (length(hit) > 0) names(df)[names(df) == hit[1]] <- "cluster"

  # ---- Household ----
  cand_hh <- c("HH2", "householdID", "hhno")
  hit <- intersect(cand_hh, names(df))
  if (length(hit) > 0) names(df)[names(df) == hit[1]] <- "householdID"

  # ---- Strata (urban/rural) ----
  if (!"strata" %in% names(df)) df$strata <- NA_character_

  if ("HH6" %in% names(df)) {
    if (length(unique(df$HH6)) > 2 && "hh6a" %in% names(df)) {
      df$strata <- ifelse(df$hh6a == 1, "urban", "rural")
    } else {
      df$strata <- ifelse(df$HH6 == 1, "urban", "rural")
    }
  }

  # ---- Weights ----
  if (type %in% c("bh", "wm") && "wmweight" %in% names(df)) df$weight <- df$wmweight
  if (type == "ch" && "chweight" %in% names(df)) df$weight <- df$chweight

  df
}

.normalize_indicator <- function(indicator) {
  if (missing(indicator) || is.null(indicator)) {
    stop(
      "process_indicator(): `indicator` must be one of NMR, ANC, or DTP3.",
      call. = FALSE
    )
  }

  if (!is.character(indicator)) {
    indicator <- deparse(substitute(indicator))
  }

  indicator <- toupper(trimws(indicator))

  if (!indicator %in% c("NMR", "ANC", "DTP3")) {
    stop(
      "process_indicator(): unsupported indicator `", indicator,
      "`. Use one of NMR, ANC, or DTP3.",
      call. = FALSE
    )
  }

  indicator
}

#' Process a supported MICS indicator into a standardized dataset
#'
#' Unified interface for processing supported indicators:
#' \code{NMR}, \code{ANC}, and \code{DTP3}.
#'
#' @param data A raw MICS survey data.frame.
#' @param indicator Indicator name. Can be supplied as
#' \code{NMR}, \code{ANC}, \code{DTP3}, or as a character string.
#'
#' @return A processed data.frame with columns
#' \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#' \code{value}, and \code{indicator}.
#' @export
process_indicator <- function(data, indicator) {
  indicator_name <- .normalize_indicator(indicator)

  out <- switch(
    indicator_name,
    "NMR"  = process_NMR(data),
    "ANC"  = process_ANC(data),
    "DTP3" = process_DTP3(data)
  )

  out
}

#' Process birth history dataset for NMR
#'
#' @param bh A birth history data.frame (bh.sav) from MICS.
#'
#' @return A data.frame with columns:
#' \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#' \code{value}, \code{indicator}
#' @export
process_NMR <- function(bh) {
  bh <- standardize_columns(bh, "bh")

  required <- c("cluster", "householdID", "weight", "strata", "BH5", "BH9C", "BH4C", "WDOI")
  missing <- setdiff(required, names(bh))

  if (length(missing) > 0) {
    stop(
      "process_NMR(): missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  bh <- bh |>
    dplyr::select(cluster, householdID, strata, weight, BH5, BH9C, BH4C, WDOI)

  names(bh)[names(bh) == "BH5"] <- "alive"
  names(bh)[names(bh) == "BH9C"] <- "age_death"
  names(bh)[names(bh) == "BH4C"] <- "birth_cmc"
  names(bh)[names(bh) == "WDOI"] <- "int_cmc"

  bh <- bh |>
    dplyr::filter(.data$birth_cmc >= .data$int_cmc - 61,
                  .data$birth_cmc <= .data$int_cmc - 1)

  bh$value <- ifelse(!is.na(bh$age_death) & bh$age_death == 0, 1, 0)
  bh$indicator <- "nmr"

  bh |>
    dplyr::select(cluster, householdID, weight, strata, value, indicator)
}

#' Process women's dataset for ANC 4+ visits
#'
#' @param wm A women's data.frame (wm.sav) from MICS.
#'
#' @return A data.frame with columns:
#' \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#' \code{value}, \code{indicator}
#' @export
process_ANC <- function(wm) {
  wm <- standardize_columns(wm, "wm")

  if ("CM17" %in% names(wm)) {
    births_var <- "CM17"
  } else if ("CM13" %in% names(wm)) {
    births_var <- "CM13"
    if (is.character(wm[[births_var]])) {
      wm[[births_var]] <- ifelse(wm[[births_var]] %in% c("Y", "1"), 1, 0)
    }
  } else {
    stop("process_ANC(): No birth recode variable (CM17 or CM13) found.", call. = FALSE)
  }

  if ("MN3" %in% names(wm)) {
    anc_var <- "MN3"
  } else if ("MN5" %in% names(wm)) {
    anc_var <- "MN5"
  } else {
    stop("process_ANC(): No ANC variable (MN3 or MN5) found.", call. = FALSE)
  }

  required <- c("cluster", "householdID", "weight", "strata", births_var, anc_var)
  missing <- setdiff(required, names(wm))

  if (length(missing) > 0) {
    stop(
      "process_ANC(): missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  out <- wm |>
    dplyr::filter(.data[[births_var]] == 1) |>
    dplyr::mutate(
      anc_visits = as.numeric(.data[[anc_var]]),
      value = dplyr::case_when(
        !is.na(.data$anc_visits) & .data$anc_visits >= 4 & .data$anc_visits < 98 ~ 1,
        !is.na(.data$anc_visits) ~ 0,
        TRUE ~ 0
      ),
      indicator = "anc"
    ) |>
    dplyr::select(cluster, householdID, weight, strata, value, indicator)

  out
}

#' Process DTP3 (Pentavalent/DTP 3rd dose) indicator from MICS child dataset
#'
#' Thailand MICS uses different variable names (IM6DH3D / IM29 / IM30).
#'
#' @param ch A MICS child-level dataset (typically ch.sav read via haven)
#'
#' @return A data.frame with columns:
#' \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#' \code{value}, \code{indicator}
#' @export
process_DTP3 <- function(ch) {
  ch <- standardize_columns(ch, "ch")
  ch <- haven::zap_labels(ch)

  is_thailand <- ("IM6DH3D" %in% names(ch)) &&
    ("IM29" %in% names(ch)) &&
    ("IM30" %in% names(ch))

  if (is_thailand) {
    age_var <- "CAGE_6"
    vax_day <- "IM6DH3D"
    vax_ever <- "IM29"
    vax_times <- "IM30"
  } else {
    age_var <- "CAGE_6"
    cand_vax_day <- c("IM6PENTA3D", "IM3PENTA3D", "IM6DTP3D")
    cand_vax_ever <- c("IM20", "IM12A")
    cand_vax_times <- c("IM21", "IM12B")

    vax_day <- cand_vax_day[cand_vax_day %in% names(ch)][1]
    vax_ever <- cand_vax_ever[cand_vax_ever %in% names(ch)][1]
    vax_times <- cand_vax_times[cand_vax_times %in% names(ch)][1]
  }

  needed <- c("cluster", "householdID", "weight", "strata", age_var, vax_day, vax_ever, vax_times)
  missing <- setdiff(needed, names(ch))

  if (length(missing) > 0) {
    stop(
      "process_DTP3(): missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  out <- ch |>
    dplyr::filter(.data[[age_var]] == 3) |>
    dplyr::mutate(
      value = dplyr::case_when(
        .data[[vax_day]] %in% c(1:31, 44) ~ 1,
        .data[[vax_ever]] == 1 & .data[[vax_times]] >= 3 & .data[[vax_times]] < 8 ~ 1,
        TRUE ~ 0
      ),
      indicator = "dtp3"
    ) |>
    dplyr::select(cluster, householdID, weight, strata, value, indicator)

  out
}
