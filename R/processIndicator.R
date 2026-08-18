#' Standardize key survey columns
#'
#' Harmonizes basic variables (cluster, household, strata, weight)
#' across different MICS survey files.
#'
#' @param df A data.frame from a MICS survey recode file.
#' @param type Character string identifying the recode file: one of
#'   \code{"bh"} (birth history), \code{"wm"} (women), \code{"ch"} (children
#'   under 5), \code{"hh"} (households), \code{"hl"} (household listing),
#'   \code{"fs"} (children age 5-17) or \code{"mn"} (men).
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
  # FLAG: Nigeria MICS6 (2021) ships two weight sets because the release merges
  # the MICS and NICS samples (Readme_Nigeria_MICS6: "the appropriate data files
  # have two sets of weights ... one for MICS and another one for MICS and NICS
  # combined"). The unsuffixed names (hhweight/chweight/...) are the combined
  # MICS-NICS weights; the *MICS-suffixed names are MICS-only. We take the
  # unsuffixed names throughout, matching both the standard MICS6 tabulation
  # syntax and the pre-existing behaviour of process_DTP3()/process_NMR().
  weight_var <- switch(
    type,
    "bh" = "wmweight",
    "wm" = "wmweight",
    "ch" = "chweight",
    "hh" = "hhweight",
    "hl" = "hhweight",
    "fs" = "fsweight",
    "mn" = "mnweight",
    NULL
  )

  if (!is.null(weight_var) && weight_var %in% names(df)) df$weight <- df[[weight_var]]

  df
}

# ---------------------------------------------------------------------------
# Country gating
# ---------------------------------------------------------------------------

# Registry mapping each supported indicator to the countries it may be run on.
# "all" means the indicator is country-agnostic (the three original indicators).
# Anything else is a character vector of lowercase country names; callers must
# then pass `country` explicitly to process_indicator().
.indicator_countries <- list(
  # --- country-agnostic (pre-existing) ---
  NMR  = "all",
  ANC  = "all",
  DTP3 = "all",

  # --- Nigeria-only: Health ---
  ANC1    = "nigeria",
  SBA     = "nigeria",
  PNCNB   = "nigeria",
  PNCMOM  = "nigeria",
  PENTA1  = "nigeria",

  # --- Nigeria-only: Child Protection ---
  BIRTHREG      = "nigeria",
  CHILDMARRIAGE = "nigeria",
  CHILDLABOUR   = "nigeria",
  FGMDAUGHTER   = "nigeria",

  # --- Nigeria-only: Education ---
  OOSPRIMARY   = "nigeria",
  OOSSECONDARY = "nigeria",
  READING      = "nigeria",
  MATH         = "nigeria",

  # --- Nigeria-only: Nutrition ---
  # Stunting and severe wasting are absent by design: Nigeria MICS6 (2021)
  # released no anthropometry. See R/indicators_nga_nutrition.R.
  EBF  = "nigeria",
  VITA = "nigeria",

  # --- Nigeria-only: WASH ---
  BASICWATER = "nigeria",
  BASICSAN   = "nigeria",
  OPENDEF    = "nigeria",

  # --- Nigeria-only: Social Policy ---
  HANDWASH    = "nigeria",
  FUNCDIFF517 = "nigeria",
  FUNCDIFF217 = "nigeria",
  HEALTHINS   = "nigeria",
  SOCTRANSFER = "nigeria",
  BANKACCT    = "nigeria",
  BORROWED    = "nigeria"
)

#' Supported indicators and the countries they may be run on
#'
#' Returns the indicator/country registry used to gate country-specific
#' indicators. Indicators mapped to \code{"all"} are country-agnostic; every
#' other indicator requires a matching \code{country} argument.
#'
#' @return A named list. Names are indicator codes; values are either the
#'   string \code{"all"} or a character vector of lowercase country names.
#' @export
indicator_countries <- function() {
  .indicator_countries
}

# Normalize a free-text country name to the registry's lowercase convention.
.normalize_country <- function(country) {
  if (is.null(country)) return(NULL)

  country <- tolower(trimws(as.character(country)))
  if (length(country) != 1L || is.na(country) || !nzchar(country)) return(NULL)

  country
}

# Check an indicator against a country. `context` names the calling function so
# the error message points at the right place.
.check_indicator_country <- function(indicator_name, country, context) {
  allowed <- .indicator_countries[[indicator_name]]

  # Unregistered indicators are country-agnostic by default.
  if (is.null(allowed) || identical(allowed, "all")) return(invisible(TRUE))

  country <- .normalize_country(country)

  if (is.null(country)) {
    stop(
      context, ": indicator `", indicator_name, "` is country-specific and ",
      "requires `country`. Supported: ",
      paste(allowed, collapse = ", "), ".",
      call. = FALSE
    )
  }

  if (!country %in% allowed) {
    stop(
      context, ": indicator `", indicator_name, "` is not supported for ",
      "country `", country, "`. Supported: ",
      paste(allowed, collapse = ", "), ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.normalize_indicator <- function(indicator) {
  supported <- names(.indicator_countries)

  if (missing(indicator) || is.null(indicator)) {
    stop(
      "process_indicator(): `indicator` must be one of ",
      paste(supported, collapse = ", "), ".",
      call. = FALSE
    )
  }

  if (!is.character(indicator)) {
    indicator <- deparse(substitute(indicator))
  }

  indicator <- toupper(trimws(indicator))

  if (!indicator %in% supported) {
    stop(
      "process_indicator(): unsupported indicator `", indicator,
      "`. Use one of ", paste(supported, collapse = ", "), ".",
      call. = FALSE
    )
  }

  indicator
}

#' Process a supported MICS indicator into a standardized dataset
#'
#' Unified interface for processing supported indicators. \code{NMR},
#' \code{ANC} and \code{DTP3} are country-agnostic; every other indicator is
#' country-specific and requires \code{country} to be supplied. See
#' \code{\link{indicator_countries}()} for the registry.
#'
#' Nigeria-only indicators (Nigeria MICS6, 2021), by sector:
#' \itemize{
#'   \item Health: \code{ANC1}, \code{SBA}, \code{PNCNB}, \code{PNCMOM},
#'     \code{PENTA1}
#'   \item Child protection: \code{BIRTHREG}, \code{CHILDMARRIAGE},
#'     \code{CHILDLABOUR}, \code{FGMDAUGHTER}
#'   \item Education: \code{OOSPRIMARY}, \code{OOSSECONDARY}, \code{READING},
#'     \code{MATH}
#'   \item Nutrition: \code{EBF}, \code{VITA}
#'   \item WASH: \code{BASICWATER}, \code{BASICSAN}, \code{OPENDEF}
#'   \item Social policy: \code{HANDWASH}, \code{FUNCDIFF517},
#'     \code{FUNCDIFF217}, \code{HEALTHINS}, \code{SOCTRANSFER},
#'     \code{BANKACCT}, \code{BORROWED}
#' }
#'
#' Call \code{indicator_countries()} for the authoritative, always-current
#' list.
#'
#' Most indicators are derived from a single recode file, passed as
#' \code{data}. A few span more than one; those take their additional recodes
#' through \code{...}. \code{FGMDAUGHTER} needs all three of \code{fg.sav}
#' (as \code{data}), \code{bh.sav} and \code{wm.sav}:
#'
#' \preformatted{
#' process_indicator(fg, "FGMDAUGHTER", country = "Nigeria", bh = bh, wm = wm)
#' }
#'
#' @param data A raw MICS survey data.frame. The recode file required depends
#'   on the indicator (see the individual \code{process_*()} functions).
#' @param indicator Indicator name, supplied bare (e.g. \code{NMR}) or as a
#'   character string.
#' @param country Character country name, e.g. \code{"Nigeria"}. Required for
#'   country-specific indicators; ignored for country-agnostic ones. Matching
#'   is case- and whitespace-insensitive.
#' @param ... Additional recode data.frames for indicators that span more than
#'   one recode file. Passed on to the underlying \code{process_*()} function,
#'   which names them; supplying them for a single-file indicator is an error.
#'
#' @return A processed data.frame with columns
#' \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#' \code{value}, and \code{indicator}.
#' @export
process_indicator <- function(data, indicator, country = NULL, ...) {
  indicator_name <- .normalize_indicator(indicator)
  .check_indicator_country(indicator_name, country, "process_indicator()")

  out <- switch(
    indicator_name,
    "NMR"    = process_NMR(data),
    "ANC"    = process_ANC(data),
    "DTP3"   = process_DTP3(data),

    # --- Nigeria-only: Health ---
    "ANC1"   = process_ANC1(data),
    "SBA"    = process_SBA(data),
    "PNCNB"  = process_PNCNB(data),
    "PNCMOM" = process_PNCMOM(data),
    "PENTA1" = process_PENTA1(data),

    # --- Nigeria-only: Child Protection ---
    "BIRTHREG"      = process_BIRTHREG(data),
    "CHILDMARRIAGE" = process_CHILDMARRIAGE(data),
    "CHILDLABOUR"   = process_CHILDLABOUR(data),
    # Multi-recode: bh and wm arrive through `...`.
    "FGMDAUGHTER"   = process_FGMDAUGHTER(data, ...),

    # --- Nigeria-only: Education ---
    # `...` carries the country-specific school-structure and reading-task
    # parameters documented on each function.
    "OOSPRIMARY"   = process_OOSPRIMARY(data, ...),
    "OOSSECONDARY" = process_OOSSECONDARY(data, ...),
    "READING"      = process_READING(data, ...),
    "MATH"         = process_MATH(data),

    # --- Nigeria-only: Nutrition ---
    "EBF"  = process_EBF(data),
    "VITA" = process_VITA(data, ...),

    # --- Nigeria-only: WASH ---
    "BASICWATER" = process_BASICWATER(data),
    "BASICSAN"   = process_BASICSAN(data),
    "OPENDEF"    = process_OPENDEF(data),

    # --- Nigeria-only: Social Policy ---
    "HANDWASH"    = process_HANDWASH(data),
    "FUNCDIFF517" = process_FUNCDIFF517(data),
    # Multi-recode: ch arrives through `...`.
    "FUNCDIFF217" = process_FUNCDIFF217(data, ...),
    "HEALTHINS"   = process_HEALTHINS(data),
    # Multi-recode: hl arrives through `...`.
    "SOCTRANSFER" = process_SOCTRANSFER(data, ...),
    "BANKACCT"    = process_BANKACCT(data),
    "BORROWED"    = process_BORROWED(data, ...)
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
