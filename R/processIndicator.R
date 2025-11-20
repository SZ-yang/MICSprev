#' Standardize key survey columns
#'
#' Harmonizes basic variables (cluster, household, strata, weight)
#' across different MICS survey files.
#'
#' @param df A data.frame from a MICS survey (bh, wm, or ch).
#' @param type Character string: one of "bh", "wm", or "ch".
#'
#' @return The modified data.frame with standardized columns:
#'   \code{cluster}, \code{householdID}, \code{strata}, and \code{weight}
#'   where available.
#' @export
standardize_columns <- function(df, type) {
  # Cluster
  cand_cluster <- c("HH1", "cluster", "clusterno")
  hit <- intersect(cand_cluster, names(df))
  if (length(hit) > 0) {
    names(df)[names(df) == hit[1]] <- "cluster"
  }

  # Household
  cand_hh <- c("HH2", "householdID", "hhno")
  hit <- intersect(cand_hh, names(df))
  if (length(hit) > 0) {
    names(df)[names(df) == hit[1]] <- "householdID"
  }

  # Area (urban/rural)
  if ("HH6" %in% names(df)) {
    if (length(unique(df$HH6)) > 2 && "hh6a" %in% names(df)) {
      df$strata <- ifelse(df$hh6a == 1, "urban", "rural")
    } else {
      df$strata <- ifelse(df$HH6 == 1, "urban", "rural")
    }
  }

  # Weights
  if (type == "bh" && "wmweight" %in% names(df)) df$weight <- df$wmweight
  if (type == "wm" && "wmweight" %in% names(df)) df$weight <- df$wmweight
  if (type == "ch" && "chweight" %in% names(df)) df$weight <- df$chweight

  df
}

#' Process birth history dataset for NMR
#'
#' Converts a MICS birth history dataset into a standardized format
#' suitable for neonatal mortality analysis.
#'
#' @param bh A birth history data.frame (bh.sav) from MICS.
#'
#' @return A data.frame with columns:
#'   \code{cluster}, \code{householdID}, \code{weight}, \code{strata}, \code{value},
#'   where \code{value = 1} indicates a neonatal death and \code{0} otherwise.
#' @export
process_NMR <- function(bh) {
  bh <- standardize_columns(bh, "bh")

  # Variables (harmonize names)
  cand_alive     <- c("BH5")
  cand_age_death <- c("BH9C")
  cand_birth_cmc <- c("BH4C")
  cand_int_cmc   <- c("WDOI")

  # Subset variables
  bh <- bh |>
    dplyr::select(
      cluster, householdID, strata, weight,
      dplyr::all_of(c(cand_alive, cand_age_death, cand_birth_cmc, cand_int_cmc))
    )

  names(bh)[names(bh) %in% cand_alive]     <- "alive"
  names(bh)[names(bh) %in% cand_age_death] <- "age_death"
  names(bh)[names(bh) %in% cand_birth_cmc] <- "birth_cmc"
  names(bh)[names(bh) %in% cand_int_cmc]   <- "int_cmc"

  # Restrict to births in last 5 years (61 months)
  bh <- bh |>
    dplyr::filter(.data$birth_cmc >= .data$int_cmc - 61,
                  .data$birth_cmc <= .data$int_cmc - 1)

  # Define outcome: neonatal death (age_death == 0)
  bh$value <- ifelse(!is.na(bh$age_death) & bh$age_death == 0, 1, 0)

  bh |>
    dplyr::select(cluster, householdID, weight, strata, value)
}

#' Process women's dataset for ANC 4+ visits
#'
#' Converts a MICS women's dataset into a standardized format
#' with an indicator for at least four ANC visits.
#'
#' @param wm A women's data.frame (wm.sav) from MICS.
#'
#' @return A data.frame with columns:
#'   \code{cluster}, \code{householdID}, \code{weight}, \code{strata}, \code{value},
#'   where \code{value = 1} indicates >= 4 ANC visits.
#' @export
process_ANC <- function(wm) {
  wm <- standardize_columns(wm, "wm")

  # --- Birth recode (CM17 vs CM13)
  if ("CM17" %in% names(wm)) {
    births_var <- "CM17"
  } else if ("CM13" %in% names(wm)) {
    births_var <- "CM13"
    # Convert Y/N or 0/1 into numeric
    if (is.character(wm[[births_var]])) {
      wm[[births_var]] <- ifelse(wm[[births_var]] %in% c("Y", "1"), 1, 0)
    }
  } else {
    stop("No birth recode variable (CM17 or CM13) found.")
  }

  # --- ANC visits (MN3 vs MN5)
  # Nigeria 2016–17 has MN3, and MN5 there is "has immunization card"
  if ("MN3" %in% names(wm)) {
    anc_var <- "MN3"
  } else if ("MN5" %in% names(wm)) {
    anc_var <- "MN5"
  } else {
    stop("No ANC variable (MN3 or MN5) found.")
  }

  wm <- wm |>
    dplyr::filter(.data[[births_var]] == 1) |>
    dplyr::mutate(
      anc_visits = as.numeric(.data[[anc_var]]),
      value = dplyr::case_when(
        !is.na(.data$anc_visits) & .data$anc_visits >= 4 & .data$anc_visits < 98 ~ 1,
        !is.na(.data$anc_visits) ~ 0,
        TRUE ~ 0
      )
    )

  wm |>
    dplyr::select(cluster, householdID, weight, strata, value)
}

#' Process children's dataset for DTP3 coverage
#'
#' Converts a MICS children's dataset into a standardized format
#' with an indicator for DTP3 vaccination.
#'
#' @param ch A children's data.frame (ch.sav) from MICS.
#'
#' @return A data.frame with columns:
#'   \code{cluster}, \code{householdID}, \code{weight}, \code{strata}, \code{value},
#'   where \code{value = 1} indicates DTP3 vaccination.
#' @export
process_DTP3 <- function(ch) {
  ch <- standardize_columns(ch, "ch")

  # Harmonize vaccine variables
  cand_age       <- c("CAGE_6")
  cand_vax_day   <- c("IM6PENTA3D", "IM3PENTA3D")
  cand_vax_ever  <- c("IM20", "IM12A")
  cand_vax_times <- c("IM21", "IM12B")

  age_var   <- cand_age[cand_age %in% names(ch)][1]
  vax_day   <- cand_vax_day[cand_vax_day %in% names(ch)][1]
  vax_ever  <- cand_vax_ever[cand_vax_ever %in% names(ch)][1]
  vax_times <- cand_vax_times[cand_vax_times %in% names(ch)][1]

  ch <- ch |>
    dplyr::filter(.data[[age_var]] == 3) |>
    dplyr::mutate(
      value = dplyr::case_when(
        .data[[vax_day]] %in% c(1:31, 44) ~ 1,
        .data[[vax_ever]] == 1 &
          .data[[vax_times]] >= 3 &
          .data[[vax_times]] < 8 ~ 1,
        TRUE ~ 0
      )
    )

  ch |>
    dplyr::select(cluster, householdID, weight, strata, value)
}
