# R/indicators_nga_wash.R
#
# Nigeria-only WASH indicators, translated from the official MICS6 SPSS
# tabulation syntax (mics.unicef.org/tools, "42 Syntax Files 20220509") against
# the Nigeria MICS6 (2021) recode variable names.
#
# Source files:
#   MICS6 - 10 - WS.1.2.sps  Basic drinking water services -> BASICWATER
#   MICS6 - 10 - WS.3.2.sps  Basic sanitation services     -> BASICSAN
#   MICS6 - 10 - WS.3.1.sps  Sanitation facility type      -> OPENDEF
#   MICS6 - 10 - WS.2.1.sps  Handwashing facility          -> HANDWASH
#                                                    (see indicators_nga_socialpolicy.R)
# Shared recodes live in define/MICS6 - 10 - WS.sps.


# Denominator shared by the WS tables: households with a completed interview,
# weighted by household size so the estimate is population-based.
#
#   select if (HH46 = 1).
#   compute hhweightHH48 = HH48 * hhweight.
#   weight by hhweightHH48.
#
# Every WS table reports "percentage of the household population", not
# percentage of households, so the household weight is multiplied by the
# number of household members (HH48). Downstream modelling therefore treats
# each household as HH48 population units, matching the published tables.
.nga_hh_population_denominator <- function(hh, fn, extra_required) {
  hh <- standardize_columns(hh, "hh")

  required <- c(
    "cluster", "householdID", "weight", "strata", "HH46", "HH48", extra_required
  )
  .require_columns(hh, required, fn)

  hh <- hh[.as_code(hh$HH46) %in% 1, , drop = FALSE]

  hh$weight <- .as_code(hh$weight) * .as_code(hh$HH48)

  # Households with no usable population weight cannot contribute.
  hh[!is.na(hh$weight), , drop = FALSE]
}


# `drinkingWater` from define/MICS6 - 10 - WS.sps:
#   1 = improved when WS1 is one of the listed codes, 2 = unimproved,
#   9 = missing. DK and missing are treated as unimproved, so for a binary
#   indicator only "improved" matters.
.nga_improved_water <- function(hh) {
  .as_code(hh$WS1) %in% c(11, 12, 13, 14, 21, 31, 41, 51, 61, 71, 72, 91, 92)
}


# `toiletType` from define/MICS6 - 10 - WS.sps:
#   1 = improved, 2 = unimproved, 3 = open defecation (WS11 = 95), 9 = missing.
.nga_toilet_type <- function(hh) {
  ws11 <- .as_code(hh$WS11)

  dplyr::case_when(
    ws11 %in% c(11, 12, 13, 18, 21, 22, 31) ~ 1L,
    ws11 %in% 95 ~ 3L,
    ws11 %in% 99 ~ 9L,
    TRUE ~ 2L
  )
}


# `sharedToilet` from define/MICS6 - 10 - WS.sps:
#   recode WS17 (1 thru 5 = 1) (97, 98, 99 = 9) (sysmis = 0) (else = 2)
#   if (WS16 = 2) sharedToilet = 3
# 0 = not shared, 1 = 5 households or less, 2 = more than 5, 3 = public
# facility, 9 = DK/missing. WS17 is left missing when the facility is not
# shared, which is what produces the 0 category.
.nga_shared_toilet <- function(hh) {
  ws17 <- .as_code(hh$WS17)
  ws16 <- .as_code(hh$WS16)

  shared <- dplyr::case_when(
    is.na(ws17) ~ 0L,
    ws17 >= 1 & ws17 <= 5 ~ 1L,
    ws17 %in% c(97, 98, 99) ~ 9L,
    TRUE ~ 2L
  )

  # WS16 = 2 ("shared with general public") overrides, applied after the
  # recode exactly as the SPSS does.
  shared[ws16 %in% 2] <- 3L
  shared
}


#' Use of at least basic drinking water services (Nigeria)
#'
#' Percentage of the household population using an improved drinking water
#' source, either on premises or with a round trip collection time of 30
#' minutes or less.
#'
#' Translated from \code{MICS6 - 10 - WS.1.2.sps} (\code{INDWS2}) and the
#' \code{drinkingWater} recode in \code{define/MICS6 - 10 - WS.sps}
#' (MICS indicator WS.1, SDG indicator 6.1.1).
#'
#' The denominator is the household \emph{population}: the household weight is
#' multiplied by household size (\code{HH48}), as the source table does.
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021)
#' household recode variable names.
#'
#' @param hh A household data.frame (hh.sav) from Nigeria MICS6.
#'
#' @return A data.frame with columns
#'   \code{cluster}, \code{householdID}, \code{weight} (population weight),
#'   \code{strata}, \code{value} (1 = at least basic service), and
#'   \code{indicator} (\code{"basicwater"}).
#' @export
process_BASICWATER <- function(hh) {
  hh <- .nga_hh_population_denominator(
    hh, "process_BASICWATER", c("WS1", "WS2", "WS3", "WS4")
  )

  # `time`: recode WS4 (0 = 2) (1 thru 30 = 2) (31 thru 990 = 3)
  #         (998, 999 = 9), then water on premises overrides to 1.
  ws4 <- .as_code(hh$WS4)

  time <- dplyr::case_when(
    ws4 %in% c(998, 999) ~ 9L,
    ws4 >= 0 & ws4 <= 30 ~ 2L,
    ws4 >= 31 & ws4 <= 990 ~ 3L,
    TRUE ~ NA_integer_
  )

  # if (any(WS1, 11, 12) or any(WS2, 11, 12) or any(WS3, 1, 2)) time = 1.
  on_premises <- .as_code(hh$WS1) %in% c(11, 12) |
    .as_code(hh$WS2) %in% c(11, 12) |
    .as_code(hh$WS3) %in% c(1, 2)
  time[on_premises] <- 1L

  # INDWS2 = 100 if (drinkingWater = 1 and time <= 2).
  basic <- .nga_improved_water(hh) & !is.na(time) & time <= 2
  hh$value <- as.integer(basic)

  .indicator_result(hh, "basicwater")
}


#' Use of at least basic sanitation services (Nigeria)
#'
#' Percentage of the household population using an improved sanitation facility
#' that is not shared with other households.
#'
#' Translated from \code{MICS6 - 10 - WS.3.2.sps} and the \code{toiletType} /
#' \code{sharedToilet} recodes in \code{define/MICS6 - 10 - WS.sps}
#' (MICS indicator WS.3, SDG indicator 6.2.1a). In the source table the
#' indicator is the "Not shared" column within "Users of improved sanitation
#' facilities", i.e. \code{toiletType = 1} and \code{sharedToilet = 0}.
#'
#' The denominator is the household \emph{population}: the household weight is
#' multiplied by household size (\code{HH48}), as the source table does.
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021)
#' household recode variable names.
#'
#' @param hh A household data.frame (hh.sav) from Nigeria MICS6.
#'
#' @return A data.frame with columns
#'   \code{cluster}, \code{householdID}, \code{weight} (population weight),
#'   \code{strata}, \code{value} (1 = at least basic service), and
#'   \code{indicator} (\code{"basicsan"}).
#' @export
process_BASICSAN <- function(hh) {
  hh <- .nga_hh_population_denominator(
    hh, "process_BASICSAN", c("WS11", "WS16", "WS17")
  )

  hh$value <- as.integer(
    .nga_toilet_type(hh) == 1L & .nga_shared_toilet(hh) == 0L
  )

  .indicator_result(hh, "basicsan")
}


#' Open defecation (Nigeria)
#'
#' Percentage of the household population practising open defecation, i.e.
#' living in households reporting no sanitation facility (bush, field).
#'
#' Translated from \code{MICS6 - 10 - WS.3.1.sps} and the \code{toiletType}
#' recode in \code{define/MICS6 - 10 - WS.sps}, where category 3 is
#' "Open defecation (no facility, bush, field)" and corresponds to
#' \code{WS11 = 95} (MICS indicator WS.4, SDG indicator 6.2.1).
#'
#' The denominator is the household \emph{population}: the household weight is
#' multiplied by household size (\code{HH48}), as the source table does.
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021)
#' household recode variable names.
#'
#' @param hh A household data.frame (hh.sav) from Nigeria MICS6.
#'
#' @return A data.frame with columns
#'   \code{cluster}, \code{householdID}, \code{weight} (population weight),
#'   \code{strata}, \code{value} (1 = practises open defecation), and
#'   \code{indicator} (\code{"opendef"}).
#' @export
process_OPENDEF <- function(hh) {
  hh <- .nga_hh_population_denominator(hh, "process_OPENDEF", "WS11")

  hh$value <- as.integer(.nga_toilet_type(hh) == 3L)

  .indicator_result(hh, "opendef")
}


# ---------------------------------------------------------------------------
# Priority two
# ---------------------------------------------------------------------------

#' Drinking water available in sufficient quantities (Nigeria)
#'
#' Percentage of the household population with drinking water available in
#' sufficient quantities, i.e. in households reporting no occasion in the last
#' month when water was insufficient.
#'
#' Translated from \code{MICS6 - 10 - WS.1.5.sps} (\code{sufficentwater}), where
#' the indicator is \code{WS7 = 2} - a "No" to whether there was any time in
#' the last month without sufficient water.
#'
#' The denominator is the household \emph{population}: the household weight is
#' multiplied by household size (\code{HH48}), as the source table does.
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021)
#' household recode variable names.
#'
#' @param hh A household data.frame (hh.sav) from Nigeria MICS6.
#'
#' @return A data.frame with the standard six indicator columns;
#'   \code{indicator} is \code{"watersuff"}.
#' @export
process_WATERSUFF <- function(hh) {
  hh <- .nga_hh_population_denominator(hh, "process_WATERSUFF", "WS7")

  # sufficentwater = 100 if WS7 = 2.
  hh$value <- as.integer(.as_code(hh$WS7) %in% 2)

  .indicator_result(hh, "watersuff")
}


#' Menstrual hygiene management, women age 15-49 (Nigeria)
#'
#' Percentage of women age 15-49 years who used appropriate menstrual hygiene
#' materials and had a private place to wash and change while at home.
#'
#' Translated from \code{MICS6 - 10 - WS.4.1.sps} (\code{materialprivateplace};
#' MICS indicator WS.10). The numerator requires both a private place
#' (\code{UN17 = 1}) and use of appropriate materials (\code{UN18 = 1}).
#'
#' The denominator is women who completed the interview and answered the
#' menstrual hygiene module, reproducing
#' \code{select if (not(sysmis(UN16)))} - women who did not menstruate in the
#' last 12 months are not asked and are excluded rather than counted as 0.
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021) women's
#' recode variable names.
#'
#' @param wm A women's data.frame (wm.sav) from Nigeria MICS6.
#'
#' @return A data.frame with the standard six indicator columns;
#'   \code{indicator} is \code{"menstrual"}.
#' @export
process_MENSTRUAL <- function(wm) {
  wm <- standardize_columns(wm, "wm")

  required <- c(
    "cluster", "householdID", "weight", "strata", "WM17", "UN16", "UN17", "UN18"
  )
  .require_columns(wm, required, "process_MENSTRUAL")

  # select if (WM17 = 1). select if (not(sysmis(UN16))).
  keep <- .as_code(wm$WM17) %in% 1 & !is.na(.as_code(wm$UN16))
  wm <- wm[keep, , drop = FALSE]

  # materialprivateplace = 100 if (privateplace = 100 and UN18 = 1).
  wm$value <- as.integer(
    .as_code(wm$UN17) %in% 1 & .as_code(wm$UN18) %in% 1
  )

  .indicator_result(wm, "menstrual")
}
