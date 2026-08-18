# R/indicators_nga_nutrition.R
#
# Nigeria-only nutrition indicators, translated from the official MICS6 SPSS
# tabulation syntax (mics.unicef.org/tools, "42 Syntax Files 20220509") against
# the Nigeria MICS6 (2021) recode variable names.
#
# Source files:
#   MICS6 - 07 - TC.7.3.sps  Breastfeeding status  -> EBF
# Shared recodes live in define/MICS6 - 07 - TC.sps.
#
# NOT IMPLEMENTED - stunting (HAZ < -2) and severe wasting (WHZ < -3):
# Nigeria MICS6 (2021) has no anthropometry. None of the eight released recode
# files (hh, hl, wm, bh, fg, ch, fs, mn) contains height, weight, oedema, or
# any z-score/flag variable, and ch.sav carries no AN module - its module
# prefixes are CA, IM, UF, EC, BD, UCD, UB, UCF, BR and ED only. These two
# indicators cannot be derived from this survey and would need another source
# (e.g. NDHS, or the MICS 2016-17 round).


#' Exclusive breastfeeding, infants under 6 months (Nigeria)
#'
#' Percentage of infants age 0-5 months who are exclusively breastfed: children
#' currently breastfeeding who received no other liquid or food the previous
#' day.
#'
#' Translated from \code{MICS6 - 07 - TC.7.3.sps} and the
#' \code{exclusivelyBreastfed} recode in \code{define/MICS6 - 07 - TC.sps}
#' (MICS indicator TC.33).
#'
#' Note that the active SPSS definition requires an explicit "No" (code 2) on
#' every liquid and food item rather than merely the absence of a "Yes"; the
#' commented-out earlier version used \code{not any(1, ...)}. A child with a
#' missing response on any item is therefore not counted as exclusively
#' breastfed. This translation follows the active definition.
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021)
#' children's recode variable names.
#'
#' @param ch A children-under-5 data.frame (ch.sav) from Nigeria MICS6.
#'
#' @return A data.frame with columns
#'   \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#'   \code{value} (1 = exclusively breastfed), and \code{indicator}
#'   (\code{"ebf"}).
#' @export
process_EBF <- function(ch) {
  ch <- standardize_columns(ch, "ch")

  liquids <- c("BD7A", "BD7B", "BD7C", "BD7D", "BD7E", "BD7X")
  foods <- c(
    "BD8A", "BD8B", "BD8C", "BD8D", "BD8E", "BD8F", "BD8G", "BD8H",
    "BD8I", "BD8J", "BD8K", "BD8L", "BD8M", "BD8N", "BD8X"
  )

  required <- c(
    "cluster", "householdID", "weight", "strata",
    "UF17", "CAGE", "BD3", liquids, foods
  )
  .require_columns(ch, required, "process_EBF")

  # Denominator: `do if (cage >= 0 and cage <= 5)` over completed under-5
  # interviews. Nigeria names the age-in-months variable CAGE.
  age <- .as_code(ch$CAGE)
  keep <- .as_code(ch$UF17) %in% 1 & !is.na(age) & age >= 0 & age <= 5
  ch <- ch[keep, , drop = FALSE]

  # exclusivelyBreastfed = 100 if BD3 = 1 and every BD7/BD8 item = 2.
  .all_no <- function(vars) {
    hits <- lapply(vars, function(v) .as_code(ch[[v]]) %in% 2)
    Reduce(`&`, hits)
  }

  ch$value <- as.integer(
    .as_code(ch$BD3) %in% 1 & .all_no(liquids) & .all_no(foods)
  )

  .indicator_result(ch, "ebf")
}


#' Vitamin A supplementation, children age 6-59 months (Nigeria)
#'
#' Percentage of children age 6-59 months who received a vitamin A supplement,
#' from a vaccination card or the mother's report.
#'
#' \strong{This indicator has no standard MICS6 tabulation syntax.} Vitamin A
#' supplementation is not among the 42 standard MICS6 tabulation files; the
#' items used here (\code{IM27B}, \code{IM6Z1*}, \code{IM6Z2*}) are specific to
#' the Nigeria MICS-NICS questionnaire. The definition below therefore applies
#' the \code{vaccdef} macro pattern from
#' \code{MICS6 - 07 - TC - RECODE.sps} - the same card-or-report logic used for
#' every vaccine - to the Nigeria vitamin A items.
#'
#' \strong{FLAG - this measures "ever received", not "received in the last 6
#' months".} The global indicator is normally defined over the six months
#' preceding the survey, but Nigeria's questionnaire asks only \code{IM27B}
#' "Child ever given Vitamin A" and records dose dates on the card. A
#' six-month version would be computable only for card holders and would
#' undercount badly. Confirm against the Nigeria survey findings report before
#' publishing, and adjust \code{dose_day_vars} if only the first dose should
#' count.
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021)
#' children's recode variable names.
#'
#' @param ch A children-under-5 data.frame (ch.sav) from Nigeria MICS6.
#' @param dose_day_vars Character vector of vitamin A dose "day" variables
#'   recorded from the card. Default \code{c("IM6Z1D", "IM6Z2D")}, i.e. a dose
#'   recorded against either card entry counts.
#'
#' @return A data.frame with columns
#'   \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#'   \code{value} (1 = received a vitamin A supplement), and \code{indicator}
#'   (\code{"vita"}).
#' @export
process_VITA <- function(ch, dose_day_vars = c("IM6Z1D", "IM6Z2D")) {
  ch <- standardize_columns(ch, "ch")

  required <- c(
    "cluster", "householdID", "weight", "strata",
    "UF17", "CAGE", "IM27B", dose_day_vars
  )
  .require_columns(ch, required, "process_VITA")

  # Denominator: children age 6-59 months with a completed interview.
  age <- .as_code(ch$CAGE)
  keep <- .as_code(ch$UF17) %in% 1 & !is.na(age) & age >= 6 & age <= 59
  ch <- ch[keep, , drop = FALSE]

  # vaccdef, specialised for vbirth = 0, vnum = 1 (the BCG-style single-dose
  # shape): a day code of 1-31, 44 or 97-99 is "given according to card",
  # 66 is "mother's report marked on card", and vgiven = 1 is mother's report.
  from_card <- Reduce(
    `|`,
    lapply(dose_day_vars, function(v) {
      day <- .as_code(ch[[v]])
      day %in% c(1:31, 44, 66, 97:99)
    })
  )

  from_report <- .as_code(ch$IM27B) %in% 1

  ch$value <- as.integer(from_card | from_report)

  .indicator_result(ch, "vita")
}
