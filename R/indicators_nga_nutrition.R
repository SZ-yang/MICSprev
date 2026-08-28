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


# ---------------------------------------------------------------------------
# Priority two
# ---------------------------------------------------------------------------
#
# NOT IMPLEMENTED - micronutrient powders (MNP), 6-59 months: Nigeria MICS6
# (2021) has no MNP item. A search of ch.sav for micronutrient, powder, sachet
# or sprinkle variables returns nothing, so like stunting and wasting this
# indicator cannot be derived from this survey.


# `solidSemiSoft` from define/MICS6 - 07 - TC.sps:
#   any(1, BD8A ... BD8N, BD8X)
.nga_solid_semi_soft <- function(ch) {
  foods <- c(
    "BD8A", "BD8B", "BD8C", "BD8D", "BD8E", "BD8F", "BD8G", "BD8H",
    "BD8I", "BD8J", "BD8K", "BD8L", "BD8M", "BD8N", "BD8X"
  )
  Reduce(`|`, lapply(foods, function(v) .as_code(ch[[v]]) %in% 1))
}


#' Age-appropriate breastfeeding, children age 0-23 months (Nigeria)
#'
#' Percentage of children age 0-23 months who were appropriately breastfed
#' during the previous day: exclusively breastfed if under 6 months, or
#' breastfed and receiving solid, semi-solid or soft foods if age 6-23 months.
#'
#' Translated from \code{MICS6 - 07 - TC.7.5.sps} (\code{appropAll}; MICS
#' indicator TC.36).
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021)
#' children's recode variable names.
#'
#' @param ch A children-under-5 data.frame (ch.sav) from Nigeria MICS6.
#'
#' @return A data.frame with the standard six indicator columns;
#'   \code{indicator} is \code{"appropbf"}.
#' @export
process_APPROPBF <- function(ch) {
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
  .require_columns(ch, required, "process_APPROPBF")

  age <- .as_code(ch$CAGE)
  keep <- .as_code(ch$UF17) %in% 1 & !is.na(age) & age >= 0 & age <= 23
  ch <- ch[keep, , drop = FALSE]
  age <- age[keep]

  breastfed <- .as_code(ch$BD3) %in% 1

  # approp0_5 = exclusivelyBreastfed (BD3 = 1 and every BD7/BD8 item = 2).
  .all_no <- function(vars) {
    Reduce(`&`, lapply(vars, function(v) .as_code(ch[[v]]) %in% 2))
  }
  exclusive <- breastfed & .all_no(liquids) & .all_no(foods)

  # approp6_23 = 100 if (BD3 = 1 and solidSemiSoft).
  complementary <- breastfed & .nga_solid_semi_soft(ch)

  ch$value <- as.integer(
    (age <= 5 & exclusive) | (age >= 6 & complementary)
  )

  .indicator_result(ch, "appropbf")
}


#' Minimum acceptable diet, children age 6-23 months (Nigeria)
#'
#' Percentage of children age 6-23 months who received a minimum acceptable
#' diet during the previous day.
#'
#' Translated from \code{MICS6 - 07 - TC.7.7.sps} (\code{minimumAcceptableDiet};
#' MICS indicator TC.39). The definition differs by breastfeeding status:
#'
#' \itemize{
#'   \item \strong{Breastfed} children need the minimum dietary diversity (at
#'     least 5 of 8 food groups, where breastmilk counts as a group) and the
#'     minimum meal frequency (solid, semi-solid or soft foods at least twice
#'     at age 6-8 months, three times at 9-23 months).
#'   \item \strong{Non-breastfed} children need the adjusted dietary diversity
#'     (at least 4 of 6 groups, dairy excluded), the minimum meal frequency
#'     (milk feeds plus food at least four times), and at least two milk feeds.
#' }
#'
#' The count variables \code{BD7D1}, \code{BD7E1}, \code{BD8A1} and \code{BD9}
#' are recoded so that system-missing, 8 and 9 become 0, and \code{BD3} is
#' recoded so that system-missing, 8 and 9 become 2 (not breastfed), exactly as
#' the source syntax does before computing anything.
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021)
#' children's recode variable names.
#'
#' @param ch A children-under-5 data.frame (ch.sav) from Nigeria MICS6.
#'
#' @return A data.frame with the standard six indicator columns;
#'   \code{indicator} is \code{"minacceptdiet"}.
#' @export
process_MINACCEPTDIET <- function(ch) {
  ch <- standardize_columns(ch, "ch")

  groups <- c(
    "BD7D", "BD7E", "BD8A", "BD8B", "BD8C", "BD8D", "BD8E", "BD8F", "BD8G",
    "BD8H", "BD8I", "BD8J", "BD8K", "BD8L", "BD8M", "BD8N"
  )
  counts <- c("BD7D1", "BD7E1", "BD8A1", "BD9")

  required <- c(
    "cluster", "householdID", "weight", "strata",
    "UF17", "CAGE", "BD3", groups, counts
  )
  .require_columns(ch, required, "process_MINACCEPTDIET")

  age <- .as_code(ch$CAGE)
  keep <- .as_code(ch$UF17) %in% 1 & !is.na(age) & age >= 6 & age <= 23
  ch <- ch[keep, , drop = FALSE]
  age <- age[keep]

  # recode BD7D1 BD7E1 BD8A1 BD9 (sysmis, 8, 9 = 0).
  .count <- function(v) {
    x <- .as_code(ch[[v]])
    x[is.na(x) | x %in% c(8, 9)] <- 0
    x
  }
  bd7d1 <- .count("BD7D1")
  bd7e1 <- .count("BD7E1")
  bd8a1 <- .count("BD8A1")
  bd9 <- .count("BD9")

  # recode BD3 (sysmis, 8, 9 = 2), i.e. treat unknown as not breastfed.
  breastfed <- .as_code(ch$BD3) %in% 1

  .yes <- function(vars) {
    Reduce(`|`, lapply(vars, function(v) .as_code(ch[[v]]) %in% 1))
  }

  # The 8 food groups.
  grains <- .yes(c("BD8B", "BD8C", "BD8E"))
  legumes <- .yes("BD8M")
  dairy <- .yes(c("BD7D", "BD7E", "BD8A", "BD8N"))
  flesh <- .yes(c("BD8I", "BD8J", "BD8L"))
  eggs <- .yes("BD8K")
  vit_a <- .yes(c("BD8D", "BD8F", "BD8G"))
  other <- .yes("BD8H")

  diversity <- rowSums(cbind(
    breastfed, grains, legumes, dairy, flesh, eggs, vit_a, other
  )) >= 5

  # Adjusted diversity for non-breastfed children: 6 groups, dairy dropped.
  diversity_adj <- rowSums(cbind(
    grains, legumes, flesh, eggs, vit_a, other
  )) >= 4

  # Minimum meal frequency.
  meal_frequency <- ifelse(
    breastfed,
    ifelse(age <= 8, bd9 >= 2, bd9 >= 3),
    bd7d1 + bd7e1 + bd9 >= 4
  )

  # Non-breastfed children additionally need at least 2 milk feeds.
  two_milk_feeds <- bd7d1 + bd7e1 + bd8a1 >= 2

  ch$value <- as.integer(ifelse(
    breastfed,
    diversity & meal_frequency,
    diversity_adj & meal_frequency & two_milk_feeds
  ))

  .indicator_result(ch, "minacceptdiet")
}


#' Severe food insecurity in the last 12 months (Nigeria)
#'
#' Percentage of the household population living in households that experienced
#' severe food insecurity in the 12 months preceding the survey.
#'
#' \strong{This indicator has no standard MICS6 tabulation syntax.} Food
#' security is not among the 42 standard MICS6 tabulation files. Nigeria
#' fielded the full eight-item Food Insecurity Experience Scale (FIES) as
#' \code{FE1}-\code{FE8}, and this function applies the standard FAO raw-score
#' cut-off: a household is severely food insecure when it answers "yes" to at
#' least seven of the eight items (SDG indicator 2.1.2, severe component).
#'
#' \strong{FLAG - raw score, not Rasch-adjusted, and the difference is large.}
#' The official SDG 2.1.2 figure is produced by fitting a Rasch model and
#' equating the national scale to the global reference scale; this function
#' does not do that. Item endorsement in the Nigeria data is very high (79
#' percent affirm \code{FE1}, falling to 33 percent for \code{FE8}), so the
#' raw-score result here is roughly 46 percent - about double typical
#' Rasch-equated published estimates for Nigeria. Treat it as a direct
#' tabulation of the survey items, not as a reproduction of SDG 2.1.2, and do
#' not compare it with FAO figures without equating first. \code{threshold} is
#' exposed so the moderate-or-severe cut-off (4 of 8) can be produced instead.
#'
#' The \code{FE*B} follow-ups ("Did this happen in the last 1 month?") are
#' deliberately not used: the base items carry the 12-month reference period
#' this indicator requires.
#'
#' The denominator is the household \emph{population}: the household weight is
#' multiplied by household size (\code{HH48}).
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021)
#' household recode variable names.
#'
#' @param hh A household data.frame (hh.sav) from Nigeria MICS6.
#' @param threshold Minimum number of affirmative FIES items. Default 7
#'   (severe); use 4 for moderate or severe.
#'
#' @return A data.frame with the standard six indicator columns;
#'   \code{indicator} is \code{"foodinsec"}.
#' @export
process_FOODINSEC <- function(hh, threshold = 7) {
  fies <- paste0("FE", 1:8)

  hh <- .nga_hh_population_denominator(hh, "process_FOODINSEC", fies)

  # Affirmative answers only; "No", DK (8) and no response (9) score 0.
  score <- rowSums(
    do.call(cbind, lapply(fies, function(v) .as_code(hh[[v]]) %in% 1))
  )

  hh$value <- as.integer(score >= threshold)

  .indicator_result(hh, "foodinsec")
}
