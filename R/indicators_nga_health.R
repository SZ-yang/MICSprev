# R/indicators_nga_health.R
#
# Nigeria-only health indicators, translated from the official MICS6 SPSS
# tabulation syntax (mics.unicef.org/tools, "42 Syntax Files 20220509") against
# the Nigeria MICS6 (2021) recode variable names.
#
# Source files:
#   MICS6 - 06 - TM.4.1.sps   Antenatal care coverage          -> ANC1
#   MICS6 - 06 - TM.6.2.sps   Assistance during delivery       -> SBA
#   MICS6 - 06 - TM.8.2.sps   Post-natal health checks newborn -> PNCNB
#   MICS6 - 06 - TM.8.7.sps   Post-natal health checks mothers -> PNCMOM
#   MICS6 - 07 - TC.1.1.sps   Vaccinations, 12-23 months       -> PENTA1
# Shared recodes live in define/MICS6 - 06 - TM.sps and
# MICS6 - 07 - TC - RECODE.sps.


# All four women's-questionnaire indicators below share the same denominator:
# women who completed the interview (WM17 = 1) and had a live birth in the two
# years preceding the survey (CM17 = 1). Source: every TM.*.sps listed above
# opens with `select if (WM17 = 1).` followed by `select if (CM17 = 1).`
.nga_recent_birth_denominator <- function(wm, fn, extra_required) {
  wm <- standardize_columns(wm, "wm")

  required <- c(
    "cluster", "householdID", "weight", "strata",
    "WM17", "CM17", extra_required
  )
  .require_columns(wm, required, fn)

  # `%in%` rather than `==` so that missing codes drop out instead of
  # producing NA rows.
  wm[.as_code(wm$WM17) %in% 1 & .as_code(wm$CM17) %in% 1, , drop = FALSE]
}


#' Antenatal care coverage: at least one visit from skilled health personnel (Nigeria)
#'
#' Percentage of women age 15-49 years with a live birth in the last 2 years
#' who were attended at least once during the pregnancy of their most recent
#' live birth by skilled health personnel (medical doctor, nurse/midwife, or
#' other qualified provider).
#'
#' Translated from \code{MICS6 - 06 - TM.4.1.sps} and the \code{anc} /
#' \code{skilledPersonnel} recodes in \code{define/MICS6 - 06 - TM.sps}
#' (MICS indicator TM.5a).
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021) women's
#' recode variable names.
#'
#' @param wm A women's data.frame (wm.sav) from Nigeria MICS6.
#'
#' @return A data.frame with columns
#'   \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#'   \code{value} (1 = attended by skilled personnel), and \code{indicator}
#'   (\code{"anc1"}).
#' @export
process_ANC1 <- function(wm) {
  mr_vars <- c("MN3A", "MN3B", "MN3C", "MN3F", "MN3G", "MN3X", "MN3NR")
  wm <- .nga_recent_birth_denominator(wm, "process_ANC1", mr_vars)

  # `anc` from define/MICS6 - 06 - TM.sps: sequential `if (anc = 0 and ...)`
  # assignments, i.e. first match wins.
  anc <- dplyr::case_when(
    .mr_yes(wm$MN3A, "A") ~ 11L,   # Medical doctor
    .mr_yes(wm$MN3B, "B") ~ 12L,   # Nurse / Midwife
    .mr_yes(wm$MN3C, "C") ~ 13L,   # Other qualified
    .mr_yes(wm$MN3F, "F") ~ 21L,   # Traditional birth attendant
    .mr_yes(wm$MN3G, "G") ~ 22L,   # Community health worker
    .mr_yes(wm$MN3X, "X") | .mr_yes(wm$MN3NR, "?") ~ 97L,  # Other/Missing
    TRUE ~ 98L                     # No antenatal care
  )

  # skilledPersonnel = 100 if any(anc, 11, 12, 13).
  wm$value <- as.integer(anc %in% c(11L, 12L, 13L))

  .indicator_result(wm, "anc1")
}


#' Skilled attendance at birth (Nigeria)
#'
#' Percentage of women age 15-49 years with a live birth in the last 2 years
#' whose most recent live birth was assisted by skilled health personnel
#' (medical doctor, nurse/midwife, or other qualified provider).
#'
#' Translated from \code{MICS6 - 06 - TM.6.2.sps} and the
#' \code{personAtDelivery} / \code{skilledAttendant} recodes in
#' \code{define/MICS6 - 06 - TM.sps} (MICS indicator TM.9a).
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021) women's
#' recode variable names.
#'
#' @param wm A women's data.frame (wm.sav) from Nigeria MICS6.
#'
#' @return A data.frame with columns
#'   \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#'   \code{value} (1 = skilled attendant at delivery), and \code{indicator}
#'   (\code{"sba"}).
#' @export
process_SBA <- function(wm) {
  mr_vars <- c(
    "MN19A", "MN19B", "MN19C", "MN19F", "MN19G", "MN19H",
    "MN19X", "MN19Y", "MN19NR"
  )
  wm <- .nga_recent_birth_denominator(wm, "process_SBA", mr_vars)

  # `personAtDelivery` from define/MICS6 - 06 - TM.sps: sequential
  # `if (personAtDelivery = 0 and ...)` assignments, i.e. first match wins.
  # The trailing `recode personAtDelivery (0 = 98)` folds unassigned cases into
  # "No attendant", which the TRUE branch below reproduces.
  person <- dplyr::case_when(
    .mr_yes(wm$MN19A, "A") ~ 11L,  # Skilled: Medical doctor
    .mr_yes(wm$MN19B, "B") ~ 12L,  # Skilled: Nurse / Midwife
    .mr_yes(wm$MN19C, "C") ~ 13L,  # Skilled: Other qualified
    .mr_yes(wm$MN19F, "F") ~ 14L,  # Traditional birth attendant
    .mr_yes(wm$MN19G, "G") ~ 15L,  # Community / village health worker
    .mr_yes(wm$MN19H, "H") ~ 16L,  # Relative / Friend
    .mr_yes(wm$MN19X, "X") | .mr_yes(wm$MN19NR, "?") ~ 97L,  # Other/Missing
    TRUE ~ 98L                     # No attendant (incl. MN19Y = "Y")
  )

  # skilledAttendant = 100 if any(personAtDelivery, 11, 12, 13).
  wm$value <- as.integer(person %in% c(11L, 12L, 13L))

  .indicator_result(wm, "sba")
}


# Shared postnatal-care logic for TM.8.2 (newborns) and TM.8.7 (mothers). The
# two tables are structurally identical; only the variable names differ.
#
# check_vars  - the two "health check before leaving facility / after delivery"
#               variables (PN4/PN8 for newborns, PN5/PN9 for mothers)
# visit_vars  - the three PNC-visit gate variables (PN6/PN10/PN11 or
#               PN17/PN19/PN20). Kept in the required set because the SPSS
#               `pncVisitC`/`pncVisitM` recode reads them.
# time_unit   - PN13U / PN22U
# time_number - PN13N / PN22N
# provider    - the five "person checking health" multi-response columns
.nga_pnc_value <- function(wm, check_vars, time_unit, time_number, provider) {
  # healthCheck = 100 if (PN4 = 1 or PN8 = 1)  [TM.8.2]
  #               100 if (PN5 = 1 or PN9 = 1)  [TM.8.7]
  health_check <-
    .as_code(wm[[check_vars[1]]]) %in% 1 | .as_code(wm[[check_vars[2]]]) %in% 1

  # pncmedC / pncmedM: person checking health, excluding relative/friend and
  # other (letters A, B, C, F, G).
  letters_ <- c("A", "B", "C", "F", "G")
  pncmed <- Reduce(
    `|`,
    Map(function(v, l) .mr_yes(wm[[v]], l), provider, letters_)
  )

  unit <- .as_code(wm[[time_unit]])
  number <- .as_code(wm[[time_number]])

  # pncVisitC / pncVisitM categories 1-3 are the checks that happened within
  # two days of the birth: 1 = same day, 2 = 1 day following birth,
  # 3 = 2 days following birth. All three require pncmed = 1.
  within_two_days <- pncmed &
    (unit %in% 1 | (unit %in% 2 & number %in% c(1, 2)))
  within_two_days[is.na(within_two_days)] <- FALSE

  # pncCheck = 100 if (healthCheck = 100 or pncVisit in 1, 2, 3).
  as.integer(health_check | within_two_days)
}


#' Post-natal health check for the newborn within 2 days of birth (Nigeria)
#'
#' Percentage of women age 15-49 years with a live birth in the last 2 years
#' whose most recent live-born child received a health check while in facility
#' or at home following delivery, or a post-natal care visit within 2 days of
#' birth.
#'
#' Translated from \code{MICS6 - 06 - TM.8.2.sps} and the \code{pncmedC} /
#' \code{pncVisitC} recodes in \code{define/MICS6 - 06 - TM.sps}
#' (MICS indicator TM.10a).
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021) women's
#' recode variable names.
#'
#' @param wm A women's data.frame (wm.sav) from Nigeria MICS6.
#'
#' @return A data.frame with columns
#'   \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#'   \code{value} (1 = newborn received a post-natal health check), and
#'   \code{indicator} (\code{"pncnb"}).
#' @export
process_PNCNB <- function(wm) {
  provider <- c("PN14A", "PN14B", "PN14C", "PN14F", "PN14G")
  extra <- c("PN4", "PN8", "PN6", "PN10", "PN11", "PN13U", "PN13N", provider)

  wm <- .nga_recent_birth_denominator(wm, "process_PNCNB", extra)

  wm$value <- .nga_pnc_value(
    wm,
    check_vars  = c("PN4", "PN8"),
    time_unit   = "PN13U",
    time_number = "PN13N",
    provider    = provider
  )

  .indicator_result(wm, "pncnb")
}


#' Post-natal health check for the mother within 2 days of birth (Nigeria)
#'
#' Percentage of women age 15-49 years with a live birth in the last 2 years
#' who received a health check while in facility or at home following delivery,
#' or a post-natal care visit within 2 days of birth.
#'
#' Translated from \code{MICS6 - 06 - TM.8.7.sps} and the \code{pncmedM} /
#' \code{pncVisitM} recodes in \code{define/MICS6 - 06 - TM.sps}
#' (MICS indicator TM.10b).
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021) women's
#' recode variable names.
#'
#' @param wm A women's data.frame (wm.sav) from Nigeria MICS6.
#'
#' @return A data.frame with columns
#'   \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#'   \code{value} (1 = mother received a post-natal health check), and
#'   \code{indicator} (\code{"pncmom"}).
#' @export
process_PNCMOM <- function(wm) {
  provider <- c("PN23A", "PN23B", "PN23C", "PN23F", "PN23G")
  extra <- c("PN5", "PN9", "PN17", "PN19", "PN20", "PN22U", "PN22N", provider)

  wm <- .nga_recent_birth_denominator(wm, "process_PNCMOM", extra)

  wm$value <- .nga_pnc_value(
    wm,
    check_vars  = c("PN5", "PN9"),
    time_unit   = "PN22U",
    time_number = "PN22N",
    provider    = provider
  )

  .indicator_result(wm, "pncmom")
}


#' Penta-1 (DTP-HepB-Hib first dose) coverage, children 12-23 months (Nigeria)
#'
#' Percentage of children age 12-23 months who received the first dose of the
#' pentavalent (DTP-HepB-Hib) vaccine by any time before the survey, from a
#' vaccination card or the mother's report.
#'
#' Translated from \code{MICS6 - 07 - TC.1.1 (12 - 23 mo).sps} and the
#' \code{vaccdef} macro in \code{MICS6 - 07 - TC - RECODE.sps}, invoked as
#' \code{vaccdef vacc=penta1 vd=IM6PENTA1D vgiven=IM20 vbirth=0 vtot=IM21
#' vnum=1}. A child counts as vaccinated when \code{vaccdef} yields 1
#' ("Vaccination card") or 2 ("Mother's report"); 3 ("Not vaccinated") and 8
#' ("DK") count as not vaccinated and stay in the denominator.
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021)
#' children's recode variable names.
#'
#' @param ch A children-under-5 data.frame (ch.sav) from Nigeria MICS6.
#'
#' @return A data.frame with columns
#'   \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#'   \code{value} (1 = received Penta-1), and \code{indicator}
#'   (\code{"penta1"}).
#' @export
process_PENTA1 <- function(ch) {
  ch <- standardize_columns(ch, "ch")

  required <- c(
    "cluster", "householdID", "weight", "strata",
    "UF17", "CAGE_6", "IM6PENTA1D", "IM20", "IM21"
  )
  .require_columns(ch, required, "process_PENTA1")

  # Denominator: completed under-5 interviews (UF17 = 1) for children age
  # 12-23 months. CAGE_6 = 3 is the "12-23" category, matching
  # `select if (cage >= 12 and cage <= 23)` in TC.1.1 (part I).
  ch <- ch[.as_code(ch$UF17) %in% 1 & .as_code(ch$CAGE_6) %in% 3, , drop = FALSE]

  day <- .as_code(ch$IM6PENTA1D)
  given <- .as_code(ch$IM20)

  # vaccdef, specialised for vbirth = 0 and vnum = 1:
  #   vacc = 1 if vd in 1..31, 44, or 97..99      (card)
  #   vacc = 2 if vd = 66                         (report, marked on card)
  #   vacc = 2 if vgiven = 1                      (mother's report; the
  #              `vbirth <> 1 & vnum = 1` branch, both true here)
  #   vacc = 8 if vgiven in 8, 9                  (DK)
  # The remaining `vaccdef` branches are unreachable at vbirth = 0, vnum = 1.
  from_card <- day %in% c(1:31, 44, 97:99)
  from_report <- day %in% 66 | given %in% 1

  ch$value <- as.integer(from_card | from_report)

  .indicator_result(ch, "penta1")
}
