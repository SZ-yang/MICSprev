# R/indicators_nga_childprotection.R
#
# Nigeria-only child protection indicators, translated from the official MICS6
# SPSS tabulation syntax (mics.unicef.org/tools, "42 Syntax Files 20220509")
# against the Nigeria MICS6 (2021) recode variable names.
#
# Source files:
#   MICS6 - 09 - PR.1.1.sps   Birth registration          -> BIRTHREG
#   MICS6 - 09 - PR.4.1W.sps  Child marriage (women)      -> CHILDMARRIAGE
#   MICS6 - 09 - PR.3.3.sps   Child labour                -> CHILDLABOUR
#   MICS6 - 09 - PR.5.3.sps   FGM among girls             -> FGMDAUGHTER
# Shared recodes live in define/MICS6 - 09 - PR.sps.


#' Birth registration, children under age 5 (Nigeria)
#'
#' Percentage of children under age 5 whose births are registered with a civil
#' authority: those with a birth certificate seen by the interviewer
#' (\code{BR1 = 1}), those reported to have a certificate that was not seen
#' (\code{BR1 = 2}), and those without a certificate but reported as registered
#' (\code{BR2 = 1}).
#'
#' Translated from \code{MICS6 - 09 - PR.1.1.sps} (\code{birthRegistrated};
#' MICS indicator PR.1, SDG indicator 16.9.1).
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021)
#' children's recode variable names.
#'
#' @param ch A children-under-5 data.frame (ch.sav) from Nigeria MICS6.
#'
#' @return A data.frame with columns
#'   \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#'   \code{value} (1 = birth registered), and \code{indicator}
#'   (\code{"birthreg"}).
#' @export
process_BIRTHREG <- function(ch) {
  ch <- standardize_columns(ch, "ch")

  required <- c(
    "cluster", "householdID", "weight", "strata", "UF17", "BR1", "BR2"
  )
  .require_columns(ch, required, "process_BIRTHREG")

  # Denominator: all children under 5 with a completed interview.
  ch <- ch[.as_code(ch$UF17) %in% 1, , drop = FALSE]

  # birthRegistrated = 100 if (BR1 = 1 or BR1 = 2 or BR2 = 1).
  ch$value <- as.integer(
    .as_code(ch$BR1) %in% c(1, 2) | .as_code(ch$BR2) %in% 1
  )

  .indicator_result(ch, "birthreg")
}


#' Women age 20-49 first married or in union before age 18 (Nigeria)
#'
#' Percentage of women age 20-49 years who were first married or in union
#' before age 18.
#'
#' Translated from \code{MICS6 - 09 - PR.4.1W.sps} (\code{before18} within
#' \code{do if (WB4 >= 20)}; MICS indicator PR.5, SDG indicator 5.3.1).
#'
#' Women who have never married have a missing \code{WAGEM}. In SPSS an
#' \code{if} whose condition is missing does not execute, so \code{before18}
#' keeps its initialised value of 0; the translation reproduces this by
#' treating missing \code{WAGEM} as not-married-before-18 rather than as a
#' missing outcome.
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021) women's
#' recode variable names.
#'
#' @param wm A women's data.frame (wm.sav) from Nigeria MICS6.
#'
#' @return A data.frame with columns
#'   \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#'   \code{value} (1 = first married/in union before age 18), and
#'   \code{indicator} (\code{"childmarriage"}).
#' @export
process_CHILDMARRIAGE <- function(wm) {
  wm <- standardize_columns(wm, "wm")

  required <- c(
    "cluster", "householdID", "weight", "strata", "WM17", "WB4", "WAGEM"
  )
  .require_columns(wm, required, "process_CHILDMARRIAGE")

  # Denominator: completed interviews (WM17 = 1) for women age 20-49
  # (`do if (WB4 >= 20)`). wm.sav covers women 15-49, so WB4 >= 20 is the
  # full 20-49 range.
  keep <- .as_code(wm$WM17) %in% 1 & .as_code(wm$WB4) >= 20
  keep[is.na(keep)] <- FALSE
  wm <- wm[keep, , drop = FALSE]

  age_at_marriage <- .as_code(wm$WAGEM)

  # before18 = 100 if (WAGEM < 18).
  # FLAG: the Nigeria recode contains 14 negative WAGEM values (an artefact of
  # the WDOM - WDOB derivation). PR.4.1W tests `WAGEM < 18` with no lower
  # bound, so these count as married before 18. Translated literally; the
  # effect is 14 of ~28,000 women in the denominator.
  married_before_18 <- age_at_marriage < 18
  married_before_18[is.na(married_before_18)] <- FALSE

  wm$value <- as.integer(married_before_18)

  .indicator_result(wm, "childmarriage")
}


# Shared components of the two child labour definitions, per
# define/MICS6 - 09 - PR.sps. Returns the filtered frame alongside the three
# logical components so each indicator can combine them as its table does.
.nga_child_labour_components <- function(fs, fn) {
  fs <- standardize_columns(fs, "fs")

  econ_vars <- c("CL1A", "CL1B", "CL1C", "CL1X")
  chore_vars <- c(
    "CL7", "CL9",
    "CL11A", "CL11B", "CL11C", "CL11D", "CL11E", "CL11F", "CL11X"
  )
  hazard_vars <- c(
    "CL4", "CL5", "CL6A", "CL6B", "CL6C", "CL6D", "CL6E", "CL6X"
  )

  required <- c(
    "cluster", "householdID", "weight", "strata",
    "FS17", "CB3", "CL3", "CL8", "CL10", "CL13",
    econ_vars, chore_vars, hazard_vars
  )
  .require_columns(fs, required, fn)

  # Denominator: children age 5-17 with a completed interview.
  fs <- fs[.as_code(fs$FS17) %in% 1, , drop = FALSE]

  age <- .as_code(fs$CB3)

  # SPSS `any(1, CL1A, ...)` is true when any listed variable equals 1.
  .any_yes <- function(vars) {
    Reduce(`|`, lapply(vars, function(v) .as_code(fs[[v]]) %in% 1))
  }

  economic_activity <- .any_yes(econ_vars)
  hh_chores <- .any_yes(chore_vars)
  econ_hours <- .as_code(fs$CL3)

  # eaMore: age-specific thresholds, each capped at CL3 < 97 to drop the
  # 97/98/99 DK/missing codes.
  ea_more <- dplyr::case_when(
    age <= 11           ~ economic_activity & econ_hours >= 1  & econ_hours < 97,
    age %in% c(12L:14L) ~ economic_activity & econ_hours >= 14 & econ_hours < 97,
    age %in% c(15L:17L) ~ economic_activity & econ_hours >= 43 & econ_hours < 97,
    TRUE ~ FALSE
  )
  ea_more[is.na(ea_more)] <- FALSE

  # hhcMore: only computed for age 5-14 (`do if (CB3 <= 14)`); children age
  # 15-17 contribute 0 via `sum(0, hhc21less)`. The 97/98/99 codes are recoded
  # to 0 before the >= 21 comparison.
  .chore_hours <- function(v) {
    x <- .as_code(fs[[v]])
    x[x %in% c(97, 98, 99)] <- 0
    x
  }

  # SPSS sum() ignores missing operands.
  chore_hours <- rowSums(
    cbind(.chore_hours("CL8"), .chore_hours("CL10"), .chore_hours("CL13")),
    na.rm = TRUE
  )

  hhc_more <- age <= 14 & hh_chores & chore_hours >= 21
  hhc_more[is.na(hhc_more)] <- FALSE

  # hazardConditions = any(1, CL4, CL5, CL6A ... CL6X).
  hazard <- .any_yes(hazard_vars)

  list(fs = fs, ea_more = ea_more, hhc_more = hhc_more, hazard = hazard)
}


#' Child labour, children age 5-17 (Nigeria)
#'
#' Percentage of children age 5-17 years engaged in child labour during the
#' previous week: economic activity above the age-specific hours threshold
#' (1 hour for age 5-11, 14 hours for age 12-14, 43 hours for age 15-17), or
#' household chores of 21 hours or more for children age 5-14.
#'
#' Translated from \code{MICS6 - 09 - PR.3.3.sps} and the
#' \code{economicActivity} / \code{eaMore} / \code{hhcMore} / \code{childLabor}
#' recodes in \code{define/MICS6 - 09 - PR.sps} (MICS indicator PR.3, SDG
#' indicator 8.7.1).
#'
#' Note that this is the SDG definition, which \strong{excludes} hazardous
#' working conditions: \code{define/MICS6 - 09 - PR.sps} computes
#' \code{childLabor = maximum(0, eaMore, hhcMore)} with the comment
#' "hazardConditions removed from the calculation of child labor". The variant
#' that includes hazardous conditions is a separate indicator
#' (\code{includeHaz} in the same file).
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021)
#' children age 5-17 recode variable names.
#'
#' @param fs A children age 5-17 data.frame (fs.sav) from Nigeria MICS6.
#'
#' @return A data.frame with columns
#'   \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#'   \code{value} (1 = engaged in child labour), and \code{indicator}
#'   (\code{"childlabour"}).
#' @export
process_CHILDLABOUR <- function(fs) {
  parts <- .nga_child_labour_components(fs, "process_CHILDLABOUR")

  # childLabor = maximum(0, eaMore, hhcMore).
  parts$fs$value <- as.integer(parts$ea_more | parts$hhc_more)

  .indicator_result(parts$fs, "childlabour")
}


#' Child labour including hazardous conditions, children age 5-17 (Nigeria)
#'
#' Percentage of children age 5-17 years engaged in economic activities or
#' household chores above the age-specific thresholds, \emph{or} working under
#' hazardous conditions.
#'
#' Translated from the \code{includeHaz} recode in
#' \code{define/MICS6 - 09 - PR.sps}:
#' \code{maximum(0, eaMore, hhcMore, hazardConditions)}, where
#' \code{hazardConditions = any(1, CL4, CL5, CL6A, CL6B, CL6C, CL6D, CL6E,
#' CL6X)}.
#'
#' This is the wider of the two child labour definitions. The narrower one,
#' \code{\link{process_CHILDLABOUR}()}, drops hazardous conditions and is the
#' basis for SDG indicator 8.7.1; see the note in that function.
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021)
#' children age 5-17 recode variable names.
#'
#' @param fs A children age 5-17 data.frame (fs.sav) from Nigeria MICS6.
#'
#' @return A data.frame with columns
#'   \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#'   \code{value} (1 = child labour including hazardous conditions), and
#'   \code{indicator} (\code{"childlabourhaz"}).
#' @export
process_CHILDLABOURHAZ <- function(fs) {
  parts <- .nga_child_labour_components(fs, "process_CHILDLABOURHAZ")

  # includeHaz = maximum(0, eaMore, hhcMore, hazardConditions).
  parts$fs$value <- as.integer(
    parts$ea_more | parts$hhc_more | parts$hazard
  )

  .indicator_result(parts$fs, "childlabourhaz")
}


#' Female genital mutilation among daughters age 0-14 (Nigeria)
#'
#' Percentage of living daughters age 0-14 years of women age 15-49 who have
#' had any form of female genital mutilation.
#'
#' Translated from \code{MICS6 - 09 - PR.5.3.sps} (\code{anyFgm}; MICS
#' indicator PR.10, SDG indicator 5.3.2).
#'
#' The denominator is assembled from two sources, as the SPSS does. Women who
#' have heard of FGM are asked the daughter roster, which forms \code{fg.sav};
#' women who have not (\code{FG2 >= 2}) skip it, so their living daughters are
#' taken from the birth history (\code{bh.sav}) and treated as not cut
#' (\code{FG17 = 2} in the syntax). Using \code{fg.sav} alone would restrict
#' the denominator to daughters of FGM-aware mothers and overstate prevalence.
#'
#' This indicator is Nigeria-only and needs three recode files; it expects the
#' Nigeria MICS6 (2021) variable names.
#'
#' @param fg A daughters-roster data.frame (fg.sav) from Nigeria MICS6, one row
#'   per daughter.
#' @param bh A birth history data.frame (bh.sav) from Nigeria MICS6.
#' @param wm A women's data.frame (wm.sav) from Nigeria MICS6, used to identify
#'   mothers who were not asked the daughter roster.
#'
#' @return A data.frame with columns
#'   \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#'   \code{value} (1 = daughter had any form of FGM), and \code{indicator}
#'   (\code{"fgmdaughter"}), one row per daughter age 0-14.
#' @export
process_FGMDAUGHTER <- function(fg, bh, wm) {
  # --- Source 1: the daughter roster (women who had heard of FGM) ----------
  # fg.sav is a woman-level file carrying wmweight, so it standardizes as "wm".
  fg <- standardize_columns(fg, "wm")
  .require_columns(
    fg,
    c("cluster", "householdID", "weight", "strata", "FG15", "FG17"),
    "process_FGMDAUGHTER"
  )

  # select if (FG15 >= 0 and FG15 <= 14).
  fg_age <- .as_code(fg$FG15)
  fg <- fg[fg_age >= 0 & fg_age <= 14 & !is.na(fg_age), , drop = FALSE]

  # anyFgm = 100 if (FG17 = 1).
  fg$value <- as.integer(.as_code(fg$FG17) %in% 1)
  roster <- .indicator_result(fg, "fgmdaughter")

  # --- Source 2: living daughters of women not asked the roster ------------
  wm <- standardize_columns(wm, "wm")
  .require_columns(
    wm, c("cluster", "householdID", "LN", "FG2"), "process_FGMDAUGHTER"
  )

  bh <- standardize_columns(bh, "bh")
  .require_columns(
    bh,
    c("cluster", "householdID", "weight", "strata", "LN", "BH3", "BH5", "BH6"),
    "process_FGMDAUGHTER"
  )

  # fgtemp.sav: `select if (FG2 >= 2)`, i.e. women who had not heard of FGM
  # and so were never asked about their daughters.
  fg2 <- .as_code(wm$FG2)
  unaware <- wm[!is.na(fg2) & fg2 >= 2, c("cluster", "householdID", "LN")]
  unaware <- unique(as.data.frame(unaware))

  # `match files /table = fgtemp.sav /by HH1 HH2 LN` then
  # `select if FG2 >= 2 and BH3 = 2 and BH5 = 1`: living daughters of those
  # women. BH6 (age of child, years) supplies FG15.
  bh_age <- .as_code(bh$BH6)
  keep <- .as_code(bh$BH3) %in% 2 &      # female
    .as_code(bh$BH5) %in% 1 &            # still alive
    bh_age >= 0 & bh_age <= 14 & !is.na(bh_age)
  bh <- bh[keep, , drop = FALSE]

  bh <- merge(
    as.data.frame(bh),
    unaware,
    by = c("cluster", "householdID", "LN")
  )

  # The syntax sets FG17 = 2 for these daughters, so anyFgm is 0 throughout;
  # they enter the denominator only.
  bh$value <- 0L
  from_births <- .indicator_result(bh, "fgmdaughter")

  rbind(as.data.frame(roster), as.data.frame(from_births))
}
