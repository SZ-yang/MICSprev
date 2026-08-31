# R/indicators_nga_socialpolicy.R
#
# Nigeria-only social policy indicators, translated from the official MICS6
# SPSS tabulation syntax (mics.unicef.org/tools, "42 Syntax Files 20220509")
# against the Nigeria MICS6 (2021) recode variable names.
#
# Source files:
#   MICS6 - 10 - WS.2.1.sps  Handwashing facility          -> HANDWASH
#   MICS6 - 11 - EQ.1.4.sps  Child functioning             -> FUNCDIFF517,
#                                                             FUNCDIFF217
#   MICS6 - 11 - EQ.2.3.sps  Health insurance, under 5     -> HEALTHINS
#   MICS6 - 11 - EQ.2.7.sps  Social transfers, children    -> SOCTRANSFER
#   (no standard table)      Financial inclusion, women    -> BANKACCT,
#                                                             BORROWED


#' Handwashing facility with water and soap (Nigeria)
#'
#' Percentage of the household population with a handwashing facility where
#' water and soap or detergent are present.
#'
#' Translated from \code{MICS6 - 10 - WS.2.1.sps} (\code{indWS7}; MICS
#' indicator WS.7, SDG indicator 6.2.1b).
#'
#' The denominator is not all households. The SPSS computes \code{indWS7}
#' inside \code{do if (fixedobserved = 100 or mobileobserved = 100 or noplace
#' = 100)}, i.e. households where a facility was observed (\code{HW1} 1-3) or
#' where there is no facility in the dwelling, yard or plot (\code{HW1 = 4}).
#' Households that did not give permission to see the facility or where it was
#' not observed for another reason (\code{HW1 >= 5}) are excluded entirely.
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
#'   \code{strata}, \code{value} (1 = water and soap present), and
#'   \code{indicator} (\code{"handwash"}).
#' @export
process_HANDWASH <- function(hh) {
  hh <- .nga_hh_population_denominator(
    hh, "process_HANDWASH", c("HW1", "HW2", "HW7A", "HW7B")
  )

  # fixedobserved (HW1 = 1, 2), mobileobserved (HW1 = 3), noplace (HW1 = 4).
  hw1 <- .as_code(hh$HW1)
  hh <- hh[hw1 %in% c(1, 2, 3, 4), , drop = FALSE]

  # indWS7 = 100 if (HW2 = 1 and (HW7A = "A" or HW7B = "B")).
  hh$value <- as.integer(
    .as_code(hh$HW2) %in% 1 &
      (.mr_yes(hh$HW7A, "A") | .mr_yes(hh$HW7B, "B"))
  )

  .indicator_result(hh, "handwash")
}


# `anyfuncdifficulty` for children age 5-17 from fs.sav, per EQ.1.4.
# Domains: seeing, hearing, walking, self-care, communication, learning,
# remembering, concentrating, accepting change, behaviour, making friends,
# anxiety, depression. All but the last two fire on codes 3 or 4 ("a lot of
# difficulty" / "cannot do at all"); anxiety and depression fire on code 1
# ("daily").
.nga_funcdiff_5_17 <- function(fs) {
  .lot <- function(v) .as_code(fs[[v]]) %in% c(3, 4)

  .lot("FCF6") |                                    # seeing
    .lot("FCF8") |                                  # hearing
    .lot("FCF10") | .lot("FCF11") |                 # walking
    .lot("FCF14") | .lot("FCF15") |
    .lot("FCF16") |                                 # self-care
    .lot("FCF17") | .lot("FCF18") |                 # communication
    .lot("FCF19") |                                 # learning
    .lot("FCF20") |                                 # remembering
    .lot("FCF21") |                                 # concentrating
    .lot("FCF22") |                                 # accepting change
    .lot("FCF23") |                                 # behaviour
    .lot("FCF24") |                                 # making friends
    .as_code(fs$FCF25) %in% 1 |                     # anxiety
    .as_code(fs$FCF26) %in% 1                       # depression
}

.nga_funcdiff_5_17_vars <- c(
  "FCF6", "FCF8", "FCF10", "FCF11", "FCF14", "FCF15", "FCF16", "FCF17",
  "FCF18", "FCF19", "FCF20", "FCF21", "FCF22", "FCF23", "FCF24", "FCF25",
  "FCF26"
)


# `anyfuncdifficulty` for children age 2-4 from ch.sav, per EQ.1.4.
# Domains: seeing, hearing, walking, fine motor, communication, learning,
# playing, controlling behaviour.
#
# FLAG: the standard syntax adds a controlling-behaviour domain keyed off
# UCF19 = 5 (note: 5, not 3 or 4). Nigeria's ch.sav has no UCF19 - its UCF
# module ends at UCF18 - so that domain is absent from the questionnaire and
# cannot be computed. The Nigeria age 2-4 figure therefore covers seven
# domains rather than eight and is, if anything, a slight undercount relative
# to the standard definition.
.nga_funcdiff_2_4 <- function(ch) {
  .lot <- function(v) .as_code(ch[[v]]) %in% c(3, 4)

  out <- .lot("UCF7") |                             # seeing
    .lot("UCF9") |                                  # hearing
    .lot("UCF11") | .lot("UCF12") | .lot("UCF13") | # walking
    .lot("UCF14") |                                 # fine motor
    .lot("UCF15") | .lot("UCF16") |                 # communication
    .lot("UCF17") |                                 # learning
    .lot("UCF18")                                   # playing

  # Only applied where the questionnaire carried the item.
  if ("UCF19" %in% names(ch)) {
    out <- out | .as_code(ch$UCF19) %in% 5          # controlling behaviour
  }

  out
}

.nga_funcdiff_2_4_vars <- c(
  "UCF7", "UCF9", "UCF11", "UCF12", "UCF13", "UCF14", "UCF15", "UCF16",
  "UCF17", "UCF18"
)


#' Functional difficulty in at least one domain, children age 5-17 (Nigeria)
#'
#' Percentage of children age 5-17 years with a functional difficulty in at
#' least one domain.
#'
#' Translated from the \code{fs.sav} branch of \code{MICS6 - 11 - EQ.1.4.sps}
#' (\code{anyfuncdifficulty}; MICS indicator EQ.1).
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021)
#' children age 5-17 recode variable names.
#'
#' @param fs A children age 5-17 data.frame (fs.sav) from Nigeria MICS6.
#'
#' @return A data.frame with columns
#'   \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#'   \code{value} (1 = functional difficulty in at least one domain), and
#'   \code{indicator} (\code{"funcdiff517"}).
#' @export
process_FUNCDIFF517 <- function(fs) {
  fs <- standardize_columns(fs, "fs")

  required <- c(
    "cluster", "householdID", "weight", "strata",
    "FS17", "CB3", .nga_funcdiff_5_17_vars
  )
  .require_columns(fs, required, "process_FUNCDIFF517")

  age <- .as_code(fs$CB3)
  keep <- .as_code(fs$FS17) %in% 1 & !is.na(age) & age >= 5 & age <= 17
  fs <- fs[keep, , drop = FALSE]

  fs$value <- as.integer(.nga_funcdiff_5_17(fs))

  .indicator_result(fs, "funcdiff517")
}


#' Functional difficulty in at least one domain, children age 2-17 (Nigeria)
#'
#' Percentage of children age 2-17 years with a functional difficulty in at
#' least one domain.
#'
#' Translated from \code{MICS6 - 11 - EQ.1.4.sps}, which builds this indicator
#' by stacking two recode files: children age 5-17 from \code{fs.sav} (the
#' \code{FCF} module, weighted by \code{fsweight}) and children age 2-4 from
#' \code{ch.sav} (the \code{UCF} module, weighted by \code{chweight}). The
#' syntax then applies \code{select if (result = 1)} and
#' \code{select if (age >= 2 and age <= 17)} to the stacked file and weights by
#' the carried-over weight.
#'
#' The two age ranges use different domain sets, as the questionnaires do: age
#' 2-4 covers seeing, hearing, walking, fine motor, communication, learning
#' and playing, while age 5-17 adds self-care, remembering, concentrating,
#' accepting change, making friends, anxiety and depression. Prevalence is
#' therefore not directly comparable between the two age bands, which is
#' inherent to the indicator.
#'
#' FLAG: the standard age 2-4 definition includes a controlling-behaviour
#' domain keyed off \code{UCF19} (code 5, not 3 or 4). Nigeria's \code{ch.sav}
#' has no \code{UCF19} - its UCF module ends at \code{UCF18} - so that domain
#' is absent from the questionnaire and the age 2-4 component covers seven
#' domains rather than eight.
#'
#' This indicator is Nigeria-only and needs two recode files; it expects the
#' Nigeria MICS6 (2021) variable names.
#'
#' @param fs A children age 5-17 data.frame (fs.sav) from Nigeria MICS6.
#' @param ch A children-under-5 data.frame (ch.sav) from Nigeria MICS6,
#'   supplying children age 2-4.
#'
#' @return A data.frame with columns
#'   \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#'   \code{value} (1 = functional difficulty in at least one domain), and
#'   \code{indicator} (\code{"funcdiff217"}), one row per child age 2-17.
#' @export
process_FUNCDIFF217 <- function(fs, ch) {
  older <- process_FUNCDIFF517(fs)
  older$indicator <- "funcdiff217"

  ch <- standardize_columns(ch, "ch")

  required <- c(
    "cluster", "householdID", "weight", "strata",
    "UF17", "UB2", .nga_funcdiff_2_4_vars
  )
  .require_columns(ch, required, "process_FUNCDIFF217")

  # ch.sav supplies the age 2-4 rows; UB2 is the age variable EQ.1.4 renames
  # to `age` before stacking.
  age <- .as_code(ch$UB2)
  keep <- .as_code(ch$UF17) %in% 1 & !is.na(age) & age >= 2 & age <= 4
  ch <- ch[keep, , drop = FALSE]

  ch$value <- as.integer(.nga_funcdiff_2_4(ch))
  younger <- .indicator_result(ch, "funcdiff217")

  rbind(as.data.frame(younger), as.data.frame(older))
}


#' Health insurance coverage, children under age 5 (Nigeria)
#'
#' Percentage of children under age 5 covered by any health insurance.
#'
#' Translated from \code{MICS6 - 11 - EQ.2.3.sps} (\code{hinsurance}; MICS
#' indicator EQ.3).
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021)
#' children's recode variable names.
#'
#' @param ch A children-under-5 data.frame (ch.sav) from Nigeria MICS6.
#'
#' @return A data.frame with columns
#'   \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#'   \code{value} (1 = covered by health insurance), and \code{indicator}
#'   (\code{"healthins"}).
#' @export
process_HEALTHINS <- function(ch) {
  ch <- standardize_columns(ch, "ch")

  required <- c("cluster", "householdID", "weight", "strata", "UF17", "UB9")
  .require_columns(ch, required, "process_HEALTHINS")

  ch <- ch[.as_code(ch$UF17) %in% 1, , drop = FALSE]

  # hinsurance = 100 if UB9 = 1.
  ch$value <- as.integer(.as_code(ch$UB9) %in% 1)

  .indicator_result(ch, "healthins")
}


#' Social transfers or benefits, children under age 18 (Nigeria)
#'
#' Percentage of children under age 18 living in households that received
#' social transfers or benefits in the last 3 months.
#'
#' Translated from \code{MICS6 - 11 - EQ.2.7.sps} (\code{any}; MICS indicator
#' EQ.5, SDG indicator 1.3.1). The indicator combines two sources:
#'
#' \itemize{
#'   \item Household programme receipt from \code{hh.sav}: a programme counts
#'     when \code{ST4U$k = 1} (unit is months) and \code{ST4N$k <= 3}.
#'   \item School support from \code{hl.sav}: any household member age 5-24
#'     attending primary or higher (\code{ED10A} 1-7) who received tuition
#'     (\code{ED12 = 1}) or material support (\code{ED14 = 1}), aggregated to
#'     the household.
#' }
#'
#' The denominator is children, not households: the syntax merges the household
#' result onto the household listing and applies \code{select if (HL6 < 18)},
#' weighting by \code{hhweight}. One row is returned per child under 18.
#'
#' FLAG: the source syntax names five programme slots
#' (\code{ST4U$1}-\code{ST4U$5}). Nigeria fielded four, so this function
#' discovers the \code{ST4U$*}/\code{ST4N$*} pairs actually present rather than
#' assuming a fixed count.
#'
#' This indicator is Nigeria-only and needs two recode files; it expects the
#' Nigeria MICS6 (2021) variable names.
#'
#' @param hh A household data.frame (hh.sav) from Nigeria MICS6.
#' @param hl A household listing data.frame (hl.sav) from Nigeria MICS6.
#'
#' @return A data.frame with columns
#'   \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#'   \code{value} (1 = household received a social transfer), and
#'   \code{indicator} (\code{"soctransfer"}), one row per child under 18.
#' @export
process_SOCTRANSFER <- function(hh, hl) {
  hh <- standardize_columns(hh, "hh")

  # Discover the programme slots present rather than hard-coding 1-5.
  unit_vars <- grep("^ST4U\\$[0-9]+$", names(hh), value = TRUE)
  num_vars <- sub("^ST4U", "ST4N", unit_vars)

  if (length(unit_vars) == 0) {
    stop(
      "process_SOCTRANSFER(): no ST4U$* programme variables found. This ",
      "indicator is Nigeria-only and expects the Nigeria MICS6 household ",
      "recode variable names.",
      call. = FALSE
    )
  }

  .require_columns(
    hh, c("cluster", "householdID", "HH46", unit_vars, num_vars),
    "process_SOCTRANSFER"
  )

  hh <- hh[.as_code(hh$HH46) %in% 1, , drop = FALSE]

  # Socialtransfer_k = 100 if (ST4U$k = 1 and ST4N$k <= 3).
  transfer <- Reduce(
    `|`,
    Map(
      function(u, n) {
        .as_code(hh[[u]]) %in% 1 & .as_code(hh[[n]]) <= 3
      },
      unit_vars, num_vars
    )
  )
  transfer[is.na(transfer)] <- FALSE

  hh_flags <- data.frame(
    cluster = hh$cluster,
    householdID = hh$householdID,
    hh_transfer = transfer
  )

  # School support, aggregated from the household listing.
  hl <- standardize_columns(hl, "hl")
  .require_columns(
    hl,
    c("cluster", "householdID", "weight", "strata",
      "HL6", "ED10A", "ED12", "ED14"),
    "process_SOCTRANSFER"
  )

  member_age <- .as_code(hl$HL6)
  level <- .as_code(hl$ED10A)

  # do if ((HL6 >= 5 and HL6 <= 24) and (ED10A >= 1 and ED10A < 8)).
  eligible <- !is.na(member_age) & member_age >= 5 & member_age <= 24 &
    !is.na(level) & level >= 1 & level < 8

  ed_support <- eligible &
    (.as_code(hl$ED12) %in% 1 | .as_code(hl$ED14) %in% 1)

  # aggregate /break = HH1 HH2 /nEDsupport = sum(EDsupport); EDsupport = 100
  # if nEDsupport > 0.
  hh_key <- paste(hl$cluster, hl$householdID, sep = "\r")
  ed_by_hh <- tapply(ed_support, hh_key, any)

  # Denominator: children under 18 in the household listing.
  hl <- hl[!is.na(member_age) & member_age < 18, , drop = FALSE]
  hl_key <- paste(hl$cluster, hl$householdID, sep = "\r")

  ed_flag <- unname(ed_by_hh[hl_key])
  ed_flag[is.na(ed_flag)] <- FALSE

  idx <- match(
    paste(hl$cluster, hl$householdID, sep = "\r"),
    paste(hh_flags$cluster, hh_flags$householdID, sep = "\r")
  )
  transfer_flag <- hh_flags$hh_transfer[idx]
  transfer_flag[is.na(transfer_flag)] <- FALSE

  # any = 100 if any programme transfer or EDsupport.
  hl$value <- as.integer(transfer_flag | ed_flag)

  .indicator_result(hl, "soctransfer")
}


#' Women age 15-49 who own a bank account (Nigeria)
#'
#' Percentage of women age 15-49 years who own a bank account or a similar
#' setup at a financial institution.
#'
#' \strong{This indicator has no standard MICS6 tabulation syntax.} Financial
#' inclusion is not among the 42 standard MICS6 tabulation files; \code{FN2}
#' ("Do you own a bank account") is specific to the Nigeria questionnaire. The
#' definition follows the MICS convention used throughout the standard tables:
#' an explicit "Yes" counts, and "No", "Don't know" (8), "No response" (9) and
#' missing all count as 0 within the full denominator of interviewed women.
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021) women's
#' recode variable names.
#'
#' @param wm A women's data.frame (wm.sav) from Nigeria MICS6.
#'
#' @return A data.frame with columns
#'   \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#'   \code{value} (1 = owns a bank account), and \code{indicator}
#'   (\code{"bankacct"}).
#' @export
process_BANKACCT <- function(wm) {
  wm <- standardize_columns(wm, "wm")

  required <- c("cluster", "householdID", "weight", "strata", "WM17", "FN2")
  .require_columns(wm, required, "process_BANKACCT")

  wm <- wm[.as_code(wm$WM17) %in% 1, , drop = FALSE]

  wm$value <- as.integer(.as_code(wm$FN2) %in% 1)

  .indicator_result(wm, "bankacct")
}


#' Women age 15-49 who borrowed money (Nigeria)
#'
#' Percentage of women age 15-49 years who borrowed money.
#'
#' \strong{This indicator has no standard MICS6 tabulation syntax.} It is
#' derived from \code{FN5} ("Where did you borrow most money from"), a Nigeria
#' questionnaire item whose response set includes an explicit
#' "Did not borrow" category (code 11). A woman counts as having borrowed when
#' she names any source: commercial bank (21), microfinance bank (22),
#' non-interest bank (23), family member (24), friend (25), co-operative group
#' (26), online platform (27), community money lender (28), or other (96).
#' "Did not borrow" (11), "No response" (99) and missing count as 0.
#'
#' \strong{FLAG - the reference period is not verifiable from the recode.} The
#' indicator is specified as "borrowed money in the last 12 months", but
#' \code{FN5} carries no time qualifier in its variable or value labels, so the
#' 12-month window rests on the questionnaire wording rather than anything in
#' the data. Confirm against the Nigeria questionnaire before publishing.
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021) women's
#' recode variable names.
#'
#' @param wm A women's data.frame (wm.sav) from Nigeria MICS6.
#' @param borrow_codes Numeric vector of \code{FN5} codes counted as having
#'   borrowed. Default \code{c(21:28, 96)}.
#'
#' @return A data.frame with columns
#'   \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#'   \code{value} (1 = borrowed money), and \code{indicator}
#'   (\code{"borrowed"}).
#' @export
process_BORROWED <- function(wm, borrow_codes = c(21:28, 96)) {
  wm <- standardize_columns(wm, "wm")

  required <- c("cluster", "householdID", "weight", "strata", "WM17", "FN5")
  .require_columns(wm, required, "process_BORROWED")

  wm <- wm[.as_code(wm$WM17) %in% 1, , drop = FALSE]

  wm$value <- as.integer(.as_code(wm$FN5) %in% borrow_codes)

  .indicator_result(wm, "borrowed")
}
