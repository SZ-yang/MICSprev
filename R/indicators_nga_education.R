# R/indicators_nga_education.R
#
# Nigeria-only education indicators, translated from the official MICS6 SPSS
# tabulation syntax (mics.unicef.org/tools, "42 Syntax Files 20220509") against
# the Nigeria MICS6 (2021) recode variable names.
#
# Source files:
#   MICS6 - 08 - LN.2.3.sps  Attendance, primary age        -> OOSPRIMARY
#   MICS6 - 08 - LN.2.4.sps  Attendance, lower secondary    -> OOSSECONDARY
#   MICS6 - 08 - LN.4.1.sps  Foundational reading skills    -> READING
#   MICS6 - 08 - LN.4.2.sps  Foundational numeracy skills   -> MATH
# School-structure parameters live in define/MICS6 - 08 - LN.sps.
#
# FLAG: define/MICS6 - 08 - LN.sps is a template whose header states the school
# structure "is country-specific and should be changed to reflect the situation
# in your country". The defaults below encode Nigeria's 6-3-3 system: primary
# is 6 years from age 6 (ages 6-11), junior secondary 3 years (ages 12-14),
# senior secondary 3 years (ages 15-17). Note the template ships
# UpSecSchoolGrades = 4, which is wrong for Nigeria; it does not affect the two
# indicators here but will matter for upper secondary completion.


# Nigeria school structure, matching the parameter names in
# define/MICS6 - 08 - LN.sps.
.nga_school_ages <- function(primary_entry_age = 6,
                             primary_grades = 6,
                             lower_sec_grades = 3) {
  low_sec_entry <- primary_entry_age + primary_grades
  list(
    primary_entry = primary_entry_age,
    primary_completion = low_sec_entry - 1,
    low_sec_entry = low_sec_entry,
    low_sec_completion = low_sec_entry + lower_sec_grades - 1,
    primary_grades = primary_grades,
    lower_sec_grades = lower_sec_grades
  )
}


# Shared body for LN.2.3 and LN.2.4. Both tables compute out-of-school
# identically; only the age window differs.
#
#   compute outOfSchool = 0.
#   if (ED4 = 2 or ED9 = 2) outOfSchool = 100.
#
# Per the LN.2.3 header notes, children whose attendance is unknown (ED9 > 2)
# are deliberately NOT counted as out of school, which `%in% 2` reproduces.
.nga_out_of_school <- function(hl, age_lo, age_hi, fn, label) {
  hl <- standardize_columns(hl, "hl")

  required <- c(
    "cluster", "householdID", "weight", "strata", "schage", "ED4", "ED9"
  )
  .require_columns(hl, required, fn)

  # select if (schage >= <entry age> and schage <= <completion age>).
  schage <- .as_code(hl$schage)
  keep <- !is.na(schage) & schage >= age_lo & schage <= age_hi
  hl <- hl[keep, , drop = FALSE]

  hl$value <- as.integer(
    .as_code(hl$ED4) %in% 2 | .as_code(hl$ED9) %in% 2
  )

  .indicator_result(hl, label)
}


#' Out-of-school rate, children of primary school age (Nigeria)
#'
#' Percentage of children of primary school age (age at the beginning of the
#' school year 6-11 in Nigeria) who did not attend any level of education or
#' early childhood education in the current school year.
#'
#' Translated from \code{MICS6 - 08 - LN.2.3.sps} (\code{outOfSchool}). Note
#' that children whose attendance is unknown (\code{ED9 > 2}) are not counted
#' as out of school, per the notes in that file.
#'
#' The age window is country-specific. Defaults encode Nigeria's structure
#' (primary entry at age 6, 6 grades); see \code{define/MICS6 - 08 - LN.sps}.
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021)
#' household listing variable names.
#'
#' @param hl A household listing data.frame (hl.sav) from Nigeria MICS6.
#' @param primary_entry_age Age of entry into primary school. Default 6.
#' @param primary_grades Number of grades in primary school. Default 6.
#'
#' @return A data.frame with columns
#'   \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#'   \code{value} (1 = out of school), and \code{indicator}
#'   (\code{"oosprimary"}).
#' @export
process_OOSPRIMARY <- function(hl, primary_entry_age = 6, primary_grades = 6) {
  ages <- .nga_school_ages(primary_entry_age, primary_grades)

  .nga_out_of_school(
    hl,
    age_lo = ages$primary_entry,
    age_hi = ages$primary_completion,
    fn = "process_OOSPRIMARY",
    label = "oosprimary"
  )
}


#' Out-of-school rate, children of lower secondary school age (Nigeria)
#'
#' Percentage of children of lower secondary school age (age at the beginning
#' of the school year 12-14 in Nigeria) who did not attend any level of
#' education or early childhood education in the current school year.
#'
#' Translated from \code{MICS6 - 08 - LN.2.4.sps} (\code{outOfSchool}), which
#' uses the same numerator definition as LN.2.3 over a different age window.
#'
#' The age window is country-specific. Defaults encode Nigeria's structure
#' (junior secondary, 3 years from age 12); see
#' \code{define/MICS6 - 08 - LN.sps}.
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021)
#' household listing variable names.
#'
#' @param hl A household listing data.frame (hl.sav) from Nigeria MICS6.
#' @param primary_entry_age Age of entry into primary school. Default 6.
#' @param primary_grades Number of grades in primary school. Default 6.
#' @param lower_sec_grades Number of grades in lower secondary school.
#'   Default 3.
#'
#' @return A data.frame with columns
#'   \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#'   \code{value} (1 = out of school), and \code{indicator}
#'   (\code{"oossecondary"}).
#' @export
process_OOSSECONDARY <- function(hl,
                                 primary_entry_age = 6,
                                 primary_grades = 6,
                                 lower_sec_grades = 3) {
  ages <- .nga_school_ages(primary_entry_age, primary_grades, lower_sec_grades)

  .nga_out_of_school(
    hl,
    age_lo = ages$low_sec_entry,
    age_hi = ages$low_sec_completion,
    fn = "process_OOSSECONDARY",
    label = "oossecondary"
  )
}


# Denominator shared by LN.4.1 and LN.4.2: children age 7-14 who completed the
# foundational learning module.
#
#   select if (FL28 = 1).
#   select if CB3 >= 7 and CB3 <= 14.
.nga_foundational_denominator <- function(fs, fn, extra_required) {
  fs <- standardize_columns(fs, "fs")

  required <- c(
    "cluster", "householdID", "weight", "strata", "FL28", "CB3", extra_required
  )
  .require_columns(fs, required, fn)

  age <- .as_code(fs$CB3)
  keep <- .as_code(fs$FL28) %in% 1 & !is.na(age) & age >= 7 & age <= 14
  fs[keep, , drop = FALSE]
}


#' Foundational reading skills, children age 7-14 (Nigeria)
#'
#' Percentage of children age 7-14 years who demonstrate foundational reading
#' skills: reading at least 90 percent of the words in a story correctly, and
#' correctly answering three literal and two inferential comprehension
#' questions.
#'
#' Translated from \code{MICS6 - 08 - LN.4.1.sps} (\code{readCorrect},
#' \code{aLiteral}, \code{aInferential}, \code{readingSkill}; MICS indicator
#' LN.22, SDG indicator 4.1.1).
#'
#' Two customisations are required for Nigeria and are exposed as arguments:
#'
#' \itemize{
#'   \item \strong{Story lengths.} The template hard-codes 72 and 61 words.
#'     Nigeria administered the reading tasks in several languages whose
#'     passages differ in length, so \code{fs.sav} carries word slots
#'     \code{FL19W1}-\code{FL19W88} and \code{FL21OW1}-\code{FL21OW63}. The
#'     slot counts are \emph{not} the story lengths: 6,861 children have data
#'     through slot 72 but only 675 reach slot 73 and 68 reach slot 79, and
#'     \code{FL20A} (words attempted) has a pronounced mode at exactly 72.
#'     Story 2 behaves the same way with a mode at 63. The defaults therefore
#'     take the modal passage lengths, 72 and 63. Since the SPSS applies one
#'     threshold to every child regardless of language, children given a longer
#'     passage are held to a slightly easier standard; this is inherent to the
#'     tabulation, not to the translation.
#'   \item \strong{Story 2 comprehension questions.} The template splits six
#'     questions as literal \code{FL22A}/\code{B}/\code{C} and inferential
#'     \code{FL22E}/\code{FL22F}. Nigeria asked only five (\code{FL22A}-
#'     \code{FL22E}), so \code{FL22F} does not exist. The defaults here mirror
#'     the 3-literal/2-inferential split of story 1 and the question wording:
#'     literal \code{FL22A} ("What class is Moses in"), \code{FL22B} ("What did
#'     Moses see on the way home"), \code{FL22D} ("Where did Moses fall");
#'     inferential \code{FL22C} ("Why did Moses start crying"), \code{FL22E}
#'     ("Why was Moses happy"). Override if the Nigeria report used a different
#'     split.
#' }
#'
#' Note also that \code{readCorrect} is a single flag across both stories, so
#' the SPSS permits the literal and inferential questions to be satisfied in
#' different stories even though the file's header comment describes the
#' composite as requiring "all of the above in the same language". The literal
#' translation of the code is kept here.
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021)
#' children age 5-17 recode variable names.
#'
#' @param fs A children age 5-17 data.frame (fs.sav) from Nigeria MICS6.
#' @param story1_words Number of words in the first reading passage.
#'   Default 72, the modal passage length in the Nigeria data.
#' @param story2_words Number of words in the second reading passage.
#'   Default 63, the modal passage length in the Nigeria data.
#' @param story2_literal Character vector of the story 2 literal comprehension
#'   variables. Default \code{c("FL22A", "FL22B", "FL22D")}.
#' @param story2_inferential Character vector of the story 2 inferential
#'   comprehension variables. Default \code{c("FL22C", "FL22E")}.
#'
#' @return A data.frame with columns
#'   \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#'   \code{value} (1 = demonstrates foundational reading skills), and
#'   \code{indicator} (\code{"reading"}).
#' @export
process_READING <- function(fs,
                            story1_words = 72,
                            story2_words = 63,
                            story2_literal = c("FL22A", "FL22B", "FL22D"),
                            story2_inferential = c("FL22C", "FL22E")) {
  story1_literal <- c("FL21BA", "FL21BB", "FL21BD")
  story1_inferential <- c("FL21BC", "FL21BE")

  extra <- c(
    "FL20A", "FL20B", "FL21PA", "FL21PB",
    story1_literal, story1_inferential, story2_literal, story2_inferential
  )
  fs <- .nga_foundational_denominator(fs, "process_READING", extra)

  # compute target1 = 0.
  # if (FL20A < 99 and FL20B < 99) target1 = FL20A - FL20B.
  .target <- function(attempted, incorrect) {
    a <- .as_code(fs[[attempted]])
    b <- .as_code(fs[[incorrect]])
    out <- rep(0, nrow(fs))
    ok <- !is.na(a) & !is.na(b) & a < 99 & b < 99
    out[ok] <- a[ok] - b[ok]
    out
  }

  target1 <- .target("FL20A", "FL20B")
  target2 <- .target("FL21PA", "FL21PB")

  # readCorrect = 100 if (target1 >= 0.9 * <story1>) or (target2 >= 0.9 * <story2>).
  read_correct <- target1 >= 0.9 * story1_words | target2 >= 0.9 * story2_words

  # All listed items answered correctly (code 1).
  .all_correct <- function(vars) {
    hits <- lapply(vars, function(v) .as_code(fs[[v]]) %in% 1)
    Reduce(`&`, hits)
  }

  # aLiteral / aInferential fire on either story, gated on readCorrect.
  a_literal <- read_correct &
    (.all_correct(story1_literal) | .all_correct(story2_literal))
  a_inferential <- read_correct &
    (.all_correct(story1_inferential) | .all_correct(story2_inferential))

  # readingSkill = readCorrect and aLiteral and aInferential.
  fs$value <- as.integer(read_correct & a_literal & a_inferential)

  .indicator_result(fs, "reading")
}


#' Foundational numeracy skills, children age 7-14 (Nigeria)
#'
#' Percentage of children age 7-14 years who demonstrate foundational numeracy
#' skills: successfully completing the number reading, number discrimination,
#' addition, and pattern recognition and completion tasks.
#'
#' Translated from \code{MICS6 - 08 - LN.4.2.sps} (\code{numberRead},
#' \code{numberDiscr}, \code{numberAdd}, \code{numberPattern}, \code{numSkill};
#' MICS indicator LN.23, SDG indicator 4.1.1). Each task requires every item in
#' its block to be answered correctly (code 1).
#'
#' This indicator is Nigeria-only; it expects the Nigeria MICS6 (2021)
#' children age 5-17 recode variable names.
#'
#' @param fs A children age 5-17 data.frame (fs.sav) from Nigeria MICS6.
#'
#' @return A data.frame with columns
#'   \code{cluster}, \code{householdID}, \code{weight}, \code{strata},
#'   \code{value} (1 = demonstrates foundational numeracy skills), and
#'   \code{indicator} (\code{"math"}).
#' @export
process_MATH <- function(fs) {
  number_read <- c("FL23A", "FL23B", "FL23C", "FL23D", "FL23E", "FL23F")
  number_discr <- c("FL24A", "FL24B", "FL24C", "FL24D", "FL24E")
  number_add <- c("FL25A", "FL25B", "FL25C", "FL25D", "FL25E")
  # FLAG: the pattern task is FL27, not FL26. Nigeria's FL26A-FL26E are
  # interviewer instructions ("Manual intro"), consistent with LN.4.2 keying
  # numberPattern off FL27A-FL27E.
  number_pattern <- c("FL27A", "FL27B", "FL27C", "FL27D", "FL27E")

  extra <- c(number_read, number_discr, number_add, number_pattern)
  fs <- .nga_foundational_denominator(fs, "process_MATH", extra)

  .all_correct <- function(vars) {
    hits <- lapply(vars, function(v) .as_code(fs[[v]]) %in% 1)
    Reduce(`&`, hits)
  }

  # numberRead uses `count numberReadTarget = FL23A ... FL23F (1)` followed by
  # `if (numberReadTarget = 6)`, i.e. all six correct.
  fs$value <- as.integer(
    .all_correct(number_read) &
      .all_correct(number_discr) &
      .all_correct(number_add) &
      .all_correct(number_pattern)
  )

  .indicator_result(fs, "math")
}
