# R/indicator_helpers.R
# Shared internals for the country-specific indicator functions.

# Standard output contract for every process_<IND>() function.
.indicator_out_cols <- c(
  "cluster", "householdID", "weight", "strata", "value", "indicator"
)

# Validate that `df` carries every column in `required`, erroring with a message
# that names the missing columns and states which country the indicator covers.
#
# fn       - name of the calling function, for the error message
# required - character vector of required column names
# country  - country the indicator is gated to, or NULL for country-agnostic
.require_columns <- function(df, required, fn, country = "Nigeria") {
  missing <- setdiff(required, names(df))

  if (length(missing) == 0) return(invisible(TRUE))

  stop(
    fn, "(): missing required column(s): ",
    paste(missing, collapse = ", "),
    if (!is.null(country)) {
      paste0(
        ". This indicator is ", country, "-only and expects the ", country,
        " MICS6 recode variable names."
      )
    } else {
      "."
    },
    call. = FALSE
  )
}

# Coerce a haven_labelled / labelled column to plain numeric without tripping
# over value labels. MICS recodes are numeric-coded throughout, and every
# indicator here keys off codes rather than labels.
.as_code <- function(x) {
  if (is.null(x)) return(NULL)
  as.numeric(haven::zap_labels(x))
}

# Assemble the standard return frame. `value` and `indicator` must already be
# columns of `df`.
.indicator_result <- function(df, indicator) {
  df$indicator <- indicator
  df[, .indicator_out_cols, drop = FALSE]
}
