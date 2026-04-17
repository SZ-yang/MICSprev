# R/process_geo_mics.R
# Exported wrapper

#' Process geo for supported MICS country-year combinations
#'
#' Users provide GPS shapefile; this constructs cluster/admin objects
#' needed by surveyPrev and SUMMER.
#'
#' Supported (currently, excluding Malawi):
#' - Honduras 2019
#' - Nigeria 2021
#' - Ghana 2017-18 (use year = 2017 or 2018)
#' - Lao PDR 2017
#' - Gambia 2018
#' - Sierra Leone 2017
#' - Madagascar 2018
#' - Zimbabwe 2019
#' - Thailand 2022
#'
#' @param country.name Character. Must match a supported handler key (aliases later).
#' @param gps_path Path to GPS .shp
#' @param out_path Optional .RData output path
#' @param year Numeric/integer or character. Used to select country-year handler.
#'
#' @return A list with cluster.info, admin.info1, admin.info2, admin0/admin1/admin2, geo
#' @export
process_geo_mics <- function(country.name, gps_path, out_path = NULL, year) {
  old_s2 <- .geo_set_s2_off()
  on.exit(.geo_restore_s2(old_s2), add = TRUE)

  key_country <- tolower(trimws(country.name))
  key_year <- as.character(year)
  key <- paste0(key_country, "__", key_year)

  handlers <- list(
    "honduras__2019"     = process_geo_hnd_2019,
    "nigeria__2021"      = process_geo_nga_2021,
    "ghana__2017"        = process_geo_gha_2017,  # Ghana 2017-18
    "ghana__2018"        = process_geo_gha_2017,  # Ghana 2017-18
    "lao pdr__2017"      = process_geo_lao_2017,  # alias support later
    "gambia__2018"       = process_geo_gmb_2018,
    "sierra leone__2017" = process_geo_sle_2017,
    "madagascar__2018"   = process_geo_mdg_2018,
    "zimbabwe__2019"     = process_geo_zwe_2019,
    "thailand__2022"     = process_geo_tha_2022
  )

  if (!key %in% names(handlers)) {
    stop(
      "Unsupported country/year: ", country.name, " / ", year, "\n\n",
      "Supported keys:\n- ", paste(names(handlers), collapse = "\n- ")
    )
  }

  handlers[[key]](gps_path = gps_path, out_path = out_path)
}
