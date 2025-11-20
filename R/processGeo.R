#' Build geo objects for any MICS country using GPS shapefile
#'
#' This function automatically downloads administrative boundaries via
#' \code{surveyPrev:::get_geoBoundaries()}, links MICS GPS cluster points,
#' and constructs \code{cluster.info}, \code{admin.info1}, \code{admin.info2}
#' needed by \pkg{surveyPrev} and \pkg{SUMMER}.
#'
#' @param country.name Character. Country name recognizable by
#'   \code{countrycode} (e.g. "Honduras", "Sierra Leone").
#' @param gps_path Path to the MICS GPS shapefile (".shp").
#' @param out_path Optional. If provided, saves all outputs into an RData file.
#' @param year Optional. Only used for messaging and filename patterns.
#'
#' @return A list containing:
#'   \item{cluster.info}{Output of surveyPrev::clusterInfo()}
#'   \item{admin.info1}{Output of surveyPrev::adminInfo()}
#'   \item{admin.info2}{Output of surveyPrev::adminInfo()}
#'   \item{admin0, admin1, admin2}{sf administrative boundaries}
#'   \item{geo}{sf cluster points}
#'
#' @importFrom rlang .data
#' @export
build_geo_mics <- function(
    country.name,
    gps_path,
    out_path = NULL,
    year = NA
) {
  # Preserve current s2 setting
  old_s2 <- sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(FALSE)

  # ---- 1) Convert country to ISO3 code ----
  iso3 <- countrycode::countrycode(
    country.name,
    origin = "country.name",
    destination = "iso3c"
  )

  if (is.na(iso3)) {
    stop("Could not convert country.name to ISO3. Check spelling.")
  }

  message("\n=== Processing geo data for ", country.name, " (", iso3, ") ===")

  # ---- 2) Download admin boundaries ----
  admin0 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM0")
  admin1 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM1")
  admin2 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM2")

  # Standardize naming
  admin0$NAME_0 <- admin0$shapeName
  admin1$NAME_1 <- admin1$shapeName
  admin2$NAME_2 <- admin2$shapeName

  # Harmonize CRS
  if (!sf::st_crs(admin1) == sf::st_crs(admin2)) {
    admin1 <- sf::st_transform(admin1, crs = sf::st_crs(admin2))
  }

  # Add NAME_1 to admin2 (assign parent)
  admin2 <- sf::st_join(
    admin2,
    admin1 %>% dplyr::select(NAME_1),
    join = sf::st_intersects,
    largest = TRUE
  )


  # ---- 3) Read GPS / cluster shapefile ----
  geo_data <- sf::st_read(gps_path, quiet = TRUE)
  geo_data <- sf::st_transform(geo_data, sf::st_crs(admin1))

  if (!"HH1" %in% names(geo_data)) {
    stop("GPS file must contain HH1 (cluster ID).")
  }

  # geo <- geo_data %>%
  #   dplyr::mutate(
  #     DHSCLUST = as.integer(.data$HH1),
  #     LONGNUM  = sf::st_coordinates(geometry)[, 1],
  #     LATNUM   = sf::st_coordinates(geometry)[, 2],
  #     ADM1NAME = dplyr::if_else("GEONAME" %in% names(geo_data),
  #                               .data$GEONAME,
  #                               NA_character_)
  #   )

  if ("GEONAME" %in% names(geo_data)) {
    geo <- geo_data %>%
      dplyr::mutate(
        DHSCLUST = as.integer(.data$HH1),
        LONGNUM  = sf::st_coordinates(geometry)[, 1],
        LATNUM   = sf::st_coordinates(geometry)[, 2],
        ADM1NAME = .data$GEONAME
      )
  } else {
    geo <- geo_data %>%
      dplyr::mutate(
        DHSCLUST = as.integer(.data$HH1),
        LONGNUM  = sf::st_coordinates(geometry)[, 1],
        LATNUM   = sf::st_coordinates(geometry)[, 2],
        ADM1NAME = NA_character_
      )
  }

  message("Loaded ", nrow(geo), " cluster points.")


  # ---- 4) Build cluster/admin info ----
  cluster.info <- surveyPrev::clusterInfo(
    geo       = geo,
    poly.adm1 = admin1,
    poly.adm2 = admin2,
    by.adm1   = "NAME_1",
    by.adm2   = "NAME_2",
    map       = TRUE
  )

  admin.info1 <- surveyPrev::adminInfo(
    poly.adm = admin1,
    admin    = 1,
    by.adm   = "NAME_1"
  )

  admin.info2 <- surveyPrev::adminInfo(
    poly.adm     = admin2,
    admin        = 2,
    by.adm       = "NAME_2",
    by.adm.upper = "NAME_1"
  )


  # ---- 5) Optional save ----
  if (!is.null(out_path)) {
    save(
      cluster.info, admin.info1, admin.info2,
      admin0, admin1, admin2, geo,
      file = out_path
    )
    message("Saved geo objects to: ", out_path)
  }

  message("\n=== Geo processing complete for ", country.name, " ===\n")

  invisible(list(
    cluster.info = cluster.info,
    admin.info1  = admin.info1,
    admin.info2  = admin.info2,
    admin0       = admin0,
    admin1       = admin1,
    admin2       = admin2,
    geo          = geo
  ))
}
