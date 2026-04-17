# R/geo_helpers.R
# Internal helper functions for geo processing

#' @importFrom rlang .data
#' @importFrom dplyr %>%
NULL

#' @keywords internal
.geo_set_s2_off <- function() {
  old_s2 <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)
  old_s2
}

#' @keywords internal
.geo_restore_s2 <- function(old_s2) {
  sf::sf_use_s2(old_s2)
}

#' @keywords internal
.read_gps_points <- function(gps_path, crs_target, require_geoname = FALSE) {
  geo_data <- sf::st_read(gps_path, quiet = TRUE)
  geo_data <- sf::st_transform(geo_data, crs_target)

  if (!"HH1" %in% names(geo_data)) {
    stop("GPS file must contain HH1 (cluster ID).")
  }

  if (!"GEONAME" %in% names(geo_data)) {
    if (require_geoname) stop("GPS file must contain GEONAME for this country.")
    geo_data$GEONAME <- NA_character_
  }

  geo_data %>%
    dplyr::mutate(
      DHSCLUST = as.integer(.data$HH1),
      LONGNUM  = sf::st_coordinates(.data$geometry)[, 1],
      LATNUM   = sf::st_coordinates(.data$geometry)[, 2],
      ADM1NAME = .data$GEONAME
    )
}

#' @keywords internal
.attach_parent_admin1 <- function(admin2, admin1, parent_col = "NAME_1") {
  if (sf::st_crs(admin1) != sf::st_crs(admin2)) {
    admin1 <- sf::st_transform(admin1, sf::st_crs(admin2))
  }

  sf::st_join(
    admin2,
    admin1 %>% dplyr::select(dplyr::all_of(parent_col)),
    join = sf::st_intersects,
    largest = TRUE
  )
}

#' @keywords internal
.build_cluster_admin_info <- function(geo, admin1, admin2, by_adm2, map = TRUE,
                                      proportion1 = NULL, proportion2 = NULL) {
  cluster.info <- surveyPrev::clusterInfo(
    geo       = geo,
    poly.adm1 = admin1,
    poly.adm2 = admin2,
    by.adm1   = "NAME_1",
    by.adm2   = by_adm2,
    map       = map
  )

  admin.info1 <- surveyPrev::adminInfo(
    poly.adm   = admin1,
    admin      = 1,
    by.adm     = "NAME_1",
    proportion = proportion1
  )

  admin.info2 <- surveyPrev::adminInfo(
    poly.adm      = admin2,
    admin         = 2,
    by.adm        = by_adm2,
    by.adm.upper  = "NAME_1",
    proportion    = proportion2
  )

  list(cluster.info = cluster.info, admin.info1 = admin.info1, admin.info2 = admin.info2)
}

#' @keywords internal
.save_geo_result <- function(res, out_path) {
  if (is.null(out_path)) return(invisible(NULL))
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  save(list = names(res), file = out_path, envir = list2env(res))
  invisible(NULL)
}
