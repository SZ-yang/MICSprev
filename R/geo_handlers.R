# R/geo_handlers.R
# Country-year specific geo processors (internal)

# Honduras 2019
#' @keywords internal
process_geo_hnd_2019 <- function(gps_path, out_path = NULL) {
  iso3 <- countrycode::countrycode("Honduras", "country.name", "iso3c")

  admin0 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM0")
  admin1 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM1")
  admin2 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM2")

  admin0$NAME_0 <- admin0$shapeName
  admin1$NAME_1 <- admin1$shapeName
  admin2$NAME_2 <- admin2$shapeName

  admin2 <- .attach_parent_admin1(admin2, admin1, parent_col = "NAME_1")
  geo <- .read_gps_points(gps_path, sf::st_crs(admin1))

  infos <- .build_cluster_admin_info(geo, admin1, admin2, by_adm2 = "NAME_2", map = TRUE)

  res <- c(infos, list(admin0 = admin0, admin1 = admin1, admin2 = admin2, geo = geo))
  .save_geo_result(res, out_path)
  res
}

# Nigeria 2021
#' @keywords internal
process_geo_nga_2021 <- function(gps_path, out_path = NULL) {
  iso3 <- countrycode::countrycode("Nigeria", "country.name", "iso3c")

  admin0 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM0")
  admin1 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM1")
  admin2 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM2")

  admin0$NAME_0 <- admin0$shapeName
  admin1$NAME_1 <- admin1$shapeName
  admin2$NAME_2 <- admin2$shapeName

  admin2 <- .attach_parent_admin1(admin2, admin1, parent_col = "NAME_1")
  admin2$admin2.name.full <- paste0(admin2$NAME_1, "_", admin2$NAME_2)

  geo <- .read_gps_points(gps_path, sf::st_crs(admin1))

  infos <- .build_cluster_admin_info(geo, admin1, admin2, by_adm2 = "NAME_2", map = TRUE)

  res <- c(infos, list(admin0 = admin0, admin1 = admin1, admin2 = admin2, geo = geo))
  .save_geo_result(res, out_path)
  res
}


# Nigeria 2016-17
#' @keywords internal
process_geo_nga_2016 <- function(gps_path, out_path = NULL) {

  iso3 <- countrycode::countrycode("Nigeria", "country.name", "iso3c")

  # Use the same boundary source as Nigeria 2021
  admin0 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM0")
  admin1 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM1")
  admin2 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM2")

  admin0$NAME_0 <- admin0$shapeName
  admin1$NAME_1 <- admin1$shapeName
  admin2$NAME_2 <- admin2$shapeName

  # Attach admin1 name to each admin2 polygon spatially
  admin2 <- .attach_parent_admin1(admin2, admin1, parent_col = "NAME_1")

  # Important for admin2 map joins later
  admin2$admin2.name.full <- paste0(admin2$NAME_1, "_", admin2$NAME_2)

  geo <- .read_gps_points(gps_path, sf::st_crs(admin1))

  infos <- .build_cluster_admin_info(
    geo,
    admin1,
    admin2,
    by_adm2 = "NAME_2",
    map = TRUE
  )

  res <- c(
    infos,
    list(
      admin0 = admin0,
      admin1 = admin1,
      admin2 = admin2,
      geo = geo
    )
  )

  .save_geo_result(res, out_path)

  res
}


# Ghana 2017-18 (use year 2017 or 2018 in wrapper)
#' @keywords internal
process_geo_gha_2017 <- function(gps_path, out_path = NULL) {
  iso3 <- countrycode::countrycode("Ghana", "country.name", "iso3c")

  admin0 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM0")
  admin1 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM1")
  admin2 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM2")

  admin0$NAME_0 <- admin0$shapeName
  admin1$NAME_1 <- admin1$shapeName
  admin2$NAME_2 <- admin2$shapeName

  admin2 <- .attach_parent_admin1(admin2, admin1, parent_col = "NAME_1")

  # dissolve duplicates like coworker's script
  dup1 <- table(admin1$NAME_1)
  if (any(dup1 > 1)) {
    admin1 <- admin1 %>%
      dplyr::group_by(.data$NAME_1) %>%
      dplyr::summarise(.groups = "drop") %>%
      sf::st_make_valid()
  }

  dup2 <- admin2 %>%
    sf::st_drop_geometry() %>%
    dplyr::count(.data$NAME_1, .data$NAME_2, name = "n") %>%
    dplyr::filter(.data$n > 1)

  if (nrow(dup2) > 0) {
    admin2 <- admin2 %>%
      dplyr::group_by(.data$NAME_1, .data$NAME_2) %>%
      dplyr::summarise(.groups = "drop") %>%
      sf::st_make_valid()
  }

  admin0 <- sf::st_make_valid(admin0)
  admin1 <- sf::st_make_valid(admin1)
  admin2 <- sf::st_make_valid(admin2)

  geo <- .read_gps_points(gps_path, sf::st_crs(admin1))

  infos <- .build_cluster_admin_info(geo, admin1, admin2, by_adm2 = "NAME_2", map = TRUE)

  res <- c(infos, list(admin0 = admin0, admin1 = admin1, admin2 = admin2, geo = geo))
  .save_geo_result(res, out_path)
  res
}

# Lao PDR 2017
#' @keywords internal
process_geo_lao_2017 <- function(gps_path, out_path = NULL) {
  iso3 <- countrycode::countrycode("Lao People's Democratic Republic", "country.name", "iso3c")

  admin0 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM0")
  admin1 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM1")
  admin2 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM2")

  admin0$NAME_0 <- admin0$shapeName
  admin1$NAME_1 <- admin1$shapeName
  admin2$NAME_2 <- admin2$shapeName

  admin2 <- .attach_parent_admin1(admin2, admin1, parent_col = "NAME_1")
  admin2$admin1.name.full <- paste0(admin2$NAME_1, "_", admin2$NAME_2)

  geo <- .read_gps_points(gps_path, sf::st_crs(admin1))

  infos <- .build_cluster_admin_info(geo, admin1, admin2, by_adm2 = "NAME_2", map = TRUE)

  res <- c(infos, list(admin0 = admin0, admin1 = admin1, admin2 = admin2, geo = geo))
  .save_geo_result(res, out_path)
  res
}

# Gambia 2018
#' @keywords internal
process_geo_gmb_2018 <- function(gps_path, out_path = NULL) {
  iso3 <- countrycode::countrycode("Gambia", "country.name", "iso3c")

  admin0 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM0")
  admin1 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM1")
  admin2 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM2")

  admin0$NAME_0 <- admin0$shapeName
  admin1$NAME_1 <- admin1$shapeName
  admin2$NAME_2 <- admin2$shapeName

  admin2 <- .attach_parent_admin1(admin2, admin1, parent_col = "NAME_1")

  geo <- .read_gps_points(gps_path, sf::st_crs(admin1))

  infos <- .build_cluster_admin_info(geo, admin1, admin2, by_adm2 = "NAME_2", map = TRUE)

  res <- c(infos, list(admin0 = admin0, admin1 = admin1, admin2 = admin2, geo = geo))
  .save_geo_result(res, out_path)
  res
}

# Sierra Leone 2017 (ADM2/ADM3)
#' @keywords internal
process_geo_sle_2017 <- function(gps_path, out_path = NULL) {
  iso3 <- countrycode::countrycode("Sierra Leone", "country.name", "iso3c")

  admin0 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM0")
  admin1 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM2")  # modeling admin1
  admin2 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM3")  # modeling admin2

  admin0$NAME_0 <- admin0$shapeName
  admin1$NAME_1 <- admin1$shapeName
  admin2$NAME_2 <- admin2$shapeName

  admin2 <- .attach_parent_admin1(admin2, admin1, parent_col = "NAME_1")

  geo_data <- sf::st_read(gps_path, quiet = TRUE)
  geo_data <- sf::st_transform(geo_data, sf::st_crs(admin1))

  if (!"HH1" %in% names(geo_data)) stop("GPS file must contain HH1 (cluster ID).")

  geo <- geo_data %>%
    dplyr::mutate(
      DHSCLUST = as.integer(.data$HH1),
      LONGNUM  = sf::st_coordinates(.data$geometry)[, 1],
      LATNUM   = sf::st_coordinates(.data$geometry)[, 2]
    )

  # derive ADM1NAME from polygons (coworker logic)
  geo <- sf::st_join(
    geo,
    admin1 %>% dplyr::select(.data$NAME_1),
    join = sf::st_intersects,
    largest = TRUE
  ) %>%
    dplyr::mutate(ADM1NAME = .data$NAME_1) %>%
    dplyr::select(-.data$NAME_1)

  infos <- .build_cluster_admin_info(geo, admin1, admin2, by_adm2 = "NAME_2", map = TRUE)

  res <- c(infos, list(admin0 = admin0, admin1 = admin1, admin2 = admin2, geo = geo))
  .save_geo_result(res, out_path)
  res
}

# Madagascar 2018
#' @keywords internal
process_geo_mdg_2018 <- function(gps_path, out_path = NULL) {
  iso3 <- countrycode::countrycode("Madagascar", "country.name", "iso3c")

  admin0 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM0")
  admin1 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM1")
  admin2 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM2")

  admin0$NAME_0 <- admin0$shapeName
  admin1$NAME_1 <- admin1$shapeName
  admin2$NAME_2 <- admin2$shapeName

  admin2 <- .attach_parent_admin1(admin2, admin1, parent_col = "NAME_1")

  geo <- .read_gps_points(gps_path, sf::st_crs(admin1))

  infos <- .build_cluster_admin_info(geo, admin1, admin2, by_adm2 = "NAME_2", map = TRUE)

  res <- c(infos, list(admin0 = admin0, admin1 = admin1, admin2 = admin2, geo = geo))
  .save_geo_result(res, out_path)
  res
}

# Zimbabwe 2019
#' @keywords internal
process_geo_zwe_2019 <- function(gps_path, out_path = NULL) {
  iso3 <- countrycode::countrycode("Zimbabwe", "country.name", "iso3c")

  admin0 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM0")
  admin1 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM1")
  admin2 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM2")

  admin0$NAME_0 <- admin0$shapeName
  admin1$NAME_1 <- admin1$shapeName
  admin2$NAME_2 <- admin2$shapeName

  admin2 <- .attach_parent_admin1(admin2, admin1, parent_col = "NAME_1")

  geo <- .read_gps_points(gps_path, sf::st_crs(admin1))

  infos <- .build_cluster_admin_info(geo, admin1, admin2, by_adm2 = "NAME_2", map = TRUE)

  res <- c(infos, list(admin0 = admin0, admin1 = admin1, admin2 = admin2, geo = geo))
  .save_geo_result(res, out_path)
  res
}

# Thailand 2022 (NAME_2_UID)
#' @keywords internal
process_geo_tha_2022 <- function(gps_path, out_path = NULL) {
  iso3 <- countrycode::countrycode("Thailand", "country.name", "iso3c")

  admin0 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM0")
  admin1 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM1")
  admin2 <- surveyPrev:::get_geoBoundaries(iso3, adm = "ADM2")

  admin0$NAME_0 <- admin0$shapeName
  admin1$NAME_1 <- admin1$shapeName
  admin2$NAME_2 <- admin2$shapeName

  if (sf::st_crs(admin1) != sf::st_crs(admin2)) {
    admin1 <- sf::st_transform(admin1, sf::st_crs(admin2))
  }

  # attach exactly one parent NAME_1 per admin2 feature
  idx_list <- sf::st_intersects(admin2, admin1)

  admin2$NAME_1 <- vapply(idx_list, function(ii) {
    if (length(ii) == 0) NA_character_ else admin1$NAME_1[ii[1]]
  }, character(1))

  admin2 <- admin2 %>%
    dplyr::mutate(
      NAME_2_UID = paste(.data$NAME_1, .data$NAME_2, sep = " | "),
      admin2.name.full = paste0(.data$NAME_1, "_", .data$NAME_2_UID)
    )

  stopifnot(anyDuplicated(admin2$NAME_2_UID) == 0)
  stopifnot(anyDuplicated(admin2$admin2.name.full) == 0)

  geo <- .read_gps_points(gps_path, sf::st_crs(admin1))

  infos <- .build_cluster_admin_info(
    geo,
    admin1,
    admin2,
    by_adm2 = "NAME_2_UID",
    map = TRUE
  )

  res <- c(
    infos,
    list(
      admin0 = admin0,
      admin1 = admin1,
      admin2 = admin2,
      geo = geo
    )
  )

  .save_geo_result(res, out_path)
  res
}
