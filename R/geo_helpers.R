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

#' @keywords internal
.geo_clean_surveyprev_name <- function(x) {
  x |>
    as.character() |>
    stringr::str_replace_all("_+", " ") |>
    stringr::str_squish()
}

#' @keywords internal
.geo_first_existing_col <- function(x, candidates, required = TRUE, label = "column") {
  hit <- candidates[candidates %in% names(x)]

  if (length(hit) == 0 && required) {
    stop(
      "Could not find ", label, ". Tried: ",
      paste(candidates, collapse = ", "),
      "\nAvailable columns: ",
      paste(names(x), collapse = ", "),
      call. = FALSE
    )
  }

  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

#' @keywords internal
.geo_coords_xy_safe <- function(x) {
  g <- sf::st_geometry(x)
  out <- matrix(NA_real_, nrow = length(g), ncol = 2)
  colnames(out) <- c("X", "Y")

  non_empty <- !sf::st_is_empty(g)

  if (any(non_empty)) {
    cc <- sf::st_coordinates(g[non_empty])
    out[non_empty, 1] <- cc[, 1]
    out[non_empty, 2] <- cc[, 2]
  }

  out
}

#' @keywords internal
.geo_pick_ocha_shp <- function(shp_files, level) {
  base <- basename(shp_files)

  pat <- paste0("(^|[_-])(adm|admin)", level, "($|[_-]|\\.)")

  hit <- shp_files[stringr::str_detect(base, stringr::regex(pat, ignore_case = TRUE))]

  hit <- hit[
    !stringr::str_detect(
      basename(hit),
      stringr::regex("line|lines|admln|point|points", ignore_case = TRUE)
    )
  ]

  if (length(hit) == 0) {
    stop(
      "Could not find OCHA admin", level, " polygon shapefile.\n",
      "Available shapefiles:\n",
      paste(basename(shp_files), collapse = "\n"),
      call. = FALSE
    )
  }

  hit[1]
}

#' @keywords internal
.geo_download_mwi_ocha_boundaries <- function(
    download_dir = file.path(tempdir(), "mwi_ocha_boundaries"),
    force_download = FALSE
) {
  if (!requireNamespace("rhdx", quietly = TRUE)) {
    stop(
      "Package `rhdx` is required to download Malawi OCHA boundaries. ",
      "Install it with remotes::install_github('dickoa/rhdx').",
      call. = FALSE
    )
  }

  dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)

  rhdx::set_rhdx_config(hdx_site = "prod", read_only = TRUE)

  ds <- rhdx::pull_dataset("cod-ab-mwi")
  resources <- rhdx::get_resources(ds)

  res_text <- vapply(
    resources,
    function(x) paste(capture.output(print(x)), collapse = " "),
    character(1)
  )

  res_fmt <- vapply(
    resources,
    function(x) {
      tryCatch(rhdx::get_resource_format(x), error = function(e) NA_character_)
    },
    character(1)
  )

  search_text <- paste(res_fmt, res_text)

  score <- rep(0, length(resources))
  score <- score + 5 * stringr::str_detect(
    search_text,
    stringr::regex("\\bSHP\\b|shp|shape|shapefile", ignore_case = TRUE)
  )
  score <- score + 4 * stringr::str_detect(
    search_text,
    stringr::regex("boundary|boundaries|admin|admbnda|COD-AB", ignore_case = TRUE)
  )
  score <- score + 2 * stringr::str_detect(
    search_text,
    stringr::regex("mwi|malawi", ignore_case = TRUE)
  )
  score <- score + 1 * stringr::str_detect(
    search_text,
    stringr::regex("\\.zip", ignore_case = TRUE)
  )

  if (all(score <= 0)) {
    stop("Could not identify a shapefile-like OCHA Malawi resource.", call. = FALSE)
  }

  boundary_resource <- resources[[which.max(score)]]

  before_files <- list.files(download_dir, recursive = TRUE, full.names = TRUE)

  invisible(
    rhdx::download_resource(
      boundary_resource,
      folder = download_dir,
      quiet = TRUE,
      force = force_download
    )
  )

  after_files <- list.files(download_dir, recursive = TRUE, full.names = TRUE)
  new_files <- setdiff(after_files, before_files)
  candidate_files <- if (length(new_files) > 0) new_files else after_files

  zip_files <- candidate_files[
    stringr::str_detect(candidate_files, stringr::regex("\\.zip$", ignore_case = TRUE))
  ]

  if (length(zip_files) == 0) {
    zip_files <- after_files[
      stringr::str_detect(after_files, stringr::regex("\\.zip$", ignore_case = TRUE))
    ]
  }

  unzip_dir <- file.path(download_dir, "unzipped")
  dir.create(unzip_dir, recursive = TRUE, showWarnings = FALSE)

  if (length(zip_files) > 0) {
    zip_files <- zip_files[order(file.info(zip_files)$mtime, decreasing = TRUE)]
    utils::unzip(zip_files[1], exdir = unzip_dir)
    shp_root <- unzip_dir
  } else {
    shp_root <- download_dir
  }

  shp_files <- list.files(
    shp_root,
    pattern = "\\.shp$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(shp_files) == 0) {
    stop("No .shp files found after downloading OCHA Malawi boundaries.", call. = FALSE)
  }

  shp_admin0 <- .geo_pick_ocha_shp(shp_files, level = 0)
  shp_admin2 <- .geo_pick_ocha_shp(shp_files, level = 2)
  shp_admin3 <- .geo_pick_ocha_shp(shp_files, level = 3)

  admin0 <- sf::st_read(shp_admin0, quiet = TRUE) |>
    sf::st_make_valid()

  admin1 <- sf::st_read(shp_admin2, quiet = TRUE) |>
    sf::st_make_valid()

  admin2 <- sf::st_read(shp_admin3, quiet = TRUE) |>
    sf::st_make_valid()

  adm0_col <- .geo_first_existing_col(
    admin0,
    c("adm0_name", "ADM0_EN", "ADM0_NAME", "admin0Name", "NAME_0"),
    required = FALSE,
    label = "admin0 name column"
  )

  admin0 <- admin0 |>
    dplyr::mutate(
      NAME_0 = if (!is.na(adm0_col)) {
        .geo_clean_surveyprev_name(.data[[adm0_col]])
      } else {
        "Malawi"
      }
    )

  adm2_name_col <- .geo_first_existing_col(
    admin1,
    c("adm2_name", "ADM2_EN", "ADM2_NAME", "admin2Name", "NAME_1"),
    label = "OCHA admin2 district name column"
  )

  admin1 <- admin1 |>
    dplyr::mutate(
      NAME_1 = .geo_clean_surveyprev_name(.data[[adm2_name_col]])
    )

  adm3_parent_col <- .geo_first_existing_col(
    admin2,
    c("adm2_name", "ADM2_EN", "ADM2_NAME", "admin2Name", "NAME_1"),
    label = "OCHA admin3 parent district column"
  )

  adm3_name_col <- .geo_first_existing_col(
    admin2,
    c("adm3_name", "ADM3_EN", "ADM3_NAME", "admin3Name", "NAME_2_raw"),
    label = "OCHA admin3 TA name column"
  )

  admin2 <- admin2 |>
    dplyr::mutate(
      NAME_1 = .geo_clean_surveyprev_name(.data[[adm3_parent_col]]),
      NAME_2_raw = .geo_clean_surveyprev_name(.data[[adm3_name_col]]),

      # Important:
      # NAME_2 is TA name only. surveyPrev will create admin2.name.full as
      # paste(admin1, admin2, sep = "_").
      NAME_2 = .data$NAME_2_raw,
      NAME_2_FULL_CHECK = paste(.data$NAME_1, .data$NAME_2, sep = "_"),
      admin2.name.full = .data$NAME_2_FULL_CHECK
    )

  list(
    admin0 = admin0,
    admin1 = admin1,
    admin2 = admin2
  )
}
