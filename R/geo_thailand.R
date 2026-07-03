# R/geo_thailand.R
# Thailand-specific geography helpers (internal)

# Internal
.tha_expected_admin1_names <- function() {
  c("Bangkok", "Central", "North", "Northeast", "South")
}

# Internal
.tha_admin2_name_aliases <- function() {
  data.frame(
    geoboundaries_name = c(
      "Buri Ram",
      "Lopburi",
      "Nong Bua Lam Phu",
      "Phra Nakhon Si Ayutthaya",
      "Saraburi",
      "Si Sa Ket"
    ),
    admin2_name = c(
      "Buriram",
      "Lop Buri",
      "Nong Bua Lamphu",
      "Ayutthaya",
      "Sara Buri",
      "Sisaket"
    ),
    stringsAsFactors = FALSE
  )
}

# Internal
.tha_normalize_admin2_name <- function(x) {
  x <- .geo_clean_surveyprev_name(x)
  x <- stringr::str_remove(
    x,
    stringr::regex("\\s+Province$", ignore_case = TRUE)
  )
  x <- stringr::str_squish(x)

  aliases <- .tha_admin2_name_aliases()
  hit <- match(x, aliases$geoboundaries_name)
  x[!is.na(hit)] <- aliases$admin2_name[hit[!is.na(hit)]]
  x
}

# Internal
.tha_work_crs <- function() {
  3857
}

# Internal
.tha_keep_polygon_geometry <- function(x, label) {
  x <- sf::st_make_valid(x)
  x <- suppressWarnings(sf::st_collection_extract(x, "POLYGON", warn = FALSE))
  x <- sf::st_make_valid(x)
  .tha_validate_polygon_geometry(x, label)
  x
}

# Internal
.tha_validate_polygon_geometry <- function(x, label) {
  geom_type <- as.character(sf::st_geometry_type(x))
  bad_type <- !geom_type %in% c("POLYGON", "MULTIPOLYGON")

  if (any(bad_type)) {
    stop(label, " must contain only polygon geometries.", call. = FALSE)
  }

  empty <- sf::st_is_empty(x)
  if (any(is.na(empty) | empty)) {
    stop(label, " contains empty geometries.", call. = FALSE)
  }

  bbox <- sf::st_bbox(x)
  if (any(is.na(bbox) | !is.finite(bbox))) {
    stop(label, " has a non-finite bounding box.", call. = FALSE)
  }

  invisible(x)
}

# Internal
.tha_admin2_admin1_crosswalk <- function() {
  crosswalk <- data.frame(
    admin2_name = c(
      "Bangkok",
      "Ang Thong",
      "Ayutthaya",
      "Chachoengsao",
      "Chai Nat",
      "Chanthaburi",
      "Chon Buri",
      "Kanchanaburi",
      "Lop Buri",
      "Nakhon Nayok",
      "Nakhon Pathom",
      "Nonthaburi",
      "Pathum Thani",
      "Phetchaburi",
      "Prachin Buri",
      "Prachuap Khiri Khan",
      "Ratchaburi",
      "Rayong",
      "Sa Kaeo",
      "Samut Prakan",
      "Samut Sakhon",
      "Samut Songkhram",
      "Sara Buri",
      "Sing Buri",
      "Suphan Buri",
      "Trat",
      "Chiang Mai",
      "Chiang Rai",
      "Kamphaeng Phet",
      "Lampang",
      "Lamphun",
      "Mae Hong Son",
      "Nakhon Sawan",
      "Nan",
      "Phayao",
      "Phetchabun",
      "Phichit",
      "Phitsanulok",
      "Phrae",
      "Sukhothai",
      "Tak",
      "Uthai Thani",
      "Uttaradit",
      "Amnat Charoen",
      "Bueng Kan",
      "Buriram",
      "Chaiyaphum",
      "Kalasin",
      "Khon Kaen",
      "Loei",
      "Maha Sarakham",
      "Mukdahan",
      "Nakhon Phanom",
      "Nakhon Ratchasima",
      "Nong Bua Lamphu",
      "Nong Khai",
      "Roi Et",
      "Sakon Nakhon",
      "Sisaket",
      "Surin",
      "Ubon Ratchathani",
      "Udon Thani",
      "Yasothon",
      "Chumphon",
      "Krabi",
      "Nakhon Si Thammarat",
      "Narathiwat",
      "Pattani",
      "Phangnga",
      "Phatthalung",
      "Phuket",
      "Ranong",
      "Satun",
      "Songkhla",
      "Surat Thani",
      "Trang",
      "Yala"
    ),
    admin1_name = c(
      "Bangkok",
      rep("Central", 25),
      rep("North", 17),
      rep("Northeast", 20),
      rep("South", 14)
    ),
    stringsAsFactors = FALSE
  )

  crosswalk$admin2_name <- .geo_clean_surveyprev_name(crosswalk$admin2_name)
  crosswalk$admin1_name <- .geo_clean_surveyprev_name(crosswalk$admin1_name)
  crosswalk$admin2_key <- .tha_normalize_admin2_name(crosswalk$admin2_name)

  .tha_validate_crosswalk(crosswalk)
  crosswalk
}

# Internal
.tha_validate_crosswalk <- function(crosswalk) {
  required <- c("admin2_name", "admin1_name", "admin2_key")
  missing_cols <- setdiff(required, names(crosswalk))

  if (length(missing_cols) > 0) {
    stop(
      "Thailand crosswalk is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(crosswalk) != 77L) {
    stop("Thailand crosswalk must contain 77 provinces.", call. = FALSE)
  }

  blank_admin2 <- is.na(crosswalk$admin2_name) | trimws(crosswalk$admin2_name) == ""
  blank_admin1 <- is.na(crosswalk$admin1_name) | trimws(crosswalk$admin1_name) == ""

  if (any(blank_admin2) || any(blank_admin1)) {
    stop("Thailand crosswalk contains missing or blank names.", call. = FALSE)
  }

  dup_admin2 <- unique(crosswalk$admin2_name[duplicated(crosswalk$admin2_name)])
  if (length(dup_admin2) > 0) {
    stop(
      "Thailand crosswalk contains duplicated province name(s): ",
      paste(dup_admin2, collapse = ", "),
      call. = FALSE
    )
  }

  dup_key <- unique(crosswalk$admin2_key[duplicated(crosswalk$admin2_key)])
  if (length(dup_key) > 0) {
    stop(
      "Thailand crosswalk contains duplicated normalized province key(s): ",
      paste(dup_key, collapse = ", "),
      call. = FALSE
    )
  }

  expected <- sort(.tha_expected_admin1_names())
  got <- sort(unique(crosswalk$admin1_name))

  if (!identical(got, expected)) {
    stop(
      "Thailand crosswalk admin1 groups must be exactly: ",
      paste(expected, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(crosswalk)
}

# Internal
.tha_stop_unmatched_admin2 <- function(geo_names, crosswalk_names) {
  unmatched_geo <- setdiff(geo_names, crosswalk_names)
  unmatched_crosswalk <- setdiff(crosswalk_names, geo_names)

  if (length(unmatched_geo) == 0 && length(unmatched_crosswalk) == 0) {
    return(invisible(NULL))
  }

  stop(
    "Thailand GeoBoundaries province names do not match the internal crosswalk.\n",
    "Unmatched GeoBoundaries names: ",
    if (length(unmatched_geo) == 0) "<none>" else paste(unmatched_geo, collapse = ", "),
    "\nUnmatched crosswalk names: ",
    if (length(unmatched_crosswalk) == 0) "<none>" else paste(unmatched_crosswalk, collapse = ", "),
    call. = FALSE
  )
}

# Internal
.tha_prepare_admin2 <- function(admin2) {
  if (!inherits(admin2, "sf")) {
    stop("Thailand admin2 input must be an sf object.", call. = FALSE)
  }

  name_col <- .geo_first_existing_col(
    admin2,
    c("shapeName", "NAME_2"),
    label = "Thailand province name column"
  )

  admin2 <- admin2 |>
    sf::st_make_valid() |>
    dplyr::mutate(
      NAME_2_SOURCE = .geo_clean_surveyprev_name(.data[[name_col]]),
      admin2_key = .tha_normalize_admin2_name(.data$NAME_2_SOURCE)
    )

  crosswalk <- .tha_admin2_admin1_crosswalk()

  .tha_stop_unmatched_admin2(
    geo_names = sort(unique(admin2$admin2_key)),
    crosswalk_names = sort(unique(crosswalk$admin2_key))
  )

  admin2 <- dplyr::left_join(
    admin2,
    crosswalk,
    by = "admin2_key"
  ) |>
    dplyr::mutate(
      NAME_2 = as.character(.data$admin2_name),
      NAME_1 = as.character(.data$admin1_name),
      admin2.name.full = paste0(.data$NAME_1, "_", .data$NAME_2)
    ) |>
    dplyr::select(-dplyr::all_of(c("admin1_name", "admin2_name", "admin2_key"))) |>
    sf::st_make_valid()

  admin2 <- .tha_keep_polygon_geometry(admin2, "Thailand admin2")
  .tha_validate_admin2(admin2)
  admin2
}

# Internal
.tha_validate_admin2 <- function(admin2) {
  required <- c("NAME_1", "NAME_2", "admin2.name.full")
  missing_cols <- setdiff(required, names(admin2))

  if (length(missing_cols) > 0) {
    stop(
      "Thailand admin2 is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (length(unique(admin2$NAME_2)) != 77L) {
    stop("Thailand admin2 must contain 77 unique provinces.", call. = FALSE)
  }

  has_missing <- is.na(admin2$NAME_1) |
    is.na(admin2$NAME_2) |
    is.na(admin2$admin2.name.full) |
    trimws(admin2$NAME_1) == "" |
    trimws(admin2$NAME_2) == "" |
    trimws(admin2$admin2.name.full) == ""

  if (any(has_missing)) {
    stop("Thailand admin2 contains missing admin names.", call. = FALSE)
  }

  if (anyDuplicated(admin2$NAME_2) > 0) {
    stop("Thailand admin2 province names must be unique.", call. = FALSE)
  }

  if (anyDuplicated(admin2$admin2.name.full) > 0) {
    stop("Thailand admin2 full names must be unique.", call. = FALSE)
  }

  got <- sort(unique(admin2$NAME_1))
  expected <- sort(.tha_expected_admin1_names())

  if (!identical(got, expected)) {
    stop(
      "Thailand admin2 must map to admin1 groups: ",
      paste(expected, collapse = ", "),
      call. = FALSE
    )
  }

  invalid <- is.na(sf::st_is_valid(admin2)) | !sf::st_is_valid(admin2)
  if (any(invalid)) {
    stop("Thailand admin2 contains invalid geometries.", call. = FALSE)
  }

  .tha_validate_polygon_geometry(admin2, "Thailand admin2")

  invisible(admin2)
}

# Internal
.tha_build_admin1_from_admin2 <- function(admin2) {
  .tha_validate_admin2(admin2)

  admin1 <- .tha_dissolve_admin1(admin2)

  .tha_validate_admin1(admin1)
  admin1
}

# Internal
.tha_dissolve_admin1 <- function(admin2) {
  original_crs <- sf::st_crs(admin2)

  admin2_work <- admin2 |>
    dplyr::select(dplyr::all_of("NAME_1")) |>
    .tha_keep_polygon_geometry("Thailand admin2")

  if (!is.na(original_crs)) {
    admin2_work <- sf::st_transform(admin2_work, .tha_work_crs())
  }

  admin1 <- admin2_work |>
    dplyr::group_by(.data$NAME_1) |>
    dplyr::summarise(.groups = "drop") |>
    .tha_keep_polygon_geometry("Thailand admin1")

  if (!is.na(original_crs)) {
    admin1 <- sf::st_transform(admin1, original_crs)
  }

  .tha_keep_polygon_geometry(admin1, "Thailand admin1")
}

# Internal
.tha_validate_admin1 <- function(admin1) {
  if (!inherits(admin1, "sf")) {
    stop("Thailand admin1 must be an sf object.", call. = FALSE)
  }

  if (!"NAME_1" %in% names(admin1)) {
    stop("Thailand admin1 must contain NAME_1.", call. = FALSE)
  }

  got <- sort(unique(admin1$NAME_1))
  expected <- sort(.tha_expected_admin1_names())

  if (!identical(got, expected) || nrow(admin1) != length(expected)) {
    stop(
      "Thailand admin1 must contain exactly: ",
      paste(expected, collapse = ", "),
      call. = FALSE
    )
  }

  invalid <- is.na(sf::st_is_valid(admin1)) | !sf::st_is_valid(admin1)
  if (any(invalid)) {
    stop("Thailand admin1 contains invalid geometries.", call. = FALSE)
  }

  .tha_validate_polygon_geometry(admin1, "Thailand admin1")

  invisible(admin1)
}

# Internal
.tha_snap_geo_to_admin2 <- function(geo, admin2) {
  if (!inherits(geo, "sf")) {
    stop("Thailand GPS input must be an sf object.", call. = FALSE)
  }

  if (!inherits(admin2, "sf")) {
    stop("Thailand admin2 input must be an sf object.", call. = FALSE)
  }

  admin2_for_join <- admin2
  if (sf::st_crs(admin2_for_join) != sf::st_crs(geo)) {
    admin2_for_join <- sf::st_transform(admin2_for_join, sf::st_crs(geo))
  }

  hit_admin2 <- lengths(sf::st_intersects(geo, admin2_for_join)) > 0
  has_geom <- !sf::st_is_empty(geo)
  outside_idx <- which(has_geom & !hit_admin2)

  geo$geo_fix_type <- "original"

  if (length(outside_idx) == 0) {
    return(geo)
  }

  geo_work <- sf::st_transform(geo, .tha_work_crs())
  admin2_work <- admin2_for_join |>
    sf::st_transform(.tha_work_crs()) |>
    .tha_keep_polygon_geometry("Thailand admin2")

  admin2_points <- admin2_work |>
    dplyr::select(dplyr::all_of(c("NAME_1", "NAME_2")))

  sf::st_geometry(admin2_points) <-
    sf::st_point_on_surface(sf::st_geometry(admin2_work))

  nearest_idx <- sf::st_nearest_feature(
    geo_work[outside_idx, ],
    admin2_work
  )

  sf::st_geometry(geo_work)[outside_idx] <-
    sf::st_geometry(admin2_points[nearest_idx, ])

  geo_work$geo_fix_type[outside_idx] <- "outside_coord_to_nearest_province"

  geo_fixed <- sf::st_transform(geo_work, sf::st_crs(geo))
  xy_fixed <- .geo_coords_xy_safe(geo_fixed)
  geo_fixed$LONGNUM <- xy_fixed[, 1]
  geo_fixed$LATNUM <- xy_fixed[, 2]

  geo_fixed
}

# Internal
.tha_prepare_admin3 <- function(admin3, admin2) {
  if (is.null(admin3) || !inherits(admin3, "sf")) {
    return(NULL)
  }

  tryCatch(
    {
      name_col <- .geo_first_existing_col(
        admin3,
        c("shapeName", "NAME_3"),
        required = FALSE,
        label = "Thailand lower-level name column"
      )

      if (is.na(name_col)) {
        return(NULL)
      }

      admin3 <- admin3 |>
        sf::st_make_valid() |>
        dplyr::mutate(
          NAME_3 = .geo_clean_surveyprev_name(.data[[name_col]])
        ) |>
        dplyr::select(-dplyr::any_of(c("NAME_1", "NAME_2", "admin3.name.full")))

      admin2_parent <- admin2 |>
        dplyr::select(dplyr::all_of(c("NAME_1", "NAME_2")))

      if (sf::st_crs(admin2_parent) != sf::st_crs(admin3)) {
        admin2_parent <- sf::st_transform(admin2_parent, sf::st_crs(admin3))
      }

      admin3 <- sf::st_join(
        admin3,
        admin2_parent,
        join = sf::st_intersects,
        largest = TRUE
      ) |>
        dplyr::mutate(
          admin3.name.full = paste0(.data$NAME_1, "_", .data$NAME_2, "_", .data$NAME_3)
        )

      required <- c("NAME_1", "NAME_2", "NAME_3", "admin3.name.full")
      if (!all(required %in% names(admin3))) {
        return(NULL)
      }

      has_missing <- is.na(admin3$NAME_1) |
        is.na(admin3$NAME_2) |
        is.na(admin3$NAME_3) |
        trimws(admin3$NAME_1) == "" |
        trimws(admin3$NAME_2) == "" |
        trimws(admin3$NAME_3) == ""

      if (any(has_missing)) {
        return(NULL)
      }

      admin3
    },
    error = function(e) NULL
  )
}
