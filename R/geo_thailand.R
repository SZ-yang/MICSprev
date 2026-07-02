# R/geo_thailand.R
# Thailand-specific geography helpers (internal)

# Internal
.tha_expected_admin1_names <- function() {
  c("Bangkok", "Central", "North", "Northeast", "South")
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

  .tha_validate_crosswalk(crosswalk)
  crosswalk
}

# Internal
.tha_validate_crosswalk <- function(crosswalk) {
  required <- c("admin2_name", "admin1_name")
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
      NAME_2 = .geo_clean_surveyprev_name(.data[[name_col]])
    )

  crosswalk <- .tha_admin2_admin1_crosswalk()

  .tha_stop_unmatched_admin2(
    geo_names = sort(unique(admin2$NAME_2)),
    crosswalk_names = sort(unique(crosswalk$admin2_name))
  )

  admin2 <- dplyr::left_join(
    admin2,
    crosswalk,
    by = c("NAME_2" = "admin2_name")
  ) |>
    dplyr::mutate(
      NAME_1 = as.character(.data$admin1_name),
      admin2.name.full = paste0(.data$NAME_1, "_", .data$NAME_2)
    ) |>
    dplyr::select(-dplyr::all_of("admin1_name")) |>
    sf::st_make_valid()

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
  admin2 |>
    dplyr::select(dplyr::all_of("NAME_1")) |>
    dplyr::group_by(.data$NAME_1) |>
    dplyr::summarise(.groups = "drop") |>
    sf::st_make_valid()
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

  invisible(admin1)
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
