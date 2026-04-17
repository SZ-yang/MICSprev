## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = TRUE,
  message = FALSE,
  warning = FALSE,
  fig.width = 10,
  fig.height = 6.5,
  fig.align = "center",
  out.width = "100%",
  dpi = 144
)

## ----eval=FALSE---------------------------------------------------------------
# install.packages(c(
#   "sf",
#   "dplyr",
#   "ggplot2",
#   "haven",
#   "patchwork",
#   "stringr",
#   "lwgeom"
# ))

## ----eval=FALSE---------------------------------------------------------------
# install.packages("remotes")
# 
# remotes::install_github("richardli/SUMMER")
# remotes::install_github("richardli/surveyPrev")

## ----eval=FALSE---------------------------------------------------------------
# install.packages(
#   "INLA",
#   repos = c(
#     getOption("repos"),
#     INLA = "https://inla.r-inla-download.org/R/stable"
#   ),
#   dep = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# # install.packages("devtools")
# devtools::install_github("SZ-yang/MICSprev")

## -----------------------------------------------------------------------------
library(MICSprev)
library(haven)
library(dplyr)

## -----------------------------------------------------------------------------
ls("package:MICSprev")
?process_geo_mics

## ----eval=FALSE---------------------------------------------------------------
# NigeriaMICS2021GPS/
#   NigeriaMICS2021GPS.shp
#   NigeriaMICS2021GPS.shx
#   NigeriaMICS2021GPS.dbf
#   NigeriaMICS2021GPS.prj
#   NigeriaMICS2021GPS.cpg
#   ...

## ----eval=FALSE---------------------------------------------------------------
# # Set paths to the unzipped survey data folder and GPS folder
# data_dir <- "path/to/unzipped/MICS_survey_data"
# gps_dir  <- "path/to/unzipped/NigeriaMICS2021GPS"
# 
# # Survey files
# bh_path <- file.path(data_dir, "bh.sav")
# wm_path <- file.path(data_dir, "wm.sav")
# ch_path <- file.path(data_dir, "ch.sav")
# 
# # GPS shapefile
# gps_path <- file.path(gps_dir, "NigeriaMICS2021GPS.shp")

## ----echo = FALSE-------------------------------------------------------------
bh_path  <- system.file("extdata", "bh.sav", package = "MICSprev")
wm_path  <- system.file("extdata", "wm.sav", package = "MICSprev")
ch_path  <- system.file("extdata", "ch.sav", package = "MICSprev")
gps_path <- system.file(
  "extdata",
  "NigeriaMICS2021GPS",
  "NigeriaMICS2021GPS.shp",
  package = "MICSprev"
)

## ----eval=FALSE---------------------------------------------------------------
# # Build geographic inputs from the GPS shapefile
# geo <- process_geo_mics(
#   country.name = "Nigeria",
#   year = 2021,
#   gps_path = gps_path
# )

