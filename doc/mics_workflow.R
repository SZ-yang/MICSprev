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
# install.packages("devtools")
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
# bh_path  <- system.file("extdata", "bh.sav", package = "MICSprev")
# wm_path  <- system.file("extdata", "wm.sav", package = "MICSprev")
# ch_path  <- system.file("extdata", "ch.sav", package = "MICSprev")
# gps_path <- system.file(
#   "extdata",
#   "NigeriaMICS2021GPS",
#   "NigeriaMICS2021GPS.shp",
#   package = "MICSprev"
# )

# data_dir <- "C:/Users/zayny/Documents/PER/OneDrive/files/MICS/RawData/Nigeria MICS6 Datasets/Nigeria MICS6 SPSS Datasets"
# gps_dir  <- "C:/Users/zayny/Documents/PER/OneDrive/files/MICS/Rawdata/NigeriaMICS2021GPS"

data_dir <- "C:/Users/Z/OneDrive/files/MICS/RawData/Nigeria MICS6 Datasets/Nigeria MICS6 SPSS Datasets"
gps_dir  <- "C:/Users/Z/OneDrive/files/MICS/Rawdata/NigeriaMICS2021GPS"

# Survey files
bh_path <- file.path(data_dir, "bh.sav")
wm_path <- file.path(data_dir, "wm.sav")
ch_path <- file.path(data_dir, "ch.sav")

# GPS shapefile
gps_path <- file.path(gps_dir, "NigeriaMICS2021GPS.shp")


## ----eval=FALSE---------------------------------------------------------------
# # Build geographic inputs from the GPS shapefile
# geo <- process_geo_mics(
#   country.name = "Nigeria",
#   year = 2021,
#   gps_path = gps_path
# )

## ----echo = FALSE-------------------------------------------------------------
geo <- process_geo_mics(
  country.name = "Nigeria",
  year = 2021,
  gps_path = gps_path
)

## -----------------------------------------------------------------------------
str(geo, max.level = 1)

## -----------------------------------------------------------------------------
bh <- read_sav(bh_path)
wm <- read_sav(wm_path)
ch <- read_sav(ch_path)

## ----eval=FALSE---------------------------------------------------------------
# dat_nmr  <- process_indicator(bh, "NMR")
# dat_anc  <- process_indicator(wm, "ANC")
# dat_dtp3 <- process_indicator(ch, "DTP3")

## ----echo=FALSE---------------------------------------------------------------
dat_nmr  <- process_NMR(bh)
dat_anc  <- process_ANC(wm)
dat_dtp3 <- process_DTP3(ch)

## -----------------------------------------------------------------------------
lapply(list(NMR = dat_nmr, ANC = dat_anc, DTP3 = dat_dtp3), function(x) {
  list(
    cols_ok   = all(c("cluster", "householdID", "weight", "strata", "value") %in% names(x)),
    n         = nrow(x),
    value_tab = table(x$value, useNA = "ifany")
  )
})

## -----------------------------------------------------------------------------
fit_nmr <- run_indicator_models(
  data = dat_nmr,
  geo = geo,
  admin_levels = c(0, 1, 2),
  models = c("direct", "fh", "cluster", "cluster_strat"),
  model = "bym2",
  ci = 0.95,
  fix_var = TRUE,
  nested = TRUE,
  verbose = TRUE
)

## -----------------------------------------------------------------------------
fit_anc <- run_indicator_models(
  data = dat_anc,
  geo  = geo,
  admin_levels = c(1, 2),
  models = c("direct", "cluster_strat"),
  verbose = TRUE
)

## -----------------------------------------------------------------------------
fit_dtp3 <- run_indicator_models(
  data = dat_dtp3,
  geo  = geo,
  admin_levels = c(0, 1, 2),
  models = c("direct"),
  verbose = TRUE
)

## -----------------------------------------------------------------------------
nmr_nat <- fit_nmr$fits[["0"]]$direct$res.natl$direct.est
cat("National NMR per 1,000:", round(nmr_nat * 1000, 2), "\n")

## -----------------------------------------------------------------------------
if ("0" %in% names(fit_anc$fits) && !is.null(fit_anc$fits[["0"]]$direct)) {
  anc_nat <- fit_anc$fits[["0"]]$direct$res.admin0$direct.est
  cat("National ANC coverage:", round(anc_nat, 4), "\n")
} else {
  cat("ANC national estimate not available because admin 0 was not requested.\n")
}

## -----------------------------------------------------------------------------
nmr_adm1_direct <- fit_nmr$fits[["1"]]$direct$res.admin1
head(nmr_adm1_direct)

## -----------------------------------------------------------------------------
nmr_adm2_fh <- fit_nmr$fits[["2"]]$fh$res.admin2
head(nmr_adm2_fh)

## -----------------------------------------------------------------------------
anc_res_admin1 <- fit_anc$fits[["1"]]$cluster_strat$res.admin1

anc_adm1_strat_full <- if ("type" %in% names(anc_res_admin1)) {
  subset(anc_res_admin1, type == "full")
} else {
  anc_res_admin1
}

head(anc_adm1_strat_full)

## -----------------------------------------------------------------------------
out_dir <- tempfile("MICSprev_plots_")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ----fig.width=13, fig.height=8-----------------------------------------------
nmr_maps_adm1 <- plot_indicator_maps(
  fit    = fit_nmr,
  geo    = geo,
  admin  = 1,
  models = c("direct", "fh", "cluster", "cluster_strat"),
  scale  = 1000,
  direction = -1,
  out_dir = out_dir,
  prefix  = "NMR_NGA2021"
)

nmr_maps_adm1$mean_map
nmr_maps_adm1$cv_map

## ----fig.width=13, fig.height=8-----------------------------------------------
nmr_maps_adm2 <- plot_indicator_maps(
  fit = fit_nmr,
  geo = geo,
  admin = 2,
  models = c("direct", "fh", "cluster", "cluster_strat"),
  scale = 1000,
  direction = -1,
  out_dir = out_dir,
  prefix = "NMR_NGA2021"
)

nmr_maps_adm2$mean_map
nmr_maps_adm2$cv_map


nmr_maps_adm2_log <- plot_indicator_maps(
  fit = fit_nmr,
  geo = geo,
  admin = 2,
  models = c("direct", "fh", "cluster", "cluster_strat"),
  scale = 1000,
  direction = -1,
  out_dir = out_dir,
  prefix = "NMR_NGA2021_log1p",
  transform = "log1p"
)

nmr_maps_adm2_log$mean_map
nmr_maps_adm2_log$cv_map

## ----fig.width=11, fig.height=16----------------------------------------------
p_nmr_ridge_adm1 <- plot_indicator_ridge(
  fit = fit_nmr,
  admin = 1,
  models = c("fh", "cluster", "cluster_strat"),
  threshold = 12 / 1000,
  scale = 1000,
  out_dir = out_dir,
  prefix = "NMR_NGA2021"
)
p_nmr_ridge_adm1

## ----fig.width=10, fig.height=7-----------------------------------------------
anc_maps_adm1 <- plot_indicator_maps(
  fit = fit_anc,
  geo = geo,
  admin = 1,
  models = c("direct", "cluster_strat"),
  scale = 100,
  direction = 1,
  out_dir = out_dir,
  prefix = "ANC4_NGA2021"
)

anc_maps_adm1$mean_map
anc_maps_adm1$cv_map

p_anc_ridge_adm1 <- plot_indicator_ridge(
  fit = fit_anc,
  admin = 1,
  models = c("cluster_strat"),
  threshold = NULL,
  scale = 100,
  out_dir = out_dir,
  prefix = "ANC4_NGA2021"
)
p_anc_ridge_adm1

## ----fig.width=10, fig.height=7-----------------------------------------------
dtp3_maps_adm1 <- plot_indicator_maps(
  fit = fit_dtp3,
  geo = geo,
  admin = 1,
  models = c("direct"),
  scale = 100,
  direction = 1,
  out_dir = out_dir,
  prefix = "DTP3_NGA2021"
)

dtp3_maps_adm1$mean_map
dtp3_maps_adm1$cv_map

## ----surveyprev-direct, fig.width=10, fig.height=6----------------------------
# Load the underlying modeling engine
library(surveyPrev)

# 1. Calculate direct estimates at Admin 1
direct_res <- directEST(
  data = dat_anc,
  cluster.info = geo$cluster.info,
  admin = 1,
  aggregation = FALSE
)

# 2. Fit a Fay-Herriot (BYM2) model at Admin 1
fh_res <- fhModel(
  data = dat_anc,
  cluster.info = geo$cluster.info,
  admin.info = geo$admin.info1,
  admin = 1,
  model = "bym2",
  aggregation = FALSE
)

# 3. Use surveyPrev's built-in plotting functions
# Compare the smoothed FH estimates against the Direct estimates
scatterPlot(
  res1 = direct_res$res.admin1,
  res2 = fh_res$res.admin1,
  value1 = "direct.est",
  value2 = "mean",
  by.res1 = "admin1.name",
  by.res2 = "admin1.name",
  title = "ANC Admin 1: Direct vs FH Estimates",
  label1 = "Direct Estimate",
  label2 = "Smoothed FH Estimate"
)

## -----------------------------------------------------------------------------
sessionInfo()

