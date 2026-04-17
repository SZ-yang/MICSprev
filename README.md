# MICSprev

**MICSprev** is an R package for processing **Multiple Indicator Cluster Survey (MICS)** data and fitting spatial small-area estimation models using the **surveyPrev** and **SUMMER** frameworks.

The package currently supports a full workflow for:

- processing MICS survey files into standardized indicator-ready datasets
- building cluster and administrative geographic objects from MICS GPS shapefiles
- fitting direct and model-based estimators at national, admin1, and admin2 levels
- comparing multiple model families, including two Fay-Herriot strategies
- generating publication-ready maps and ridge plots

## Features

- **Built-in indicator processing**
  - Neonatal mortality rate (**NMR**)
  - Antenatal care 4+ visits (**ANC**)
  - DTP3 vaccination coverage (**DTP3**)

- **Geographic preprocessing**
  - Converts GPS shapefiles into the `cluster.info`, `admin.info1`, `admin.info2`, `admin0`, `admin1`, `admin2`, and `geo` objects needed downstream

- **Model fitting wrapper**
  - `direct` direct estimates
  - `fh` standard Fay-Herriot model
  - `fh_nested` alternative Fay-Herriot model fit on the full data with `var.fix = TRUE` and `nested = TRUE`
  - `cluster` cluster-level model
  - `cluster_strat` stratified cluster-level model

- **Admin2 FH safeguards**
  - For the standard `fh` admin2 fit, the wrapper can optionally identify and remove areas with near-zero direct variance before fitting
  - Dropped admin2 areas and clusters are recorded in `fit$dropped$fh_admin2`

- **Visualization**
  - Mean maps
  - CV maps
  - Ridge plots
---

## Installation

```r
# install.packages("remotes")
remotes::install_github("SZ-yang/MICSprev")
library(MICSprev)
```

---

## Core workflow

MICSprev follows a simple end-to-end workflow:

1. Build geographic inputs from the GPS shapefile
2. Read the raw MICS survey files
3. Process each indicator into a standardized format
4. Run direct and model-based estimators
5. Extract results
6. Generate maps and ridge plots

---

## Main exported functions

| Function | Description |
|---|---|
| `process_geo_mics()` | Process a MICS GPS shapefile and construct cluster/admin geographic objects |
| `process_NMR()` | Process birth-history data for neonatal mortality rate estimation |
| `process_ANC()` | Process women's survey data for ANC 4+ estimation |
| `process_DTP3()` | Process child survey data for DTP3 estimation |
| `run_indicator_models()` | Run direct, Fay-Herriot, nested Fay-Herriot, and cluster-level models |
| `plot_indicator_maps()` | Plot mean and CV maps for selected models |
| `plot_indicator_ridge()` | Plot ridge plots for posterior model comparisons |

---

## Example: full pipeline

```r
library(MICSprev)
library(haven)

# 1) Build geographic objects
geo <- process_geo_mics(
  country.name = "Nigeria",
  year = 2021,
  gps_path = "/path/to/NigeriaMICS2021GPS.shp"
)

# 2) Read survey files
bh <- read_sav("bh.sav")
wm <- read_sav("wm.sav")
ch <- read_sav("ch.sav")

# 3) Process indicators
dat_nmr  <- process_NMR(bh)
dat_anc  <- process_ANC(wm)
dat_dtp3 <- process_DTP3(ch)

# 4) Run models
fit_nmr <- run_indicator_models(
  data = dat_nmr,
  geo = geo,
  admin_levels = c(0, 1, 2),
  models = c("direct", "fh", "fh_nested", "cluster", "cluster_strat"),
  model = "bym2",
  ci = 0.95,
  drop_lowvar_admin2 = TRUE,
  lowvar_cutoff = 1e-30,
  verbose = TRUE
)

# 5) Access selected outputs
fit_nmr$fits[["0"]]$direct$res.natl
fit_nmr$fits[["1"]]$direct$res.admin1
fit_nmr$fits[["2"]]$fh$res.admin2
fit_nmr$fits[["2"]]$fh_nested$res.admin2
fit_nmr$dropped$fh_admin2

# 6) Plot maps
nmr_maps_adm1 <- plot_indicator_maps(
  fit = fit_nmr,
  geo = geo,
  admin = 1,
  models = c("direct", "fh", "fh_nested", "cluster", "cluster_strat"),
  scale = 1000,
  prefix = "NMR_NGA2021"
)

nmr_maps_adm2 <- plot_indicator_maps(
  fit = fit_nmr,
  geo = geo,
  admin = 2,
  models = c("direct", "fh", "fh_nested", "cluster", "cluster_strat"),
  scale = 1000,
  prefix = "NMR_NGA2021"
)

# 7) Ridge plot
p_nmr_ridge_adm1 <- plot_indicator_ridge(
  fit = fit_nmr,
  admin = 1,
  models = c("fh", "fh_nested", "cluster", "cluster_strat"),
  threshold = 12 / 1000,
  scale = 1000,
  prefix = "NMR_NGA2021"
)
```

---

## Output structure

`run_indicator_models()` returns a nested list indexed first by admin level and then by model name.

```r
fit$fits[["0"]]   # national level
fit$fits[["1"]]   # admin1
fit$fits[["2"]]   # admin2
```

Each admin level may contain some or all of:

- `direct`
- `fh`
- `fh_nested`
- `cluster`
- `cluster_strat`

For admin2 standard FH fits, dropped low-variance areas are stored in:

```r
fit$dropped$fh_admin2
```

with components such as:

- `bad_admin2`
- `bad_clusters`

---

## Notes on the two FH strategies

### `fh`
The standard Fay-Herriot route. At admin2, this can optionally drop areas with extremely small direct variance before model fitting.

### `fh_nested`
An alternative Fay-Herriot route fit on the full data using:

- `var.fix = TRUE`
- `nested = TRUE`

This is useful when users want to compare the default dropped-variance FH strategy against a nested/fixed-variance FH fit.

---

## Supported country-year combinations

Currently supported:

- Honduras 2019
- Nigeria 2021
- Ghana 2017–18
- Lao PDR 2017
- Gambia 2018
- Sierra Leone 2017
- Madagascar 2018
- Zimbabwe 2019
- Thailand 2022

---

## Example data and vignette

The package vignette demonstrates a complete Nigeria 2021 workflow using packaged example files stored in `inst/extdata/`. After installation, these can be located with `system.file()`.

```r
bh_path <- system.file("extdata", "bh.sav", package = "MICSprev")
wm_path <- system.file("extdata", "wm.sav", package = "MICSprev")
ch_path <- system.file("extdata", "ch.sav", package = "MICSprev")
gps_path <- system.file(
  "extdata",
  "NigeriaMICS2021GPS",
  "NigeriaMICS2021GPS.shp",
  package = "MICSprev"
)
```

To browse installed vignettes:

```r
browseVignettes("MICSprev")
```

---

## Dependencies

Core dependencies include:

- `surveyPrev`
- `SUMMER`
- `sf`
- `dplyr`
- `ggplot2`
- `haven`

---

## Development notes

For local development, a typical workflow is:

```r
devtools::load_all()
rmarkdown::render("vignettes/mics_workflow.Rmd")
devtools::document()
devtools::build_vignettes()
```

---

## License

Please add your preferred license information in `DESCRIPTION` and this README if not already finalized.
