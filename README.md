# MICSprev

**MICSprev** is an R package for processing **Multiple Indicator Cluster Survey (MICS)** data and fitting spatial small-area estimation models using the **surveyPrev** and **SUMMER** frameworks.

The package provides an end-to-end workflow for:

- processing MICS survey files into standardized indicator datasets  
- constructing geographic inputs from MICS GPS shapefiles  
- fitting direct and model-based estimators  
- comparing multiple model families  
- generating publication-ready maps and ridge plots  

---

## Features

### Indicators
- Neonatal mortality rate (**NMR**)  
- Antenatal care 4+ visits (**ANC**)  
- DTP3 vaccination coverage (**DTP3**)  

### Geographic preprocessing
- Builds `cluster.info`, `admin.info1`, `admin.info2`, `admin0`, `admin1`, `admin2`, and `geo`

### Model fitting
- `direct`
- `fh` (Fay–Herriot)
- `fh_nested` (fixed-variance nested FH)
- `cluster`
- `cluster_strat`

### Visualization
- Mean maps  
- CV maps  
- Ridge plots  
- Optional log-transformed maps for visual contrast  

---

## Installation

### Install dependencies

```r
install.packages(c(
  "sf", "dplyr", "ggplot2", "haven",
  "patchwork", "stringr", "lwgeom"
))
```

### Install SUMMER and surveyPrev

```r
install.packages("remotes")

remotes::install_github("richardli/SUMMER")
remotes::install_github("richardli/surveyPrev")
```

### Install INLA

```r
install.packages(
  "INLA",
  repos = c(
    getOption("repos"),
    INLA = "https://inla.r-inla-download.org/R/stable"
  ),
  dep = TRUE
)
```

### Install MICSprev

```r
install.packages("devtools")
devtools::install_github("SZ-yang/MICSprev")
```

---

## Data access

This package **does not include MICS data**.

Users must download:

- Survey data (`bh.sav`, `wm.sav`, `ch.sav`)  
- GPS shapefiles  

from the MICS website:

https://mics.unicef.org/

---

## Core workflow

1. Build geographic inputs  
2. Read survey data  
3. Process indicators  
4. Fit models  
5. Extract results  
6. Visualize outputs  

---

## Example

```r
library(MICSprev)
library(haven)

geo <- process_geo_mics(
  country.name = "Nigeria",
  year = 2021,
  gps_path = "/path/to/NigeriaMICS2021GPS.shp"
)

bh <- read_sav("bh.sav")
dat_nmr <- process_NMR(bh)

fit_nmr <- run_indicator_models(
  data = dat_nmr,
  geo = geo,
  admin_levels = c(1, 2),
  models = c("direct", "fh", "cluster", "cluster_strat")
)
```

---

## Visualization example

```r
# Standard scale
maps <- plot_indicator_maps(
  fit = fit_nmr,
  geo = geo,
  admin = 2,
  scale = 1000
)

# Log-transformed scale (for visual contrast)
maps_log <- plot_indicator_maps(
  fit = fit_nmr,
  geo = geo,
  admin = 2,
  scale = 1000,
  transform = "log1p"
)
```

The log transformation improves contrast when a few large values dominate the color scale.  
The underlying estimates are unchanged; only the plotted values are transformed.

---

## Output structure

```r
fit$fits[["0"]]  # national
fit$fits[["1"]]  # admin1
fit$fits[["2"]]  # admin2
```

---

## Notes

- Users must supply their own MICS data  
- Model fitting may take time depending on complexity  
- Log-transformed maps are for visualization only, not interpretation  

---

## Development

```r
devtools::load_all()
rmarkdown::render("vignettes/mics_workflow.Rmd")
devtools::document()
```

---

## License

Please specify your license in the `DESCRIPTION` file.
