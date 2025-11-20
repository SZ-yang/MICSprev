# MICSprev

MICSprev is a small R package that helps process MICS household survey data and run spatial analyses using the **surveyPrev** and **SUMMER** frameworks.  
The package provides:

- Functions to clean survey data:
  - Neonatal Mortality (NMR)
  - Antenatal Care (ANC 4+)
  - DTP3 vaccination
- Automatic downloading of administrative boundaries
- GPS cluster processing for any MICS country
- A full NMR analysis pipeline (direct estimation + FH + unit models + plots)

---

## Installation

```r
remotes::install_github("SZ-yang/MICSprev")
library(MICSprev)
```

## Main Functions

| Function          | Description                                         |
|------------------|-----------------------------------------------------|
| `process_NMR()`   | Preprocess birth history data for NMR estimation    |
| `process_ANC()`   | Preprocess women’s health data (ANC 4+)             |
| `process_DTP3()`  | Preprocess child vaccination data                   |
| `build_geo_mics()`| Download boundaries + process GPS cluster shapefiles|
| `run_nmr_analysis()` | Run full NMR analysis and generate plots/tables |
