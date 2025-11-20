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
