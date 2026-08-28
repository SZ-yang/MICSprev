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

**Country-agnostic** — run on any supported MICS survey:

- Neonatal mortality rate (**NMR**)  
- Antenatal care 4+ visits (**ANC**)  
- DTP3 vaccination coverage (**DTP3**)  

#### Nigeria-only indicators

The following are translated from the official MICS6 SPSS tabulation syntax
against the **Nigeria MICS6 (2021)** recode variable names, and are gated to
Nigeria (see [Country gating](#country-gating)). The MICS6 source table for
each is given in brackets.

| Sector | Indicator | Code | Recode |
|---|---|---|---|
| Health | ANC, 1+ visit from a skilled provider [TM.4.1] | `ANC1` | wm |
| Health | Skilled attendance at birth [TM.6.2] | `SBA` | wm |
| Health | Post-natal check, newborn, <2 days [TM.8.2] | `PNCNB` | wm |
| Health | Post-natal check, mother, <2 days [TM.8.7] | `PNCMOM` | wm |
| Health | Penta-1 coverage, 12–23 months [TC.1.1] | `PENTA1` | ch |
| Child protection | Birth registration, under 5 [PR.1.1] | `BIRTHREG` | ch |
| Child protection | Women 20–49 married before 18 [PR.4.1W] | `CHILDMARRIAGE` | wm |
| Child protection | Child labour, 5–17 [PR.3.3] | `CHILDLABOUR` | fs |
| Child protection | Child labour incl. hazardous conditions | `CHILDLABOURHAZ` | fs |
| Child protection | FGM among daughters 0–14 [PR.5.3] | `FGMDAUGHTER` | fg + bh + wm |
| Education | Out-of-school, primary age [LN.2.3] | `OOSPRIMARY` | hl |
| Education | Out-of-school, lower secondary age [LN.2.4] | `OOSSECONDARY` | hl |
| Education | Foundational reading skills, 7–14 [LN.4.1] | `READING` | fs |
| Education | Foundational numeracy skills, 7–14 [LN.4.2] | `MATH` | fs |
| Education | Net intake rate, grade 1 [LN.2.2] | `NETINTAKE` | hl |
| Education | Primary completion rate [LN.2.7] | `PRIMARYCOMPL` | hl |
| Education | Lower secondary completion rate [LN.2.7] | `LOWSECCOMPL` | hl |
| Education | Upper secondary completion rate [LN.2.7] | `UPSECCOMPL` | hl |
| Nutrition | Exclusive breastfeeding, <6 months [TC.7.3] | `EBF` | ch |
| Nutrition | Vitamin A supplementation, 6–59 months | `VITA` | ch |
| Nutrition | Age-appropriate breastfeeding, 0–23 m [TC.7.5] | `APPROPBF` | ch |
| Nutrition | Minimum acceptable diet, 6–23 m [TC.7.7] | `MINACCEPTDIET` | ch |
| Nutrition | Severe food insecurity, 12 months (FIES) | `FOODINSEC` | hh |
| WASH | At least basic drinking water [WS.1.2] | `BASICWATER` | hh |
| WASH | At least basic sanitation [WS.3.2] | `BASICSAN` | hh |
| WASH | Open defecation [WS.3.1] | `OPENDEF` | hh |
| WASH | Sufficient drinking water [WS.1.5] | `WATERSUFF` | hh |
| WASH | Menstrual hygiene, women 15–49 [WS.4.1] | `MENSTRUAL` | wm |
| Social policy | Handwashing facility, water and soap [WS.2.1] | `HANDWASH` | hh |
| Social policy | Functional difficulty, 5–17 [EQ.1.4] | `FUNCDIFF517` | fs |
| Social policy | Functional difficulty, 2–17 [EQ.1.4] | `FUNCDIFF217` | fs + ch |
| Social policy | Health insurance, under 5 [EQ.2.3] | `HEALTHINS` | ch |
| Social policy | Social transfers, children <18 [EQ.2.7] | `SOCTRANSFER` | hh + hl |
| Social policy | Women 15–49 with a bank account | `BANKACCT` | wm |
| Social policy | Women 15–49 who borrowed money | `BORROWED` | wm |

`indicator_countries()` returns the authoritative registry at runtime.

##### Not available in Nigeria MICS6 (2021)

Three requested indicators **cannot be derived from this survey** because the
modules were not fielded or not released:

| Indicator | Reason |
|---|---|
| Stunting (HAZ < −2) | No anthropometry in any of the eight recode files — no height, weight, oedema, or z-score variables, and `ch.sav` has no `AN` module |
| Severe wasting (WHZ < −3) | As above |
| Micronutrient powders, 6–59 months | No MNP item anywhere in `ch.sav` |

These would need another source, such as the NDHS or the MICS 2016–17 round.

##### Indicators needing confirmation before publication

Several indicators required a judgement the source syntax does not settle.
Each is marked with a `# FLAG:` comment at the relevant line and exposed as a
function argument where possible:

- **`VITA`** measures *ever received*, not *received in the last 6 months* —
  Nigeria asks only `IM27B` "Child ever given Vitamin A".
- **`BORROWED`** — `FN5` carries no time qualifier, so the "last 12 months"
  reference period rests on questionnaire wording, not on the recode.
- **`FOODINSEC`** is a raw FIES score, not Rasch-equated; it is roughly double
  typical published SDG 2.1.2 estimates for Nigeria and should not be compared
  with FAO figures without equating first.
- **`READING`** — story lengths and the story 2 literal/inferential question
  split are country-specific; see `?process_READING`.
- **Completion rates** — Nigeria's education level codes differ entirely from
  the MICS6 template, and the grouping of VEI/IEI and Secondary Technical is a
  judgement call; see `?process_PRIMARYCOMPL`.

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

## Country gating

Country-specific indicators are gated so a Nigeria-only definition cannot be
run silently against another country's recode — MICS6 variable names are
shared across countries, so a column check alone would not catch it.

The gate has three layers:

1. **`process_indicator()` requires `country`** for any indicator that is not
   registered as `"all"`. `NMR`, `ANC` and `DTP3` stay country-agnostic, so
   existing calls that pass no `country` are unaffected.
2. **`run_indicator_models()` re-checks** the indicator against `geo$country`,
   which every geo handler now stamps. This catches a processed dataset that
   is saved and later paired with the wrong `geo`. A `geo` object built before
   this existed has no `country` and produces a warning, not an error.
3. **Required-column checks** in each `process_*()` function name every
   missing variable and state that the indicator is Nigeria-only.

```r
# Country-agnostic — unchanged
process_indicator(bh, "NMR")

# Nigeria-only — country is required
process_indicator(wm, "ANC1", country = "Nigeria")

# Errors: indicator `ANC1` is not supported for country `ghana`
process_indicator(wm, "ANC1", country = "Ghana")

# The registry, at runtime
indicator_countries()
```

---

## Data access

This package **does not include MICS data**.

Users must download:

- Survey data (`hh.sav`, `hl.sav`, `wm.sav`, `bh.sav`, `fg.sav`, `ch.sav`,
  `fs.sav`, `mn.sav` — which files you need depends on the indicator)  
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

### Nigeria-only indicators

Pass `country`, and read whichever recode the indicator needs:

```r
wm <- read_sav("wm.sav")
hh <- read_sav("hh.sav")

# Skilled attendance at birth
dat_sba <- process_indicator(wm, "SBA", country = "Nigeria")

fit_sba <- run_indicator_models(
  data = dat_sba,
  geo = geo,                       # geo$country == "nigeria" is re-checked
  admin_levels = c(1, 2),
  models = c("direct", "fh")
)

# Open defecation. WASH indicators are population-based: `weight` is the
# household weight multiplied by household size, as the MICS6 tables do.
dat_od <- process_indicator(hh, "OPENDEF", country = "Nigeria")
```

#### Indicators spanning several recode files

A few indicators genuinely need more than one file. Pass the first as `data`
and the rest by name — `process_indicator()` forwards them:

```r
fg <- read_sav("fg.sav"); ch <- read_sav("ch.sav"); hl <- read_sav("hl.sav")

# fg.sav rosters only the daughters of women who had heard of FGM; the rest
# come from the birth history. Using fg.sav alone overstates prevalence.
dat_fgm <- process_indicator(fg, "FGMDAUGHTER", country = "Nigeria",
                             bh = bh, wm = wm)

# Ages 5-17 from fs.sav stacked with ages 2-4 from ch.sav
dat_fd <- process_indicator(fs, "FUNCDIFF217", country = "Nigeria", ch = ch)

# Household programme receipt, with children under 18 as the denominator
dat_st <- process_indicator(hh, "SOCTRANSFER", country = "Nigeria", hl = hl)
```

#### Country-specific parameters

Indicators whose definition depends on national policy expose those choices as
arguments rather than burying them. Defaults encode Nigeria:

```r
# Nigeria 6-3-3: primary ages 6-11, junior secondary 12-14
process_indicator(hl, "OOSSECONDARY", country = "Nigeria",
                  primary_entry_age = 6, primary_grades = 6,
                  lower_sec_grades = 3)

# Reading passage lengths differ by language; defaults are the modal lengths
process_indicator(fs, "READING", country = "Nigeria",
                  story1_words = 72, story2_words = 63)

# FIES: 7 of 8 items is severe, 4 of 8 is moderate-or-severe
process_indicator(hh, "FOODINSEC", country = "Nigeria", threshold = 4)
```

---

## Visualization example

```r
# NMR: darker color = higher mortality (direction = -1)
maps_nmr <- plot_indicator_maps(
  fit = fit_nmr,
  geo = geo,
  admin = 1,
  scale = 1000,
  direction = -1
)

# ANC / DTP3: darker color = lower coverage (direction = 1)
maps_anc <- plot_indicator_maps(
  fit = fit_anc,
  geo = geo,
  admin = 1,
  scale = 100,
  direction = 1
)

# Log-transformed mean map (for visual contrast when a few areas dominate)
maps_nmr_log <- plot_indicator_maps(
  fit = fit_nmr,
  geo = geo,
  admin = 2,
  scale = 1000,
  direction = -1,
  transform = "log1p"
)
```

### Color direction

The `direction` parameter controls the color scale of the **mean map only**:

| Indicator | `direction` | Interpretation |
|-----------|-------------|----------------|
| NMR | `-1` | darker = higher mortality |
| ANC | `1` | darker = lower coverage |
| DTP3 | `1` | darker = lower coverage |

The **CV map** always uses `direction = -1` (darker = higher uncertainty) regardless of indicator.

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
