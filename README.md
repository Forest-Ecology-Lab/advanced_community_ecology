# Advanced Community Ecology

Repository for **hands-on analysis in community ecology and dendrochronology**, developed for the *Advanced Community Ecology* course at the **University of Bern**.

This repository provides a **fully reproducible pipeline in R** to analyse tree-ring data, climate–growth relationships, and forest resilience.

------------------------------------------------------------------------

## 📌 Overview

This project guides users through a complete dendroecological workflow:

-   Processing tree-ring (`.rwl`) data
-   Building chronologies
-   Analysing climate–growth relationships
-   Detecting pointer years
-   Calculating forest resilience metrics

The workflow is modular and designed for **learning, reproducibility, and scalability**.

------------------------------------------------------------------------

## 📂 Repository Structure

```         
advanced_community_ecology/

├── 01_scripts/                # Main analysis pipeline
│   ├── 00_setup.R
│   ├── 01_import_ITRDB.R
│   ├── 02_download_expvariables.R
│   ├── 03_age_mean_growth.R
│   ├── 04_growth_trends.R
│   ├── 05_interseries_correlation.R
│   ├── 06_chronology_cv.R
│   ├── 07_climate_extractions.R
│   ├── 08_climate_correlations.R
│   ├── 09_spei_correlations.R
│   ├── 10_pointeryear_calculations.R
│   ├── 11_resilience_calculations.R
│   └── 12_export_final_dataframes.R
│
├── 02_data/
│   ├── 01_tree_data/          # ITRDB data (.rwl)
│   ├── 02_spatial_data/       # Climate, soil, topography
│   └── 03_derived_data/       # Outputs from scripts
│
├── outputs/                  # Figures and results
├── renv/                     # Reproducible R environment
├── renv.lock
└── README.md
```

------------------------------------------------------------------------

## ⚙️ Setup

Run this **first**:

``` r
source("01_scripts/00_setup.R")
```

This installs and loads all required packages using `renv`.

------------------------------------------------------------------------

## 🔄 Workflow (Step-by-step)

### 1. Data preparation

-   Import ITRDB data
-   Filter metadata (species, country, years)
-   Select `.rwl` files

------------------------------------------------------------------------

### 2. Growth metrics

-   Mean growth
-   Tree age
-   Growth variability

Script:

```         
03_age_mean_growth.R
```

------------------------------------------------------------------------

### 3. Growth trends

-   Detrending (Spline)
-   Chronology building
-   Basal Area Increment (BAI)

------------------------------------------------------------------------

### 4. Interseries correlation

-   Signal strength across trees
-   Site-level synchrony

------------------------------------------------------------------------

### 5. Chronology variability

-   Coefficient of variation (CV)

------------------------------------------------------------------------

### 6. Climate data extraction

-   CRU climate grids (precipitation, temperature)
-   SPEI drought index
-   Spatial extraction per site

------------------------------------------------------------------------

### 7. Climate–growth relationships

-   Correlation analysis using `treeclim::dcc()`
-   Monthly climate sensitivity

------------------------------------------------------------------------

### 8. SPEI correlations

-   Drought–growth relationships
-   Spatial and temporal patterns

------------------------------------------------------------------------

### 9. Pointer years

-   Detection of extreme growth years
-   Filtering by threshold and period

------------------------------------------------------------------------

### 10. Resilience components

Based on **Lloret et al. (2011)**:

-   Resistance
-   Recovery
-   Resilience
-   Relative resilience

------------------------------------------------------------------------

### 11. Final data integration

-   Merge all outputs
-   Add environmental variables (soil, climate, topography)
-   Generate analysis-ready datasets

------------------------------------------------------------------------

## 📊 Outputs

The pipeline generates:

-   Site-level metrics
-   Tree-level growth data
-   Climate correlations
-   Resilience indicators
-   Spatial datasets

All outputs are saved in:

```         
02_data/03_derived_data/
```

------------------------------------------------------------------------

## 🎓 Teaching Use

This repository is designed for:

-   Forest ecology courses
-   Community ecology training
-   Dendrochronology practicals

Each script corresponds to a **clear analytical step**, making it ideal for teaching.

------------------------------------------------------------------------

## 🧠 Learning Goals

By completing this workflow, users will:

-   Understand tree-ring data processing
-   Analyse climate–growth relationships
-   Interpret ecological stability metrics
-   Build reproducible R workflows

------------------------------------------------------------------------

## ⚠️ Notes

-   Large datasets (e.g., climate grids) are not included
-   Users must download external data (CRU, SPEI, WorldClim)
-   Computation time depends on dataset size

------------------------------------------------------------------------

## 🤝 Contributions

Feel free to:

-   Suggest improvements for scripts
-   Suggest new analyses
-   Report issues

------------------------------------------------------------------------

## 👤 Authors

-   Zabdi López
-   Rubén D. Manzanedo

------------------------------------------------------------------------

## ⭐ Acknowledgements

-   ITRDB (International Tree-Ring Data Bank)
-   OpenRing Project
-   Forest Ecology Lab (Bern)

------------------------------------------------------------------------
