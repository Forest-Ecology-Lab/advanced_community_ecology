# ----------------------------- Climate Correlations ------------------------ *
### 2025.03.20: R Script: Extracting climate variables used in the
### practice of Forest Ecology from the course Advance Community Ecology.
###
###
### Elaborated by Zabdi López and Rubén D. Manzanedo, Forest Ecology Lab,
### University of Bern

## R script to download:
## - 01 Prepare the data
## - 02 Calculate Climate Correlations
## - 03 Export maps

# ---- Load libraries ----
library(dplyr)
library(tibble)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

# ---- 01 Load Data ----

## High Resolution map of Europe
europe_map <- ne_countries(scale = "large",
                           continent = "Europe",
                           returnclass = "sf")

## Set Europe limits
europe_limits <- list(xlim = c(-10, 35), ylim = c(35, 75))

europe_crop <- st_crop(europe_map,
                       xmin = min(europe_limits$xlim),
                       xmax = max(europe_limits$xlim),
                       ymin = min(europe_limits$ylim),
                       ymax = max(europe_limits$ylim))

## Set europe extent
europe_extent <- ext(europe_limits$xlim[1], europe_limits$xlim[2],
                     europe_limits$ylim[1], europe_limits$ylim[2])

## ITRDB Metadata
itrdb_metadata <- read.csv(file.path("02_data", "01_tree_data",
                                     "01_ITRDB_dendroecology",
                                     "itrdb_site_metadata.csv")) %>%
  select(rwl = studyCode, first_year, last_year, country,
         latitude, longitude, elevation, species_code, species_name) %>%
  mutate(rwl = tolower(rwl)) %>% 
  filter(!is.na(latitude), !is.na(longitude))

## Load Climate Data
derived_in <- file.path("02_data", "03_derived_data")

itrdb_pre <- read.csv(file.path(derived_in, "precipitation_itrdb.csv"),
                      check.names = FALSE)
itrdb_tmp <- read.csv(file.path(derived_in, "temperature_itrdb.csv"),
                      check.names = FALSE)
itrdb_tmn <- read.csv(file.path(derived_in, "mintemp_itrdb.csv"),
                      check.names = FALSE)
itrdb_tmx <- read.csv(file.path(derived_in, "maxtemp_itrdb.csv"),
                      check.names = FALSE)
itrdb_spei <- read.csv(file.path(derived_in, "spei6month_itrdb.csv"),
                       check.names = FALSE)

# Tree Ring Derived Data
itrdb_treeage <- read.csv(file.path(derived_in, "tree_age.csv")) %>%
  tibble() %>%
  filter(first >= 1901)
itrdb_rwi <- read.csv(file.path(derived_in, "rwi_calculations.csv")) %>%
  tibble() %>%
  left_join(itrdb_metadata, by = "rwl") %>%
  filter(year >= 1901)
itrdb.bai <- read.csv(file.path(derived_in, "bai_calculations.csv")) %>%
  tibble() %>%
  left_join(itrdb_metadata, by = "rwl") %>%
  filter(year >= 1901)
itrdb_chr_cv <- read.csv(file.path(derived_in, "chronology_cv.csv")) %>%
  tibble() %>%
left_join(itrdb_metadata, by = "rwl")

itrdb_interseries_cor <- read.csv(file.path(derived_in, "mean_interseries.csv")) %>%
  tibble() %>%
  left_join(itrdb_metadata, by = "rwl")
  
# ---- 02 Climate Correlation Analysis ----

# Convert to long format
itrdb_pre_l <- itrdb_pre %>%
  pivot_longer(-rwl,
               names_to = "date",
               values_to = "pre") %>%
  separate(date, into = c("year","month"), sep = "_") %>%
  mutate(year = as.integer(year)) %>% 
  left_join(itrdb_metadata %>% select(rwl,elevation,species_code))


# Prepare rwi data
rwi_clim <- itrdb_pre_l %>%
  inner_join(itrdb_rwi %>% select(-first_year, -last_year, -elevation,-species_code, -species_name),
            by = c("rwl","year"))


