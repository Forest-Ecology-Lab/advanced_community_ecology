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

# Summarise climate data
pre_gs <- itrdb_pre_l %>%
  filter(month %in% c("May","Jun","Jul","Aug")) %>%
  group_by(rwl, year) %>%
  summarise(pre_gs = mean(pre, na.rm = TRUE), .groups = "drop")

# Prepare rwi data
rwi_clim <- itrdb_rwi %>%
  left_join(pre_gs, by = c("rwl","year"))

# Correlation Analysis
cor_pre <- rwi_clim %>%
  group_by(rwl) %>%
  summarise(
    n_pairs = sum(complete.cases(rwi, pre_gs)),
    cor_pre_gs = ifelse(
      n_pairs >= 2,
      cor(rwi, pre_gs, use = "complete.obs"),
      NA_real_
    ),
    .groups = "drop"
  )

# Merge again with ITRDB Metadata
cor_pre_map <- cor_pre %>%
  left_join(itrdb_metadata, by = "rwl")

# 03 Plot the Correlations ----
plot1 <- ggplot() +
  geom_sf(data = europe_crop, fill = "#d0d9db", color = "gray30",
          linewidth = 0.2) +
  geom_point(data = cor_pre_map %>%
               #YOU CAN FILTER TO PLOT SPECIFIC GROUPS (species_code, elevation)
               filter(species_code == "PISY"), 
             aes(longitude, latitude, colour = cor_pre_gs)) +
  scale_colour_gradient2(
    low = "brown",
    mid = "white",
    high = "blue",
    midpoint = 0) +
  labs(title = "ITRDB Site correlation to Precipitation",
       color = "Country",
       x = "Longitude",
       y = "Latitude") +
  theme_classic() +
  theme(panel.grid = element_line(colour = "gray90", linewidth = 0.3))

# 04 Export maps ----
maps_out <- file.path("03_output", "02_maps")

ggsave(plot = plot1,
       path = file.path(maps_out),
       filename = "itrdb_pisy_climcorr.png",
       width = 20, height = 20, units = "cm",
       dpi = 300)
