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
library(purrr)
library(stringr)
library(tidyverse)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(treeclim)

# --------------------------------------------------------------------------- *
# --------------------------------------------------------------------------- *

# ---- 01 Load Data ----
# ----Map of Europe ---- *
## Load the map
europe_map <- ne_countries(scale = "large",
                           continent = "Europe",
                           returnclass = "sf")

## Set Europe limits
europe_limits <- list(xlim = c(-10, 35), ylim = c(35, 75))

# ---- Load ITRDB filtered Metadata ---- *
derived_in <- file.path("02_data", "03_derived_data")

# Remember to load the metadata previously created.
metadata <- read.csv(file.path(derived_in, "metadata_age.csv")) %>%
  tibble()

# ---- Load the chronologies ----- *
itrdb_rwi <- readRDS(file.path(derived_in, "chron_list.rds"))

# ---- Load the Climate Data ----- *
itrdb_pre <- read.csv(file.path(derived_in, "precipitation_itrdb.csv"),
                      check.names = FALSE)
itrdb_tmp <- read.csv(file.path(derived_in, "temperature_itrdb.csv"),
                      check.names = FALSE)
itrdb_tmn <- read.csv(file.path(derived_in, "mintemp_itrdb.csv"),
                      check.names = FALSE)
itrdb_tmx <- read.csv(file.path(derived_in, "maxtemp_itrdb.csv"),
                      check.names = FALSE)

## Convert to long format

### Precipitation
pre_l <- itrdb_pre %>%
  pivot_longer(-rwl,
               names_to = "date",
               values_to = "pre") %>%
  separate(date, into = c("year", "month"), sep = "_") %>%
  mutate(year = as.integer(year))

### Temperature
tmp_l <- itrdb_tmp %>%
  pivot_longer(-rwl,
               names_to = "date",
               values_to = "tmp") %>%
  separate(date, into = c("year", "month"), sep = "_") %>%
  mutate(year = as.integer(year))

### Minimum Temperature
tmn_l <- itrdb_tmn %>%
  pivot_longer(-rwl,
               names_to = "date",
               values_to = "tmn") %>%
  separate(date, into = c("year", "month"), sep = "_") %>%
  mutate(year = as.integer(year))

### Maximum Temperature
tmx_l <- itrdb_tmx %>%
  pivot_longer(-rwl,
               names_to = "date",
               values_to = "tmx") %>%
  separate(date, into = c("year", "month"), sep = "_") %>%
  mutate(year = as.integer(year))

### Join all together
clim_l <-  list(pre_l, tmp_l, tmn_l, tmx_l) %>%
  reduce(left_join, by = c("rwl", "year", "month"))
# Tree Ring Derived Data

# --------------------------------------------------------------------------- *
# --------------------------------------------------------------------------- *

# ---- 02 Loop of climate correlations through sites ----

# Set variables for loop
site_names <- names(itrdb_rwi)
clim_vars <- c("pre", "tmp", "tmn", "tmx")

results_list <- list()
k <- 1

# Loop through sites and variables
pb <- txtProgressBar(min = 0, max = length(site_names), style = 3)

for (i in seq_along(site_names)) {

  site <- site_names[i]

  chrono_site <- itrdb_rwi[[site]] %>%
    select(year, std) %>%
    as.data.frame() %>%
    column_to_rownames("year")

  for (v in clim_vars) {

    out <- tryCatch({

      climate_site <- clim_l %>%
        filter(rwl == site) %>%
        select(year, month, all_of(v)) %>%
        mutate(month = match(month, month.abb)) %>%
        as.data.frame()

      resp <- dcc(
                  chrono = chrono_site,
                  climate = climate_site,
                  selection = -6:9,
                  verbose = FALSE)


      data.frame(rwl = site,
                 clim_var = v,
                 window = rownames(resp$coef),
                 coef = resp$coef$coef,
                 significant = as.logical(resp$coef$significant),
                 month = resp$coef$month,
                 stringsAsFactors = FALSE)

    }, error = function(e) {

      data.frame(rwl = site,
        clim_var = v,
        window = NA,
        coef = NA_real_,
        stringsAsFactors = FALSE
      )
    })

    results_list[[k]] <- out
    k <- k + 1
  }

  setTxtProgressBar(pb, i)
}

close(pb)

# --------------------------------------------------------------------------- *
# --------------------------------------------------------------------------- *

# ---- 03 Plot the data ----

clim_coefs <- bind_rows(results_list) %>%
  mutate(
    month_plot = ifelse(
      grepl("prev", window),
      paste0("p", month),
      month
    ),
    month_plot = factor(
      month_plot,
      levels = c("pJun", "pJul", "pAug", "pSep", "pOct", "pNov", "pDec",
                 "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP")
    )
  )

saveRDS(object = clim_coefs,
        file = file.path(derived_in, "clim_coefs.rds"))

clim_plot <- clim_coefs %>%
  select(rwl, month_plot, clim_var, coef, significant) %>%
  left_join(metadata, by = "rwl")


# mean line and proportion of significant sites
clim_mean <- clim_plot %>%
  filter(species_code == "PISY") %>%
  group_by(clim_var, month_plot) %>%
  summarise(
    mean_coef = mean(coef, na.rm = TRUE),
    prop_sig = mean(significant, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(mark_sig = prop_sig >= 0.5)

ggplot() +
  geom_line(
    data = clim_plot,
    aes(x = month_plot, y = coef, group = rwl),
    colour = "grey70",
    alpha = 0.08,
    linewidth = 0.3
  ) +
  geom_line(
    data = clim_mean,
    aes(x = month_plot, y = mean_coef, group = 1),
    colour = "black",
    linewidth = 1
  ) +
  geom_point(
    data = clim_mean,
    aes(x = month_plot, y = mean_coef),
    colour = "black",
    size = 1.8
  ) +
  geom_point(
    data = clim_mean %>% filter(mark_sig),
    aes(x = month_plot, y = mean_coef),
    colour = "red",
    size = 3
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ clim_var, scales = "free_y") +
  labs(
    x = "Months",
    y = "Climate correlation",
    title = "Climate correlation signature across sites for PISY"
  ) +
  theme_minimal()
