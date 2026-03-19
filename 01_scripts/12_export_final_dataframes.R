# ----------------------------- Final Output -------------------------------- *
### In this final script you will export the data sets that will help you
### analize the remaining data.
###
### R Script to:
### - Load all the data frames created beforeç
### - Merge them into a single data
### - Export the final data

###-------------------------------------------------------------------------- *

# ---- Load libraries ---- *

library(dplyr)
library(terra)
library(sf)
###-------------------------------------------------------------------------- *

europe_extent <- ext(-10, 35, 35, 75)
ref_raster <- rast(ext(europe_extent),
                   resolution = 1 / 12,
                   crs = "EPSG:4326")

# ---- 01 Load previously prepared data ----
derived_in <- file.path("02_data", "03_derived_data")
spatial_in <- file.path("02_data", "02_spatial_data")
tree_data <- file.path("02_data", "01_tree_data")

# Metadata
con_broad <- read.csv(file.path(tree_data, "tree_species_con_broad.csv"))
meta <- read.csv(file.path(derived_in, "metadata_age.csv")) %>% 
  left_join(con_broad, by = "species_code")

coords <- meta %>%
  dplyr::select(longitude, latitude)

# Site level data
interseries_site <- read.csv(file.path(derived_in, "interseries_site.csv"))
chron_cv <- read.csv(file.path(derived_in, "chronology_cv.csv"))
clim_coefs <- readRDS(file.path(derived_in, "clim_coefs.rds"))
spei_coefs <- readRDS(file.path(derived_in, "spei_coefs.rds")) %>%
  rename(spei = varname,
         window_spei = window,
         spei_coef = coef,
         spei_significant = significant,
         spei_month = month,
         spe_ci_lower = ci_lower,
         spei_ci_upper = ci_upper)

itrdb_chronologies <- read.csv(file.path(derived_out, "rwi_calculations.csv"))

itrdb_resilience <- read.csv(file.path(derived_in, "itrdb_resilience.csv"))

# Tree level data
tree_age <- read.csv(file.path(derived_in, "tree_age.csv"))
bai <- read.csv(file.path(derived_in, "bai_calculations.csv"))
interseries_tree <- read.csv(file.path(derived_out, "interseries_tree.csv"))


# Other variables
clim_sites <- read.csv(file.path(derived_out, "clim_long.csv"))
spei_sites <- read.csv(file.path(derived_in, "spei_long.csv"))

# You can check which variables you from world clim you want to use here:
# https://www.worldclim.org/data/bioclim.html
annual_tmp <- rast(file.path(spatial_in, "01_climate", "01_worldclim",
                             "01_climate", "climate", "wc2.1_5m",
                             "wc2.1_5m_bio_1.tif"))
annual_pre <- rast(file.path(spatial_in, "01_climate", "01_worldclim",
                             "01_climate", "climate", "wc2.1_5m",
                             "wc2.1_5m_bio_12.tif"))

# Diversity
birds <- rast(file.path(spatial_in, "02_diversity", "01_animal_diver", "Birds",
                        "Richness_10km_Birds_v7_EckertIV_breeding_no_seabirds.tif"))
amphibians <- rast(file.path(spatial_in, "02_diversity", "01_animal_diver",
                             "Richness_10km_AMPHIBIANS_dec2017_EckertIV_Anura.tif"))
mammals <- rast(file.path(spatial_in, "02_diversity", "01_animal_diver", "Mammals",
                          "Richness_10km_MAMMALS_mar2018_EckertIV.tif"))

tree_rich <- rast(file.path(spatial_in, "02_diversity", "02_treediver",
                            "tree_sp_richness.tif"))

# Ecoregion
ecoregions_csv <- read.csv(file.path(spatial_in, "03_ecoregions",
                                     "KG_1986-2010_5m.csv")) %>%
  tibble() %>%
  mutate(climate = as.factor(climate),
         main_climate = as.factor(main.climate),
         mainc_code = as.numeric(mainc.code),
         KG = as.numeric(KG))

# Convert ecoregion to Spatial format
ecoregions_sf <- ecoregions_csv %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Rasterize the ecoregion data
# Fist for main climate (general one)
ecoregions_raster <- terra::rasterize(ecoregions_sf,
                                      rast(ext(ecoregions_sf),
                                           resolution = 1 / 10),
                                      field = "mainc_code")

ecoregions <- crop(ecoregions_raster,
                   europe_extent)

# Soil characteristics 
soil_grids <- rast(file.path(spatial_in, "04_soilgrids",
                             "soilgrids_europe.tif"))

# Topography
elevation <- rast(file.path(spatial_in, "05_topography",
                            "elevation_10KMmd_GMTEDmd.tif")) %>%
  crop(europe_extent)

slope <- rast(file.path(spatial_in, "05_topography",
                        "slope_10KMmd_GMTEDmd.tif")) %>%
  crop(europe_extent)

eastness <- rast(file.path(spatial_in, "05_topography",
                           "eastness_10KMmd_GMTEDmd.tif")) %>%
  crop(europe_extent)

northness <- rast(file.path(spatial_in, "05_topography",
                            "northness_10KMmd_GMTEDmd.tif")) %>%
  crop(europe_extent)

mountains <- rast(file.path(spatial_in, "06_mountains",
                            "k3classes.tif")) %>%
  crop(europe_extent)

# Bind it to metadata
metadata <- meta %>%
  bind_cols(terra::extract(ecoregions, coords) %>%
              dplyr::select(ecoregion = last)) %>%
  bind_cols(terra::extract(soil_grids$`soc_mean_0-5cm`, coords) %>%
              dplyr::select(soil_carbon = `soc_mean_0-5cm`)) %>%
  bind_cols(terra::extract(soil_grids$`phh2o_mean_0-5cm`, coords) %>%
              dplyr::select(soil_ph = `phh2o_mean_0-5cm`)) %>%
  bind_cols(terra::extract(soil_grids$`nitrogen_mean_0-5cm`, coords) %>%
              dplyr::select(soil_nitrogen = `nitrogen_mean_0-5cm`)) %>%
  bind_cols(terra::extract(elevation, coords) %>%
              dplyr::select(derived_elev = elevation_10KMmd_GMTEDmd)) %>%
  bind_cols(terra::extract(slope, coords) %>%
              dplyr::select(slope = slope_10KMmd_GMTEDmd)) %>%
  bind_cols(terra::extract(eastness, coords) %>%
              dplyr::select(eastness = eastness_10KMmd_GMTEDmd)) %>%
  bind_cols(terra::extract(northness, coords) %>%
              dplyr::select(northness = northness_10KMmd_GMTEDmd)) %>%
  bind_cols(terra::extract(mountains, coords) %>%
              dplyr::select(mountains = Rowid))
###-------------------------------------------------------------------------- *

# ---- 02 Prepare final Data frames ----

## Site level metadata
### You will have the interseries correlation, the chronology coefficent of 
### variation along with the metadata.
site_growth_data <- list(metadata, interseries_site, chron_cv) %>%
  reduce(left_join, by = c("rwl"))

### You will have the monthly climate correlations (temperature, max temperature,
### min, temperature along with the metadata.
site_clim_correlation <- metadata %>%
  left_join(clim_coefs, by = "rwl")

### You will have the SPEI correlations along with the metadata.
site_spei_correlation <- metadata %>%
  left_join(spei_coefs, by = "rwl")

### You will get the site chronologies along with the metadata.
site_year_chronologies <- itrdb_chronologies %>%
  left_join(metadata, by = "rwl")

### You will get the site resilience (Resistance, Recovery, and Resilience),
### along with the metadata
site_year_resilience <- itrdb_resilience %>%
  left_join(metadata, by = "rwl")

### You will get the site monthly information along with the metadata
site_month_climate <- clim_sites %>%
  left_join(spei_sites, by = c("rwl", "year", "month")) %>%
  left_join(metadata, by = "rwl")

### You will get the site growth trends along with the metadata. MAKE SURE THAT
### YOU UPDATE THE TIME WINDOW.
site_growth_trend <- itrdb_chronologies %>%
  group_by(rwl) %>%
  # ADJUST THE TIME WINDOW YOU WANT TO CALCULATE
  filter(year <= 2025,
         year >= 1975) %>%
  summarise(
    n_years = n(),
    growth_slope = coef(lm(std ~ year))[2],
    .groups = "drop") %>%
  left_join(metadata, by = "rwl")

## You will get the tree growth (BAI) level metadata along with the metadata.
tree_year_growth <- bai %>%
  left_join(tree_age, by = c("rwl", "tree")) %>% 
  left_join(interseries_tree, by = c("rwl", "tree"))

## You will get the tree growth trend along with the metadata. MAKE SURE THAT
## YOU UPDATE THE TIME WINDOW ACCORDINGLY
tree_growth_trend <- bai %>%
  group_by(rwl, tree) %>%
  filter(is.finite(bai)) %>% 
  # ADJUST THE TIME WINDOW YOU WANT TO CALCULATE
  filter(year <= 2025,
         year >= 1975) %>%
  summarise(
    growth_slope = coef(lm(bai ~ year))[2],
    .groups = "drop") %>%
  left_join(metadata, by = "rwl") %>% 
  left_join(interseries_tree, by = c("rwl", "tree"))


###-------------------------------------------------------------------------- *

# ---- 03 Export them ----
data_out <- file.path("03_output", "03_dataframes")

write.csv(x = site_growth_data,
          file = file.path(data_out, "site_growth_data.csv"),
          row.names = FALSE)
write.csv(x = site_clim_correlation,
          file = file.path(data_out, "site_clim_correlation.csv"),
          row.names = FALSE)
write.csv(x = site_spei_correlation,
          file = file.path(data_out, "site_spei_correlation.csv"),
          row.names = FALSE)
write.csv(x = site_year_chronologies,
          file = file.path(data_out, "site_year_chronologies.csv"),
          row.names = FALSE)
write.csv(x = site_year_resilience,
          file = file.path(data_out, "site_year_resilience.csv"),
          row.names = FALSE)
write.csv(x = site_month_climate,
          file = file.path(data_out, "site_month_climate.csv"),
          row.names = FALSE)
write.csv(x = site_growth_trend,
          file = file.path(data_out, "site_growth_trend.csv"),
          row.names = FALSE)
write.csv(x = tree_year_growth,
          file = file.path(data_out, "tree_year_growth.csv"),
          row.names = FALSE)
write.csv(x = tree_growth_trend,
          file = file.path(data_out, "tree_growth_trend.csv"),
          row.names = FALSE)
