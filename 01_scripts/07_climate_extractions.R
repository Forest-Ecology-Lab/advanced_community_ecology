# ----------------------------- Climate Extractions ------------------------- *
### 2025.03.20: R Script: Extracting climate variables used in the
### practice of Forest Ecology from the course Advance Community Ecology.
###
###
### Elaborated by Zabdi López and Rubén D. Manzanedo, Forest Ecology Lab,
### University of Bern

## R script to download:
## - 01 Prepare the data
## - 02 Extract data from Climate Grids
## - 03 Export data frames

# ---- Load libraries ----
library(terra)
library(dplyr)

# --------------------------------------------------------------------------- *

# ---- 01: Prepare the data ----

derived_in <- file.path("02_data", "03_derived_data")

# Remember to load the metadata previously created.
metadata <- read.csv(file.path(derived_in, "metadata_age.csv")) %>%
  tibble()

# load and stack the climate variables with monthly resolution,
# pre = precipitation, tmp = mean temperatue tmn = minimum temperature,
# tmx = maximum temperature, spei= 6month spei

climgrid_in <- file.path("02_data", "02_spatial_data", "01_climate",
                         "02_climgrid")
spei_in <- file.path("02_data", "02_spatial_data", "01_climate", "03_spei")

europe_extent <- ext(-10, 35, 35, 75)

pre <-  rast(file.path(climgrid_in, "cru_ts4.09.1901.2024.pre.dat.nc"),
             subds = "pre") %>%
  crop(europe_extent)
tmp <- rast(file.path(climgrid_in, "cru_ts4.09.1901.2024.tmp.dat.nc"),
            subds = "tmp") %>%
  crop(europe_extent)
tmn <- rast(file.path(climgrid_in, "cru_ts4.09.1901.2024.tmn.dat.nc"),
            subds = "tmn") %>%
  crop(europe_extent)
tmx <- rast(file.path(climgrid_in, "cru_ts4.09.1901.2024.tmx.dat.nc"),
            subds = "tmx") %>%
  crop(europe_extent)
spei <- rast(file.path(spei_in, "spei06.nc")) %>%
  crop(europe_extent)

coordinates <- metadata %>%
  select(rwl,
         lon = longitude,
         lat = latitude) %>%
  distinct(lon, lat, .keep_all = TRUE)

pts <- vect(coordinates, geom = c("lon", "lat"), crs = "EPSG:4326")

years  <- 1901:2024
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
            "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

clim_names <- paste(rep(years, each = 12),
                    rep(months, times = length(years)),
                    sep = "_")

# --------------------------------------------------------------------------- *

# ---- 02: Extract data from CLimatic Grids ----
pre_sites <- as.data.frame(terra::extract(pre, pts)) %>%
  mutate(rwl = coordinates$rwl, .before = 1) %>%
  select(-ID) %>%
  setNames(c("rwl", clim_names))
tmp_sites <- as.data.frame(terra::extract(tmp, pts)) %>%
  mutate(rwl = coordinates$rwl, .before = 1) %>%
  select(-ID) %>%
  setNames(c("rwl", clim_names))
tmn_sites <- as.data.frame(terra::extract(tmn, pts)) %>%
  mutate(rwl = coordinates$rwl, .before = 1) %>%
  select(-ID) %>%
  setNames(c("rwl", clim_names))
tmx_sites <- as.data.frame(terra::extract(tmx, pts)) %>%
  mutate(rwl = coordinates$rwl, .before = 1) %>%
  select(-ID) %>%
  setNames(c("rwl", clim_names))
spei_sites <- as.data.frame(terra::extract(spei, pts)) %>%
  mutate(rwl = coordinates$rwl, .before = 1) %>%
  select(-ID) %>%
  select(-c(spei_1:spei_12)) %>%
  setNames(c("rwl", clim_names[-(1:12)]))

# --------------------------------------------------------------------------- *

# ---- 03: Export data sets ----
write.csv(x = pre_sites,
          file = file.path(derived_in, "precipitation_itrdb.csv"),
          row.names = FALSE)
write.csv(x = tmp_sites, file = file.path(derived_in, "temperature_itrdb.csv"),
          row.names = FALSE)
write.csv(x = tmn_sites, file = file.path(derived_in, "mintemp_itrdb.csv"),
          row.names = FALSE)
write.csv(x = tmx_sites, file = file.path(derived_in, "maxtemp_itrdb.csv"),
          row.names = FALSE)
write.csv(x = spei_sites, file = file.path(derived_in, "spei6month_itrdb.csv"),
          row.names = FALSE)
# --------------------------------------------------------------------------- *
# --------------------------------------------------------------------------- *
# --------------------------------------------------------------------------- *
