# ---------------------------- download climatic data ----------------------- *
### 2025.03.20: R Script: Downloading explanatory variables used in the
### practice of Forest Ecolgoy from the course Advance Community Ecology.
###
###
###
### Elaborated by Zabdi López and Rubén D. Manzanedo, Forest Ecology Lab,
### University of Bern

## R script to download:
## - 01 WorldClim Data
## - 02 WWF Diversity Maps (Amphibians, Birds, Mammals)
## - 03 Climatic Ecoregions
## - 04 Soil Grid Maps (Carbon, pH, and Nitrogen)
## - 05 Topography Maps ()
## - 06 Mountain Classification

# ---- Load libraries ----
library(geodata)
library(sf)
library(terra)
library(soilDB)
library(ncdf4)

# ---- Define Europe limits ----*

europe_extent <- ext(-10, 35, 35, 75)
ref_raster <- rast(ext(europe_extent),
                   resolution = 1 / 12,
                   crs = "EPSG:4326")

# ---- 01 Climate data ----

## ---- WorldClim ---*
clim_dir <- file.path("02_data", "02_spatial_data", "01_climate",
                      "01_worldclim")
dir.create(clim_dir, recursive = TRUE, showWarnings = FALSE)

bioclim <- worldclim_global(var = "bio", res = 5,
                            path = file.path(clim_dir, "01_climate"))

## ---- Climatic Grid ----*
cgrid_out <- file.path("02_data", "02_spatial_data", "01_climate",
                       "02_climgrid")
dir.create(cgrid_out, recursive = TRUE, showWarnings = FALSE)
tmpcgrid <- tempfile(fileext = ".gz")

# Defining the variables precipitation ("pre"), temperature ("tmp"),
# minimun temperature ("tmn"), maximum temperature ("tmx").
vars <- c("pre", "tmp", "tmn", "tmx")

base_url <- paste0("https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.09/",
                   "cruts.2503051245.v4.09")

for (v in vars) {

  gz_file <- file.path(tempdir(),
                       paste0("cru_ts4.09.1901.2024.", v, ".dat.nc.gz"))
  out_file <- file.path(cgrid_out,
                        paste0("cru_ts4.09.1901.2024.", v, ".dat.nc"))

  url <- paste0(base_url, "/", v, "/", "cru_ts4.09.1901.2024.", v, ".dat.nc.gz")

  download.file(
    url = url,
    destfile = gz_file,
    mode = "wb",
    method = "libcurl"
  )

  R.utils::gunzip(
    filename = gz_file,
    destname = out_file,
    overwrite = TRUE,
    remove = TRUE
  )
}

## ---- Standardised Precipitation-Evapotranspiration Index ----*
spei_out <- file.path("02_data", "02_spatial_data", "01_climate", "03_spei")
dir.create(spei_out, recursive = TRUE, showWarnings = FALSE)

download.file("https://spei.csic.es/spei_database/nc/spei06.nc",
              destfile = file.path(spei_out, "spei06.nc"),
              mode = "wb", method = "libcurl")

# ---- 02 Diversity ----

## ---- Animal Diversity ----*

bio_dir <- file.path("02_data", "02_spatial_data", "02_diversity",
                     "01_animal_diver")
dir.create(bio_dir, recursive = TRUE, showWarnings = FALSE)
tmpbio <- tempfile(fileext = ".zip")

download.file(url = paste0("https://biodiversitymapping.org/wp-content/",
                           "uploads/2020/02/",
                           "BiodiversityMapping_TIFFs_2019_03d14_revised.zip"),
              destfile = tmpbio,
              mode = "wb", method = "libcurl")

unzip(tmpbio, exdir = bio_dir)

created_folder <- list.dirs(bio_dir, recursive = FALSE)[1]
file.rename(
            from = list.files(created_folder, full.names = TRUE),
            to   = file.path(bio_dir,
                             basename(list.files(created_folder,
                                                 full.names = TRUE))))

unlink(created_folder, recursive = TRUE)

## ---- Tree diversity ----*
tree_dir <- file.path("02_data", "02_spatial_data", "02_diversity",
                      "02_treediver")
dir.create(tree_dir, recursive = TRUE, showWarnings = FALSE)


download.file(url = "https://ndownloader.figshare.com/files/35133814",
              destfile = file.path(tree_dir, "tree_sp_richness.tif"),
              mode = "wb", method = "libcurl")

# ---- 03 Ecoregions ----

eco_dir <- file.path("02_data", "02_spatial_data", "03_ecoregions")
dir.create(eco_dir, recursive = TRUE, showWarnings = FALSE)
tmpeco <- tempfile(fileext = ".zip")

download.file(url = paste0("https://koeppen-geiger.vu-wien.ac.at/data/",
                           "Koeppen-Geiger-ASCII.zip"),
              destfile = tmpeco,
              mode = "wb", method = "libcurl")

unzip(tmpeco, exdir = eco_dir)

# ---- 04 Soil Grids ----
soil_dir <- file.path("02_data", "02_spatial_data", "04_soilgrids")
dir.create(soil_dir, recursive = TRUE, showWarnings = FALSE)

bbox_sf <- st_as_sfc(st_bbox(c(
  xmin = xmin(europe_extent),
  xmax = xmax(europe_extent),
  ymin = ymin(europe_extent),
  ymax = ymax(europe_extent)
), crs = crs(ref_raster)))

soilgrids_europe <- fetchSoilGrids(
  x = bbox_sf,
  variables = c("soc", "phh2o", "nitrogen"),
  depth_intervals = "0-5",
  summary_type = "mean",
  target_resolution = c(10000, 10000),
  grid = TRUE,
  progress = TRUE
)

soilgrids_europe_proj <- terra::project(
  soilgrids_europe,
  ref_raster,
  method = "bilinear"
)

writeRaster(soilgrids_europe_proj,
            file.path(soil_dir, "soilgrids_europe.tif"),
            overwrite = TRUE)

# ---- 05 Topography ----
topo_out <- file.path("02_data", "02_spatial_data", "05_topography")
dir.create(topo_out, recursive = TRUE, showWarnings = FALSE)

vars <- c("elevation", "slope", "eastness", "northness")

base_url <- "https://data.earthenv.org/topography"

for (v in vars) {

  file <- paste0(v, "_10KMmd_GMTEDmd.tif")
  url  <- paste0(base_url, "/", file)

  download.file(url,
                destfile = file.path(topo_out, file),
                mode = "wb", method = "libcurl")
}

# ---- 06 Mountains Classification ----
mount_out <- file.path("02_data", "02_spatial_data", "06_mountains")
dir.create(mount_out, recursive = TRUE, showWarnings = FALSE)
tmpmount <- tempfile(fileext = ".zip")

download.file(paste0("https://www.sciencebase.gov/catalog/file/get",
                     "/638fbf72d34ed907bf7d3080?f=__disk__93%2Fcc%",
                     "2Fb1%2F93ccb1d82335b2b475f4f728bcca635ba355d4eb"),
              destfile = tmpmount,
              mode = "wb", method = "libcurl")

unzip(tmpmount, exdir = mount_out)

tree_out <- file.path("02_data", "01_tree_data")
dir.create(tree_out, recursive = TRUE, showWarnings = FALSE)

download.file(url = paste0("https://www.ncei.noaa.gov/pub/data/paleo/",
                           "templates/tree-species-code.csv"),
              destfile = file.path(tree_out, "tree-species-code.csv"),
              mode = "wb", method = "libcurl")
