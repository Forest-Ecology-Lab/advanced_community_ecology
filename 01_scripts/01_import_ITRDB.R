# -------------------------- Download the ITRDB ----------------------------- *
# Script for
# 1. ITRDB study metadata
# 2. ITRDB list of raw measurement files
# 3. Dataframe of rwl files
# 4. Output and save metadata and rwl files

library(dplR)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(readr)
library(admisc)
library(lubridate)
library(future)
library(furrr)
library(httr2)

plan(multisession, workers = parallel::detectCores() - 1)

# Ping WDS-Paleo API####
url <- request("https://www.ncei.noaa.gov/paleo-search/study/search.json") |>
  req_url_query(
    metadataOnly = "true",
    dataPublisher = "NOAA",
    dataTypeId = 18
  )

## This step can take some time depending on your internet connection and the
## status of the NOAA NCEI webpages
resp <- url |> req_perform()
stopifnot(resp_status(resp) == 200)

parsed <- resp_body_json(resp, simplifyVector = TRUE)

p_dat <- structure(list(
  url = url,
  content = parsed,
  response = resp
))

in_list <- unnest(p_dat$content$study, site)


# Filter for ITRDB sites #UPDATED to include just Europe and only > 1900 ###
# If you want to check other Parameters you can check the link:
# https://www.ncei.noaa.gov/access/paleo-search/study/params.json

itrdb_list <- in_list %>%
  filter(str_detect(studyName, "ITRDB"),
         str_detect(locationName, "Continent>Europe"),
         !str_detect(locationName, "Continent>Europe>Eastern Europe>Russia"),
         earliestYearCE <= 2026,
         mostRecentYearCE >= 1900,
         mostRecentYearCE <= 2026)

range(itrdb_list$earliestYearCE)

# 1. Create ITRDB site metadata table -------------------------------------####

# Helper function
get_itrdb_meta <- function(itrdb_list) {
  out_df <- itrdb_list %>%
    transmute(NOAAStudyId, siteName, studyCode, investigators,
              locationName,
              first_year = earliestYearCE, last_year = mostRecentYearCE)
  out_noaa_params <- itrdb_list %>%
    transmute(NOAAStudyId, doi, contr_year = year(contributionDate),
              NOAASiteId, url = onlineResourceLink)
  
  # Extract Lat-Long coordinates
  out_coords <- out_df %>%
    select(NOAAStudyId) %>%
    mutate(coords = pluck(itrdb_list, "geo", "geometry", "coordinates"),
           c_n = map_dbl(coords, length),
           latitude = map2_dbl(c_n, coords,
                               ~ {
                                 if (.x == 2) {
                                   as.numeric(.y[1])
                                 } else {
                                   mean(as.numeric(.y[1]), as.numeric(.y[2]))
                                 }
                               }
           ),
           longitude = map2_dbl(c_n, coords,
                                ~ {
                                  if (.x == 2) {
                                    as.numeric(.y[2])
                                  } else {
                                    mean(as.numeric(.y[3]), as.numeric(.y[4]))
                                  }
                                 }
           ),
           elevation = as.numeric(pluck(itrdb_list,
                                        "geo",
                                        "properties",
                                        "minElevationMeters"))
    ) %>%
    select(-coords, -c_n)
  # Extract species codes and names
  out_spp <- out_df %>%
    select(NOAAStudyId) %>%
    mutate(species_list = map(itrdb_list[["paleoData"]], "species"),
           #extract the species code
           species_code_df = map(species_list, ~.x[[1]]$speciesCode),
           #extract the species sci name
           species_name_df = map(species_list, ~.x[[1]]$scientificName),
           #combine codes if multiple
           species_code = map_chr(species_code_df, str_c, collapse = ", "),
           #combine names if multiple
           species_name = map_chr(species_name_df, str_c, collapse = ", ")
    ) %>%
    select(NOAAStudyId, species_code, species_name)
  out_age <- out_df %>%
    select(NOAAStudyId) %>%
    mutate(date_type = map_chr(itrdb_list[["paleoData"]], "timeUnit")
    )
  
  # Extract published references
  pub_list <- unnest(itrdb_list, publication)
  if (nrow(pub_list) > 0) {
    out_pub <- pub_list  %>%
      select(NOAAStudyId, citation) %>%
      group_by(NOAAStudyId) %>%
      summarize(reference = str_c(citation, collapse = "; "))
  } else {
    out_pub <- data.frame(
      NOAAStudyId = NA,
      reference = NA)
  }
  # Combine
  itrdb_site_meta <- out_df %>%
    mutate(country = map_chr(strsplit(locationName, ">", fixed = TRUE), ~ tail(.x, 1))) %>%
    left_join(out_age, by = "NOAAStudyId") %>%
    left_join(out_coords, by = "NOAAStudyId") %>%
    left_join(out_spp, by = "NOAAStudyId") %>%
    left_join(out_pub, by = "NOAAStudyId") %>%
    left_join(out_noaa_params, by = "NOAAStudyId") %>%
    filter(!is.na(studyCode))
  
  itrdb_site_meta
}

itrdb_site_meta <- get_itrdb_meta(itrdb_list)

# 2. ITRDB list of raw measurement files --------------------------------------

itrdb_rawmeas_files <- itrdb_list %>%
  #select the columns for Study ID and first year (earliestYearCE)
  select(NOAAStudyId, first_year = earliestYearCE) %>%
  #Gets columns "dataFile" and "paleoData"
  mutate(data = map(itrdb_list[["paleoData"]], "dataFile"),
         files = map(data, ~.x[[1]][c("linkText", "urlDescription", "fileUrl")]),
         vars = map(data, ~.x[[1]]$variables)
  ) %>%
  select(-data) %>%
  unnest(cols = c(files, vars)) %>%
  filter(urlDescription == "Raw Measurements",
         !str_detect(linkText, "-noaa")) %>%
  # found an error so added this filter to take out data that is not a df
  filter(map_lgl(vars, is.data.frame)) %>%
  unnest(cols = vars) %>%
  filter(cvWhat == "physical property>width>total ring width") %>%
  select(NOAAStudyId, first_year, linkText, fileUrl, cvWhat, cvUnit)


# 3. ITRDB RWL files ----------------------------------------------------------

# Helper function
read_itrdb_rwl <- function(first_year, fileUrl) {
  require(dplR)
  if (any(c(length(first_year) > 1,
            length(fileUrl) > 1))
  ) abort("the function can only handle one case at a time")
  # read file
  if (first_year < -999) {
    read.rwl(fileUrl, long = TRUE, format = "tucson")
  } else {
    read.rwl(fileUrl, format = "tucson")
  }
}


# clean up
rm(in_list, itrdb_list, resp, p_dat, parsed)

# Read rwls from the HTTPS server. This is SLOW
all_rwl <- itrdb_rawmeas_files %>%
  transmute(NOAAStudyId, linkText,
            RWL = future_map2(first_year, fileUrl,
                              ~ try(read_itrdb_rwl(.x, .y)),
                              .options = furrr_options(seed = TRUE),
                              .progress = TRUE))

# Sometimes sites are loaded on the ITRDB before the data are released,
# triggering access errors in the script above. This check_files will identify
# those so they can be removed from further use, as below

check_files <- all_rwl %>%
  mutate(check = map_chr(RWL, ~{
    if ("rwl" %in% class(.x)) {
      out <- "Good"
    } else out = "Error"
    out
  })
  ) %>%
  filter(check == "Error")

#will give the results of the previous code
cat("Files with errors:", nrow(check_files), "\n")

# 4. Output data --------------------------------------------------------------

# save rwl files in folder

out_dir <- file.path("02_data", "01_tree_data", "01_ITRDB_dendroecology", "rwl")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_dat <- all_rwl %>%
  # Filters out files that have problems from step before and check presicion
  filter(!NOAAStudyId %in% check_files$NOAAStudyId) %>%
  # Check precision from 1st series
  mutate(ndec = map_dbl(RWL, ~ numdec(.x[, 1]))) %>%
  inner_join(itrdb_rawmeas_files, by = c("NOAAStudyId", "linkText"))

pwalk(list(out_dat$RWL, out_dat$linkText, out_dat$ndec),
      ~ write.tucson(..1,
                     file.path(out_dir, ..2),
                     prec = if_else(..3 < 3, 0.01, 0.001)))

out_meta <- file.path("02_data", "01_tree_data", "01_ITRDB_dendroecology")
                     
# Save metadata files
write_csv(itrdb_rawmeas_files,
          file.path(out_meta,
                    "itrdb_raw_measurement_files.csv"))

write_csv(itrdb_site_meta,
          file.path(out_meta,
                    "itrdb_site_metadata.csv"))

future::plan(sequential)
# --------------------------------------------------------------------------- *
# --------------------------------------------------------------------------- *
# --------------------------------------------------------------------------- *
