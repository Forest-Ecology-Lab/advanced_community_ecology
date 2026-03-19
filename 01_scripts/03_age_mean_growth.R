# -------------------------------- Mean Growth ------------------------------ *
### 2025.03.20: R Script: Calculate Mean growth from the ITRDB downloaded data.
### Practice of forest ecology from the course Advance Community Ecology.
###
###
###
### Elaborated by Zabdi López and Rubén D. Manzanedo, Forest Ecology Lab,
### University of Bern

## R script to:
## - 01 Loop to calculate growth variables
## - 02 Generate dataframes
## - 03 Export dataframes
# --------------------------------------------------------------------------- *
# --------------------------------------------------------------------------- *

# ---- Load libraries ----
library(dplyr)
library(dplR)

# 01 Prepare the data ----
# Set the folders path to locate the files.
itrdb_in <- file.path("02_data", "01_tree_data",
                      "01_ITRDB_dendroecology")

derived_in <- file.path("02_data", "03_derived_data")

# Load the metadata
metadata <- read.csv(file.path(itrdb_in, "itrdb_site_metadata.csv")) %>%
  tibble() %>%
  select(rwl = studyCode, first_year, last_year, country,
         latitude, longitude, elevation, species_code, species_name) %>%
  mutate(rwl = tolower(rwl))

# ----- FILTER THE METADATA ------- *
# DEPENDING ON YOUR RESEARCH QUESTION YOU WILL HAVE TO APPLY THIS FILTER!!!
# Filter the data set depending on which variables you need.

# You can filter by one or several of these:
# "rwl", "first_year", "last_year", "country", "latitude", "longitude",
# "elevation", "species_code", "species_name".
# Check the variables from the column you want to filter:
# unique(metadata$country); range(metadata$elevation, na.rm = TRUE)
# unique(metadata$species_name); range(metadata$latitude, na.rm = TRUE)

# Here are some examples you can use to filter your data:
# NOTE: if you want to work with the whole data set remember that it
# will take more time to run.

# == BY SPECIES ==
## Check the unique values
table(metadata$species_code)
## Apply the filter
metadata_filter <- metadata %>%
  filter(species_code == "PISY") %>%
  distinct(rwl)

# == BY COUNTRY ==
## Check the unique values
unique(metadata$country)
## Apply the filter
metadata_filter <- metadata %>%
  filter(country == "Norway") %>%
  distinct(rwl)

# == BY YEAR ==
## Check the unique values
range(metadata$last_year, na.rm = TRUE)
## Apply the filter OLDER THAN
metadata_filter <- metadata %>%
  filter(last_year >= 1980) %>%
  distinct(rwl)

## Apply the filter YOUNGER THAN
metadata_filter <- metadata %>%
  filter(last_year <= 1980) %>%
  distinct(rwl)

## Apply the filter YEAR RANGE
metadata_filter <- metadata %>%
  filter(last_year >= 1980,
         last_year <= 2000) %>%
  distinct(rwl)

# == OR A MIX ==
## Apply the filter for country and years
metadata_filter <- metadata %>%
  filter(country == "Norway",
         last_year >= 1980,
         last_year <= 2000) %>%
  distinct(rwl)

# Check that your filter worked
cat("metadata files:", nrow(metadata),
    "\nfiltered metadata files:", nrow(metadata_filter), "\n")
print(metadata_filter)

# Export the filter to use it in following scripts
saveRDS(object = metadata_filter,
        file = file.path(derived_in, "metadata_filter.rds"))

# This will set the path for each rwl file
rwl_files <- list.files(file.path(itrdb_in, "rwl"),
                        pattern = "\\.rwl$",
                        full.names = TRUE)

# Here you apply the filter to the RWL files by the parameter you established
# before.

# NOTE! If you do not want to run a filter and want to use the whole data set,
# please run metadata_filter <- metadata so you will not have to change all the
# script.

rwl_files <- rwl_files[
  tools::file_path_sans_ext(basename(rwl_files)) %in% metadata_filter$rwl
]

# Check how many rwl files you end up with.
message(length(rwl_files), " .rwl files selected from metadata filter")

# --------------------------------------------------------------------------- *
# --------------------------------------------------------------------------- *

# 02 Loop to calculate growth variables ----

# This loop will go through all the RWL files it will read the RWL files, then
# run the function rwl.stats() from the packge dplR to calculate the growth 
# stats of each tree and site.

# Set a list that will contain all the results from each iteration.
age_list <- vector("list", length(rwl_files))
# This is a proogress bar that will help visualize
pb <- txtProgressBar(min = 0, max = length(rwl_files), style = 3)

# Loop starts
# i = interation
# Each iteration is one RWL so it will run the entire loop fully once per RWL
# file. So the time it takes will depend on how long is your filter.
for (i in seq_along(rwl_files)) {

  # Set the path of the file and file names
  rwl_path <- rwl_files[i]
  rwl_name <- tools::file_path_sans_ext(basename(rwl_path))
  
  # Update the progress bar
  if (i %% 20 == 0) setTxtProgressBar(pb, i)

  # reads the rwl file without printing the results
  rwl <- tryCatch(
    {
      capture.output(
        rwl <- read.rwl(rwl_path),
        file = NULL
      )
      rwl
    },
    # If there is an error reading the RWL this will give you a message
    error = function(e) {
      message("Could not read ", rwl_name, ": ", e$message)
      NULL
    }
  )

  # It there is an error the output will be NAs and have a Flag
  if (is.null(rwl)) {
    age_list[[i]] <- tibble(
      rwl = rwl_name,
      tree = NA_character_,
      first = NA_real_,
      last = NA_real_,
      age = NA_real_,
      mean_growth = NA_real_,
      growth_sd = NA_real_,
      ar1 = NA_real_,
      rwl_meanage = NA_real_,
      rwl_totage = NA_real_,
      flag = "reading rwl file failed"
    )
    next
  }

  # Once the RWL is read it will calculate the rwl.stats which will generate
  # the information of growth
  rwl_stat <- tryCatch(
    dplR::rwl.stats(rwl),
    error = function(e) {
      message("rwl.stats failed ", rwl_name, ": ", e$message)
      NULL
    }
  )

  # If not, then it will give you NAs and a Flag
  if (is.null(rwl_stat)) {
    age_list[[i]] <- tibble(
      rwl = rwl_name,
      tree = NA_character_,
      first = NA_real_,
      last = NA_real_,
      age = NA_real_,
      mean_growth = NA_real_,
      growth_sd = NA_real_,
      ar1 = NA_real_,
      rwl_meanage = NA_real_,
      rwl_totage = NA_real_,
      flag = "rwl.stats failed"
    )
    next
  }

  # Will save the output as a tibble and store it in a list that will be used
  # later
  age_list[[i]] <- tibble(rwl_stat) %>%
    mutate(rwl = rwl_name) %>%
    select(
      rwl, tree = series, first, last, age = year,
      mean_growth = mean, growth_sd = stdev, ar1
    ) %>%
    mutate(
      rwl_meanage = mean(age, na.rm = TRUE),
      rwl_totage = max(last, na.rm = TRUE) - min(first, na.rm = TRUE) + 1,
      flag = NA_character_
    )
}

# Closes the Progress bar
close(pb)

# Note: Some ITRDB .rwl files do not strictly follow Tucson format.
# dplR may issue warnings but usually rereads the file successfully.

# 03 Generate the data frames ----
tree_age <- bind_rows(age_list)

# 04 Export the data frames ----
derived_out <- file.path("02_data", "03_derived_data")

# Data Frame at the tree level
write.csv(x = tree_age,
          file = file.path(derived_out, "tree_age.csv"),
          row.names = FALSE)

# Data Frame at the Site Level
metadata_age <- tree_age %>%
  group_by(rwl) %>%
  summarise(rwl_meanage = first(rwl_meanage),
            rwl_totage = first(rwl_totage),
            mean_AR1 = mean(ar1, na.rm = TRUE),
            max_age = max(age, na.rm = TRUE),
            min_age = min(age, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(metadata, by = "rwl")

write.csv(x = metadata_age,
          file =  file.path(derived_out, "metadata_age.csv"),
          row.names = FALSE)
# --------------------------------------------------------------------------- *
# --------------------------------------------------------------------------- *
# --------------------------------------------------------------------------- *
