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

# ---- Load libraries ----
library(dplyr)
library(dplR)

# 01 Prepare the data ----
itrdb_in <- file.path("02_data", "01_tree_data",
                      "01_ITRDB_dendroecology")

derived_out <- file.path("02_data", "03_derived_data")

metadata <- read.csv(file.path(itrdb_in, "itrdb_site_metadata.csv")) %>%
  tibble() %>%
  select(rwl = studyCode, first_year, last_year, country,
         latitude, longitude, elevation, species_code, species_name) %>%
  mutate(rwl = tolower(rwl))

# ----- FILTER ------- *
# Filter the data set depending on your research question by several variables.

# First filter by one or several of these:
# "rwl", "first_year", "last_year", "country", "latitude", "longitude",
# "elevation", "species_code", "species_name".
# Check the variables from the column you want to filter:
# unique(metadata$country); range(metadata$elevation, na.rm = TRUE)
# unique(metadata$species_name); range(metadata$latitude, na.rm = TRUE)

metadata_filter <- metadata %>%
  # In this case we will add a filter to only use "Pinus sylvestris L." (PISY)
  filter(species_code == "PISY") %>%
  distinct(rwl)

rwl_files <- list.files(file.path(itrdb_in, "rwl"),
                        pattern = "\\.rwl$",
                        full.names = TRUE)

# Here you apply the filter to the RWL files by the parameter you established
# before. NOTE: if you want to work with the whole data set remember that it
# will take more time to run.

rwl_files <- rwl_files[
  tools::file_path_sans_ext(basename(rwl_files)) %in% metadata_filter$rwl
]

# Check how many rwl files you end up with.
message(length(rwl_files), " .rwl files selected from metadata filter")

age_list <- vector("list", length(rwl_files))
pb <- txtProgressBar(min = 0, max = length(rwl_files), style = 3)

# 02 Loop to calculate growth variables ----
for (i in seq_along(rwl_files)) {

  rwl_path <- rwl_files[i]
  rwl_name <- tools::file_path_sans_ext(basename(rwl_path))

  if (i %% 20 == 0) setTxtProgressBar(pb, i)

  rwl <- tryCatch(
    {
      capture.output(
        rwl <- read.rwl(rwl_path),
        file = NULL
      )
      rwl
    },
    error = function(e) {
      message("Could not read ", rwl_name, ": ", e$message)
      NULL
    }
  )

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

  rwl_stat <- tryCatch(
    dplR::rwl.stats(rwl),
    error = function(e) {
      message("rwl.stats failed ", rwl_name, ": ", e$message)
      NULL
    }
  )

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

close(pb)

# Note: Some ITRDB .rwl files do not strictly follow Tucson format.
# dplR may issue warnings but usually rereads the file successfully.

# 03 Generate the data frames ----
tree_age <- bind_rows(age_list)

# 04 Export the data frames ----

# Tree level
write.csv(x = tree_age,
          file = file.path(derived_out, "tree_age.csv"),
          row.names = FALSE)

# Site Level
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
