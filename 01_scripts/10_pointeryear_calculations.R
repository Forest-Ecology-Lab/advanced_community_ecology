###-------------------------Resilience components---------------------------- *

### 2025.02.14: R script to calculate the resilience components proposed by *
### Lloret et al. (2011) for the course of Advanced Community Ecology.
###
###
###
### Elaborated by Zabdi López and Rubén D. Manzanedo, Forest Ecology Lab,
### University of Bern

### R Script to:
### -Load files and verify it has more that 3 trees and more than 26 years.
### -Detrend the rwl file
### -calculate pointer years (using a normalization in a moving window)
### -calculate missing rings identified as 0.01 and 0.001
### -extract the information for the pointer years
### -Flags:
### -- 1. Data set with less than 3 series
### -- 2. Data set longer than 26 years
### -- 3. Detrended data set longer than 26 years
### -- 4. Detrended data set not older than 1900
### -- 5. No pointer year calculated for detrended data set

###-------------------------------------------------------------------------- *

#---------------------Load libraries-------------------------------------------
library(data.table)
library(dplyr)
library(dplR)
library(pointRes)
library(tibble)
library(furrr)
library(future)

###-------------------------------------------------------------------------- *

# ---- 01 Prepare the data ----
itrdb_in <- file.path("02_data", "01_tree_data",
                      "01_ITRDB_dendroecology")

derived_in <- file.path("02_data", "03_derived_data")

# Remember to load the metadata previously created.
metadata <- read.csv(file.path(derived_in, "metadata_age.csv")) %>%
  tibble()

metadata_filter <- readRDS(file.path(derived_in, "metadata_filter.rds"))

rwl_files <- list.files(file.path(itrdb_in, "rwl"),
                        pattern = "\\.rwl$",
                        full.names = TRUE)

# Filter the rwl files again
rwl_files <- rwl_files[
  tools::file_path_sans_ext(basename(rwl_files)) %in% metadata_filter$rwl
]

##-- Step 2: Calculate ITRDB resilience components--------------------------
## Calculate the resilience components using 'pointRes' package. Loads
## individually each df and do the analysis and then saves it in a general df.
## Note: Modify the (ncol(rwl_data_check) and (nrow(rwl_data_check) thresholds
## depending on the analysis you want to do. They are set in relation to the
## pointer year calculation and the minimum amount of trees to calculate the
## components.


pointeryear_calc <- function(i) {
  file <- rwl_files[i]
  file_name <- tools::file_path_sans_ext(basename(file))

  # Read as rwl (matrix-like: rows=years, cols=series)
  rwl_py <- tryCatch(
    {
      capture.output(
        rwl <- read.rwl(file),
        file = NULL
      )
      rwl
    },
    error = function(e) {
      message("Could not read ", file_name, ": ", e$message)
      NULL
    }
  )

  # Missing rings counts (if your codes are inside the rwl)
  missing_years_0 <- sum(rwl_py == 0, na.rm = TRUE)
  missing_years_2 <- sum(rwl_py == 0.01, na.rm = TRUE)
  missing_years_3 <- sum(rwl_py == 0.001, na.rm = TRUE)

  # Check 1: at least 3 series
  if (ncol(rwl_py) < 3) {
    return(tibble(
  rwl = file_name, pointer_year = NA_real_,
      flag = "less than 3 trees",
      series_count = ncol(rwl_py),
      total_year = nrow(rwl_py),
      missing_years_0 = missing_years_0,
      missing_years_2 = missing_years_2,
      missing_years_3 = missing_years_3,
      total_missing_years = missing_years_0 + missing_years_2 + missing_years_3
    ))
  }

  # Check 2: at least 26 years
  if (nrow(rwl_py) < 26) {
    return(tibble(
      rwl = file_name, pointer_year = NA_real_,
      flag = "Less than 26 years",
      series_count = ncol(rwl_py),
      total_year = nrow(rwl_py),
      missing_years_0 = missing_years_0,
      missing_years_2 = missing_years_2,
      missing_years_3 = missing_years_3,
      total_missing_years = missing_years_0 + missing_years_2 + missing_years_3
    ))
  }

  rwl_py_detrend <- lowpass13(rwl_py)

  # Check 3: Check again after detrend
  if (nrow(rwl_py_detrend) < 26) {
    return(tibble(
      rwl = file_name, pointer_year = NA_real_,
      flag = "Less than 26 years after detrend",
      series_count = ncol(rwl_py_detrend),
      total_year = nrow(rwl_py_detrend),
      missing_years_0 = missing_years_0,
      missing_years_2 = missing_years_2,
      missing_years_3 = missing_years_3,
      total_missing_years = missing_years_0 + missing_years_2 + missing_years_3
    ))
  }

  pointer_years <- pointer.norm(
    rwl_py_detrend,
    series.thresh = 60,
    method.thresh = "Neuwirth",
    make.plot = FALSE,
    window = 13
  )

  py <- as.data.frame(pointer_years$out) %>% filter(year >= 1900)

  # Check 4
  if (nrow(py) == 0) {
    return(tibble(
      rwl = file_name, pointer_year = NA_real_,
      flag = "older than 1900",
      series_count = ncol(rwl_py_detrend),
      total_year = nrow(rwl_py_detrend),
      missing_years_0 = missing_years_0,
      missing_years_2 = missing_years_2,
      missing_years_3 = missing_years_3,
      total_missing_years = missing_years_0 + missing_years_2 + missing_years_3
    ))
  }

  negative_py <- py %>%
    filter(nature == -1) %>%
    dplyr::select(-perc.pos.extreme, -perc.pos.strong, -perc.pos.weak) %>%
    rename(pointer_year = year, py_series = nb.series) %>%
    mutate(
      rwl = file_name,
      series_count = ncol(rwl_py_detrend),
      total_year = nrow(rwl_py_detrend),
      missing_years_0 = missing_years_0,
      missing_years_2 = missing_years_2,
      missing_years_3 = missing_years_3,
      total_missing_years = missing_years_0 + missing_years_2 + missing_years_3
    )

  # Check 5
  if (nrow(negative_py) == 0) {
    return(tibble(
      rwl = file_name, pointer_year = NA_real_,
      flag = "no pointer years calculated",
      series_count = ncol(rwl_py_detrend),
      total_year = nrow(rwl_py_detrend),
      missing_years_0 = missing_years_0,
      missing_years_2 = missing_years_2,
      missing_years_3 = missing_years_3,
      total_missing_years = missing_years_0 + missing_years_2 + missing_years_3
    ))
  }

  negative_py
}

# run in parallel
plan(multisession, workers = max(1, parallel::detectCores() - 1))

itrdb_pointeryear <- future_map_dfr(seq_along(rwl_files), pointeryear_calc,
                                    .options = furrr_options(seed = TRUE))

# Export the pointer years ----
write.csv(x = itrdb_pointeryear, file = file.path("02_data",
                                                  "03_derived_data",
                                                  "itrdb_pointeryear.csv"),
          row.names = FALSE)

###-------------------------------------------------------------------------- *

# Plot the pointer years ----
plots_out <- file.path("03_output", "01_figures")

pointer_plot <- ggplot() +
  geom_bar(data = itrdb_pointeryear, aes(x = pointer_year)) +
  labs(title = "Pointer Year frequency in the ITRDB",
       x = "Year",
       y = "Pointer Years Count") +
  theme_minimal()

ggsave(plot = pointer_plot,
       filename = "pointer_years_itrdb.png",
       path = plots_out, device = "png",
       width = 1920, height = 1080, units = "px")
###-------------------------------------------------------------------------- *
###-------------------------------------------------------------------------- *
###-------------------------------------------------------------------------- *
