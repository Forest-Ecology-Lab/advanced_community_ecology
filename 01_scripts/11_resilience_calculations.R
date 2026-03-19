###-------------------------Resilience components---------------------------- *

### 2025.02.14: R script to calculate the resilience components proposed by *
### Lloret et al. (2011) from the ITRDB data set.
### © Zabdi López, Forest Ecology Lab, Uni. Bern

### R Script to:
### -Load each file and verify it has more that 3 trees and more than 26 years.
### -Convert it into rwl format.
### -#13-year weighted low-pass filter to improve detection (detrend)
### -calculate pointer years (using a normalization in a moving window)
### -calculate the resilience components (resilience, relative resilience,
###  resistance, recovery)
### -calculate missing rings identified as 0.01 and 0.001
### -extract the information for the pointer years
### -Merge it with ITRDB metadata
### -Flags:
### -- 1. Data set with less than 3 series
### -- 2. Data set longer than 26 years
### -- 3. Detrended data set longer than 26 years
### -- 4. Detrended data set not older than 1900
### -- 5. No pointer year calculated for detrended data set


###-------------------------------------------------------------------------- *

# ---- Load libraries ---- *
library(tidyverse)
library(dplyr)
library(dplR)
library(pointRes)
library(tibble)
library(furrr)

# --------------------------------------------------------------------------- *

# ---- Step 1: Set the folders and file names ----
itrdb_in <- file.path("02_data", "01_tree_data", "01_ITRDB_dendroecology",
                      "rwl")
derived_in <- file.path("02_data", "03_derived_data")

rwl_files <- list.files(file.path(itrdb_in),
                        pattern = "\\.rwl$",
                        full.names = TRUE)

metadata <- read.csv(file.path(derived_in, "metadata_age.csv")) %>%
  tibble()

# Load the filter
metadata_filter <- readRDS(file = file.path(derived_in,
                                            "metadata_filter.rds"))

# Filter the rwl files again
rwl_files <- rwl_files[
  tools::file_path_sans_ext(basename(rwl_files)) %in% metadata_filter$rwl
]

rwl_files_names <- basename(rwl_files)

# Load the pointer year
itrdb_py <- read.csv(file = file.path(derived_in, "itrdb_pointeryear.csv"))

cat("Pointer years frequency", capture.output(
  itrdb_py %>%
    count(pointer_year) %>%
    arrange(desc(n))
),
sep = "\n"
  )

# --------------------------------------------------------------------------- *

# YOU CAN MODIFY THIS FILTER AS YOU WANT
## The filter uses the
py_filter <- itrdb_py %>%
  count(pointer_year) %>%
  filter(n > 10)

# Creates a list of only the pointer years
yrs <- sort(unique(py_filter$pointer_year))

# --------------------------------------------------------------------------- *

##-- Step 2: Calculate ITRDB resilience components--------------------------
## Calculate the resilience components using 'pointRes' package. Loads
## individually each rwl calculates the resilience components and the filters
## by pointer years.YOU CAN CHANGE THE POINTER YEARS TO USE ACCORDING TO THE
## QUESTION YOU WANT TO ANSWER.

# Function 1: This function helps to create a mean and sd per year.
mean_sd <- function(x) {
  mat <- as.matrix(x)

  # Replace NA, NaN, Inf, -Inf with NA
  mat[!is.finite(mat)] <- NA_real_

  df <- as.data.frame(mat) %>%
    mutate(
      average = rowMeans(., na.rm = TRUE),
      sd = apply(., 1, sd, na.rm = TRUE)
    )

  df$average[is.nan(df$average)] <- NA_real_
  df$sd[is.nan(df$sd)] <- NA_real_

  df
} # End of function 1

###-------------------------------------------------------------------------- *

# Function 2: Loop to calculate resilience components for each rwl.
process_one_rwl <- function(i) {
  file_rwl <- rwl_files[i]
  file_name <- tools::file_path_sans_ext(basename(file_rwl))

  # Read rwls
  rwl_py <- tryCatch(
    {
      capture.output(
        rwl <- read.rwl(file_rwl),
        file = NULL
      )
      rwl
    },
    error = function(e) {
      message("Could not read ", file_name, ": ", e$message)
      NULL
    }
  )

  # Check 1: at least 3 series
  if (ncol(rwl_py) < 3) {
    return(tibble(
      rwl = file_name, pointer_year = NA_real_,
      resilience = NA_real_, res_sd = NA_real_,
      rel_resilience = NA_real_, relres_sd = NA_real_,
      resistance = NA_real_, resis_sd = NA_real_,
      recovery = NA_real_, recov_sd = NA_real_,
      flag = "less than 3 trees"
    ))
  }

  # Check 2: at least 26 years
  if (nrow(rwl_py) < 26) {
    return(tibble(
      rwl = file_name, pointer_year = NA_real_,
      resilience = NA_real_, res_sd = NA_real_,
      rel_resilience = NA_real_, relres_sd = NA_real_,
      resistance = NA_real_, resis_sd = NA_real_,
      recovery = NA_real_, recov_sd = NA_real_,
      flag = "Less than 26 years"
    ))
  }

  # Detrend the rwl file
  rwl_detrend <- lowpass13(rwl_py)

  # Check 3: Recheck if there are 26 years after detrend
  if (nrow(rwl_detrend) < 26) {
    return(tibble(
      rwl = file_name, pointer_year = NA_real_,
      resilience = NA_real_, res_sd = NA_real_,
      rel_resilience = NA_real_, relres_sd = NA_real_,
      resistance = NA_real_, resis_sd = NA_real_,
      recovery = NA_real_, recov_sd = NA_real_,
      flag = "Less than 26 years after detrend"
    ))
  }

  # Calculate resilience components
  resilience_comp <- res.comp(rwl_detrend, nb.yrs = c(4, 4), max.yrs.rec = 10)

  # Extract the resilience data frames
  ## Resilience
  resil_df <- mean_sd(resilience_comp$resil) %>%
    rownames_to_column("pointer_year") %>%
    mutate(rwl = file_name, pointer_year = as.numeric(pointer_year))

  ## Relative Resilience
  rel_df <- mean_sd(resilience_comp$rel.resil) %>%
    rownames_to_column("pointer_year") %>%
    mutate(rwl = file_name, pointer_year = as.numeric(pointer_year))

  ## Resistance
  resist_df <- mean_sd(resilience_comp$resist) %>%
    rownames_to_column("pointer_year") %>%
    mutate(rwl = file_name, pointer_year = as.numeric(pointer_year))

  ## Recovery
  recov_df <- mean_sd(resilience_comp$recov) %>%
    rownames_to_column("pointer_year") %>%
    mutate(rwl = file_name, pointer_year = as.numeric(pointer_year))

  # Filtered resilience data frames by pointer years
  resilience_pointy <- resil_df %>%
    filter(pointer_year %in% yrs) %>%
    transmute(rwl = file_name, pointer_year, resilience = average,
              resis_sd = sd)

  # Filtered relative resilience data frames by pointer years
  re_resilience_pointy <- rel_df %>%
    filter(pointer_year %in% yrs) %>%
    transmute(rwl = file_name, pointer_year, rel_resilience = average,
              relres_sd = sd)

  # Filtered resistance data frames by pointer years
  resistance_pointy <- resist_df %>%
    filter(pointer_year %in% yrs) %>%
    transmute(rwl = file_name, pointer_year, resistance = average,
              resist_sd = sd)

  # Filtered recovery data frames by pointer years
  recovery_pointy <- recov_df %>%
    filter(pointer_year %in% yrs) %>%
    transmute(rwl = file_name, pointer_year, recovery = average,
              recov_sd = sd)

  # Calculate the number of series used
  res_nb_series <- as.data.frame(resilience_comp$nb.series) %>%
    rownames_to_column("pointer_year") %>%
    mutate(pointer_year = as.numeric(pointer_year)) %>%
    filter(pointer_year %in% yrs) %>%
    dplyr::select(pointer_year, resil, rel.resil, resist, recov) %>%
    rename(
      resilience_series = resil,
      rel_resilience_series = rel.resil,
      resist_series = resist,
      recovery_series = recov
    )

  # Join the data frames together
  resilience_out <- resilience_pointy %>%
    full_join(re_resilience_pointy, by = c("rwl", "pointer_year")) %>%
    full_join(resistance_pointy, by = c("rwl", "pointer_year")) %>%
    full_join(recovery_pointy, by = c("rwl", "pointer_year")) %>%
    left_join(res_nb_series, by = "pointer_year") %>%
    mutate(
      flag = NA_character_
    )

  resilience_out
} # End of function 2

###-------------------------------------------------------------------------- *

# run in parallel
plan(multisession, workers = max(1, parallel::detectCores() - 1))

# Run the function
resilience_df <- future_map_dfr(seq_along(rwl_files), process_one_rwl,
                                .options = furrr_options(seed = TRUE))

# Final Resilience data frame
itrdb_resilience <- resilience_df %>%
  left_join(metadata, by = "rwl")

# Save the resilience info
write.csv(x = itrdb_resilience, file = file.path("02_data", "03_derived_data",
                                                 "itrdb_resilience.csv"),
          row.names = FALSE)

###-------------------------------------------------------------------------- *
###-------------------------------------------------------------------------- *
###-------------------------------------------------------------------------- *
