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

### Last update: 2026.02.23 ### Linted code

###-------------------------------------------------------------------------- *

#---------------------Load libraries-------------------------------------------
library(data.table)
library(dplyr)
library(dplR)
library(pointRes)
library(tibble)
library(furrr)
library(future)

##-- Step 1: Set the folders and file names------------------------------------
rwl_data_input <- file.path("02_data", "01_EU_ITRDB_dendroecology", "rwl")
rwl_files <- list.files(path = rwl_data_input,
                        pattern = "\\.rwl$", full.names = TRUE)
rwl_files_names <- basename(rwl_files)

##-- Step 2: Calculate ITRDB resilience components--------------------------
## Calculate the resilience components using 'pointRes' package. Loads
## individually each df and do the analysis and then saves it in a general df.
## Note: Modify the (ncol(rwl_data_check) and (nrow(rwl_data_check) thresholds
## depending on the analysis you want to do. They are set in relation to the
## pointer year calculation and the minimum amount of trees to calculate the
## components.


pointeryear_calc <- function(i) {
  file <- rwl_files[i]
  file_name <- gsub("\\.rwl$", "", basename(file))
  
  # Read as rwl (matrix-like: rows=years, cols=series)
  rwl_py <- read.rwl(file)
  
  # Missing rings counts (if your codes are inside the rwl)
  missing_years_0 <- sum(rwl_py == 0, na.rm = TRUE)
  missing_years_2 <- sum(rwl_py == 0.01, na.rm = TRUE)
  missing_years_3 <- sum(rwl_py == 0.001, na.rm = TRUE)
  
  # Check 1: at least 3 series
  if (ncol(rwl_py) < 3) {
    return(tibble(
      rwl = file_name, pointer_year = NA_real_,
      resilience = NA_real_, res_sd = NA_real_,
      rel_resilience = NA_real_, relres_sd = NA_real_,
      resistance = NA_real_, resis_sd = NA_real_,
      recovery = NA_real_, recov_sd = NA_real_,
      flag = "less than 3 trees",
      series_count = ncol(rwl_py),
      total_year = nrow(rwl_py),
      missing_years_2 = missing_years_2,
      missing_years_3 = missing_years_3,
      total_missing_years = missing_years_2 + missing_years_3
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
      flag = "Less than 26 years",
      series_count = ncol(rwl_py),
      total_year = nrow(rwl_py),
      missing_years_2 = missing_years_2,
      missing_years_3 = missing_years_3,
      total_missing_years = missing_years_2 + missing_years_3
    ))
  }
  
  rwl_py_detrend <- lowpass13(rwl_py)
  
  # Check 3
  if (nrow(rwl_py_detrend) < 26) {
    return(tibble(
      rwl = file_name, pointer_year = NA_real_,
      resilience = NA_real_, res_sd = NA_real_,
      rel_resilience = NA_real_, relres_sd = NA_real_,
      resistance = NA_real_, resis_sd = NA_real_,
      recovery = NA_real_, recov_sd = NA_real_,
      flag = "Less than 26 years after detrend",
      series_count = ncol(rwl_py_detrend),
      total_year = nrow(rwl_py_detrend),
      missing_years_2 = missing_years_2,
      missing_years_3 = missing_years_3,
      total_missing_years = missing_years_2 + missing_years_3
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
      resilience = NA_real_, res_sd = NA_real_,
      rel_resilience = NA_real_, relres_sd = NA_real_,
      resistance = NA_real_, resis_sd = NA_real_,
      recovery = NA_real_, recov_sd = NA_real_,
      flag = "older than 1900",
      series_count = ncol(rwl_py_detrend),
      total_year = nrow(rwl_py_detrend),
      missing_years_2 = missing_years_2,
      missing_years_3 = missing_years_3,
      total_missing_years = missing_years_2 + missing_years_3
    ))
  }
  
  negative_py <- py %>%
    filter(nature == -1) %>%
    dplyr::select(-perc.pos.extreme, -perc.pos.strong, -perc.pos.weak) %>%
    rename(pointer_year = year, py_series = nb.series) %>%
    mutate(
      rwl = file_name,
      missing_years_2 = missing_years_2,
      missing_years_3 = missing_years_3,
      total_missing_years = missing_years_2 + missing_years_3
    )
  
  # Check 5
  if (nrow(negative_py) == 0) {
    return(tibble(
      rwl = file_name, pointer_year = NA_real_,
      resilience = NA_real_, res_sd = NA_real_,
      rel_resilience = NA_real_, relres_sd = NA_real_,
      resistance = NA_real_, resis_sd = NA_real_,
      recovery = NA_real_, recov_sd = NA_real_,
      flag = "no pointer years calculated",
      series_count = ncol(rwl_py_detrend),
      total_year = nrow(rwl_py_detrend),
      missing_years_2 = missing_years_2,
      missing_years_3 = missing_years_3,
      total_missing_years = missing_years_2 + missing_years_3
    ))
  }

  
  out
}

# run in parallel
plan(multisession, workers = max(1, parallel::detectCores() - 1))

itrdb_pointeryear <- future_map_dfr(seq_along(rwl_files), pointeryear_calc)

#Check the data frame
View(itrdb_pointeryear)

#Save the resilience info
write.csv(x = itrdb_pointeryear, file = file.path("03_data",
                                                    "04_resilience_analysis",
                                                    "itrdb_pointeryear.csv"))

###-------------------------------------------------------------------------- *
###-------------------------------------------------------------------------- *
###-------------------------------------------------------------------------- *
