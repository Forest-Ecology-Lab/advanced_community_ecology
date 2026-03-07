# ------------------------------- Growth Trends ----------------------------- *
### 2025.03.20: R Script: Calculate growth trends from the ITRDB downloaded
### data. Practice of forest ecology from the course Advance Community Ecology.
###
###
###
### Elaborated by Zabdi López and Rubén D. Manzanedo, Forest Ecology Lab,
### University of Bern

## R script to:
## - 01 Loop to calculate growth trends
## - 02 Generate dataframes
## - 03 Export dataframes

# ---- Load libraries ----
library(dplR)
library(dplyr)
library(tidyr)
library(tibble)

# 01 Prepare the data ----
itrdb_in <- file.path("02_data", "01_tree_data",
                      "01_ITRDB_dendroecology")

metadata <- read.csv(file.path(itrdb_in, "itrdb_site_metadata.csv")) %>%
  tibble() %>%
  select(rwl = studyCode, first_year, last_year, country,
         latitude, longitude, elevation, species_code, species_name) %>%
  mutate(rwl = tolower(rwl))

rwl_files <- list.files(file.path(itrdb_in, "rwl"),
                        pattern = "\\.rwl$",
                        full.names = TRUE)

chron_list <- vector("list", length(rwl_files))
bai_list <- vector("list", length(rwl_files))
flag_list <- vector("list", length(rwl_files))

pb <- txtProgressBar(min = 0, max = length(rwl_files), style = 3)

# 02 Loop to calculate growth trends ----
for (i in seq_along(rwl_files)) {

  rwl_path <- rwl_files[i]
  rwl_name <- gsub("\\.rwl$", "", basename(rwl_path))

  setTxtProgressBar(pb, i)

  rwl <- tryCatch(
    {
      capture.output(
        rwl <- read.rwl(rwl_path)
      )
      rwl
    },
    error = function(e) {
      message("Could not read ", rwl_name, ": ", e$message)
      NULL
    }
  )

  if (is.null(rwl)) {
    flag_list[[i]] <- tibble(
      rwl = rwl_name,
      flag = "reading rwl file failed"
    )
    next
  }

  rwi <- tryCatch(
    {
      detrend(rwl = rwl, method = "Spline")
    },
    error = function(e) {
      message("Could not detrend ", rwl_name, ": ", e$message)
      NULL
    }
  )
  if (is.null(rwi)) {
    flag_list[[i]] <- tibble(
      rwl = rwl_name,
      flag = "detrend failed"
    )
    next
  }

  site_chron <- tryCatch(
    {
      chron(rwi) %>%
        tibble::rownames_to_column("year") %>%
        mutate(year = as.integer(year))
    },
    error = function(e) {
      message("Could not build Chronology ", rwl_name, ": ", e$message)
      NULL
    }
  )
  if (is.null(site_chron)) {
    flag_list[[i]] <- tibble(
      rwl = rwl_name,
      flag = "Chronology failed"
    )
    next
  }

  site_bai_w <- tryCatch(
    {
      bai.out(rwl) %>%
        tibble::rownames_to_column("year") %>%
        mutate(year = as.integer(year))
    },
    error = function(e) {
      message("Could not calculate BAI ", rwl_name, ": ", e$message)
      NULL
    }
  )
  if (is.null(site_bai_w)) {
    flag_list[[i]] <- tibble(
      rwl = rwl_name,
      flag = "BAI failed"
    )
    next
  }

  site_bai <- site_bai_w %>%
    pivot_longer(cols = -year,
                 names_to = "tree",
                 values_to = "bai")

  chron_list[[i]] <- site_chron %>%
    mutate(rwl = rwl_name) %>%
    select(rwl, year, std, samp.depth) %>%
    rename(rwi = std)

  bai_list[[i]] <- site_bai %>%
    mutate(rwl = rwl_name) %>%
    select(rwl, tree, year, bai)

  flag_list[[i]] <- NULL
}
close(pb)

# 03 Generate the data frames ----
chron_df <- bind_rows(chron_list)
bai_df <- bind_rows(bai_list)
flags <- bind_rows(flag_list)

# 04 Export the data frames ----
derived_out <- file.path("02_data", "03_derived_data")

write_csv(x = chron_df, file = file.path(derived_out, "rwi_calculations.csv"))
write_csv(x = bai_df, file = file.path(derived_out, "bai_calculations.csv"))
write_csv(x = flags, file = file.path(derived_out, "growth_trend_flags.csv"))
# --------------------------------------------------------------------------- *
# --------------------------------------------------------------------------- *
# --------------------------------------------------------------------------- *
