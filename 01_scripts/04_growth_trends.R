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

# --------------------------------------------------------------------------- *

# ---- Load libraries ----
library(dplR)
library(dplyr)
library(tidyr)
library(tibble)

# --------------------------------------------------------------------------- *

# 01 Prepare the data ----
itrdb_in <- file.path("02_data", "01_tree_data",
                      "01_ITRDB_dendroecology")
derived_in <- file.path("02_data", "03_derived_data")

# LOAD THE FILTERED METADATA
metadata <- read.csv(file.path(derived_in, "metadata_age.csv")) %>%
  tibble()

# LOAD THE PREVIOUS FILTER
metadata_filter <- readRDS(file = file.path(derived_in,
                                            "metadata_filter.rds"))

# Set the path of the RWL files
rwl_files <- list.files(file.path(itrdb_in, "rwl"),
                        pattern = "\\.rwl$",
                        full.names = TRUE)

# Filter the rwl files again
rwl_files <- rwl_files[
  tools::file_path_sans_ext(basename(rwl_files)) %in% metadata_filter$rwl
]

# --------------------------------------------------------------------------- *
# 02 Loop to calculate growth trends ----

# This loop will read the RWL file, then proceed to detrend the RWL,then it will
# generate the chronology. After that it will proceed to calculate the BAI
# growth. We use bai.out() because we are calculating it from in to out.

# Create the empty list that will be filled by the loop.
chron_list <- vector("list", length(rwl_files))
bai_list <- vector("list", length(rwl_files))
flag_list <- vector("list", length(rwl_files))

# Set up the progress bar
pb <- txtProgressBar(min = 0, max = length(rwl_files), style = 3)

# Loop starts
for (i in seq_along(rwl_files)) {

  # set rwl paths and names
  rwl_path <- rwl_files[i]
  rwl_name <- tools::file_path_sans_ext(basename(rwl_path))

  # Set the progress bar
  setTxtProgressBar(pb, i)

  # Read the RWL file
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
    flag_list[[i]] <- tibble(
      rwl = rwl_name,
      flag = "reading rwl file failed"
    )
    next
  }

  # Detrend the RWL
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

  # Calculate the chronology from the site
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

  # Calculate the BAI growth
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

  # Convert to long format
  site_bai <- site_bai_w %>%
    pivot_longer(cols = -year,
                 names_to = "tree",
                 values_to = "bai")

  # Save the chronology into the list
  chron_list[[i]] <- site_chron %>%
    mutate(rwl = rwl_name) %>%
    select(rwl, year, std, samp.depth)

  # Save the BAI into the list 
  bai_list[[i]] <- site_bai %>%
    mutate(rwl = rwl_name) %>%
    select(rwl, tree, year, bai)

  flag_list[[i]] <- NULL
} # Loop end

# Close the progress bar
close(pb)
# --------------------------------------------------------------------------- *

# 03 Generate the data frames ----
# Bind the list into Dataframes
chron_df <- bind_rows(chron_list)
bai_df <- bind_rows(bai_list)
flags <- bind_rows(flag_list)

# Name the chronologies in the list
chron_list <- chron_list[!vapply(chron_list, is.null, logical(1))]
bai_list   <- bai_list[!vapply(bai_list, is.null, logical(1))]

names(chron_list) <- vapply(chron_list, function(x) unique(x$rwl), character(1))

# --------------------------------------------------------------------------- *
# 04 Export the data frames ----
derived_out <- file.path("02_data", "03_derived_data")

# Chronologies as R object (It is faster to load)
saveRDS(chron_list, file = file.path(derived_out, "chron_list.rds"))

# Chronologies and other as csv files
write.csv(x = chron_df, file = file.path(derived_out, "rwi_calculations.csv"),
          row.names = FALSE)
write.csv(x = bai_df, file = file.path(derived_out, "bai_calculations.csv"),
          row.names = FALSE)
write.csv(x = flags, file = file.path(derived_out, "growth_trend_flags.csv"),
          row.names = FALSE)
# --------------------------------------------------------------------------- *
# --------------------------------------------------------------------------- *
# --------------------------------------------------------------------------- *
