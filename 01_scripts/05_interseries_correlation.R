# ------------------------- Interseries Correlations ------------------------ *
### 2025.03.20: R Script: Calculate growth trends from the ITRDB downloaded
### data. Practice of forest ecology from the course Advance Community Ecology.
###
###
###
### Elaborated by Zabdi López and Rubén D. Manzanedo, Forest Ecology Lab,
### University of Bern

## R script to:
## - 01 Loop to calculate Interseries Correlations
## - 02 Generate dataframes
## - 03 Export dataframes

# ---- Load libraries ---- *
library(dplR)
library(dplyr)
library(tibble)

# --------------------------------------------------------------------------- *

# ---- 01 Prepare the data ----
itrdb_in <- file.path("02_data", "01_tree_data",
                      "01_ITRDB_dendroecology")
derived_in <- file.path("02_data", "03_derived_data")


# Load the metadata previously created.
metadata <- read.csv(file.path(derived_in, "metadata_age.csv")) %>%
  tibble()

# Load the filter
metadata_filter <- readRDS(file = file.path(derived_in,
                                            "metadata_filter.rds"))

# RWL paths
rwl_files <- list.files(file.path(itrdb_in, "rwl"),
                        pattern = "\\.rwl$",
                        full.names = TRUE)

# Filter the rwl files again
rwl_files <- rwl_files[
  tools::file_path_sans_ext(basename(rwl_files)) %in% metadata_filter$rwl
]

# --------------------------------------------------------------------------- *

# ---- 02 Loop to calculate Interseries Correlations ----

# The loop will read the rwl files, then detrend it to calculate the 
# interseries correlation.
interseries_cor_list <- vector("list", length(rwl_files))
flag_list <- vector("list", length(rwl_files))

# Set the progress bar
pb <- txtProgressBar(min = 0, max = length(rwl_files), style = 3)

# Start the loop
for (i in seq_along(rwl_files)) {

  # Set paths and names of RWL files
  rwl_path <- rwl_files[i]
  rwl_name <- tools::file_path_sans_ext(basename(rwl_path))

  # Start the progress bar
  setTxtProgressBar(pb, i)

  # Read the rwl files
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

  # Calculate the interseries correlation
  interseries_cor <- tryCatch(
    {
      interseries.cor(rwi)
    },
    error = function(e) {
      message("Could not calculate Interseries Correlation ",
              rwl_name, ": ", e$message)
      NULL
    }
  )
  if (is.null(interseries_cor)) {
    flag_list[[i]] <- tibble(
      rwl = rwl_name,
      flag = "Interseries Correlation failed"
    )
    next
  }

  # Save the interseries correlation results in to the list
  interseries_cor_list[[i]] <- interseries_cor %>% 
    mutate(rwl = rwl_name,
           res.cor = interseries_cor$res.cor,
           p.val = p.val) %>% 
    tibble::rownames_to_column("tree")

  # Save the flags into the flag list
  flag_list[[i]] <- NULL
}

# Close the progress bar
close(pb)

# --------------------------------------------------------------------------- *

# ---- 03 Generate the data frames ----
interseries_tree <- bind_rows(interseries_cor_list)
flags <- bind_rows(flag_list)

interseries_site <- interseries_tree %>% 
  group_by(rwl) %>% 
  summarise(mean_cor = mean(res.cor, na.rm = TRUE),
            mean_pval = mean(p.val, na.rm = TRUE))

# --------------------------------------------------------------------------- *

# ---- 04 Export the data frames ----
derived_out <- file.path("02_data", "03_derived_data")

write.csv(x = interseries_tree,
          file = file.path(derived_out, "interseries_tree.csv"),
          row.names = FALSE)
write.csv(x = interseries_site,
          file = file.path(derived_out, "interseries_site.csv"),
          row.names = FALSE)
write.csv(x = flags,
          file = file.path(derived_out, "interseries_cor_flags.csv"),
          row.names = FALSE)
# --------------------------------------------------------------------------- *
# --------------------------------------------------------------------------- *
# --------------------------------------------------------------------------- *
