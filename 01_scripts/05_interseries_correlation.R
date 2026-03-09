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

# ---- 01 Prepare the data ----
itrdb_in <- file.path("02_data", "01_tree_data",
                      "01_ITRDB_dendroecology")

rwl_files <- list.files(file.path(itrdb_in, "rwl"),
                        pattern = "\\.rwl$",
                        full.names = TRUE)

interseries_cor_list <- vector("list", length(rwl_files))
flag_list <- vector("list", length(rwl_files))


pb <- txtProgressBar(min = 0, max = length(rwl_files), style = 3)

# ---- 02 Loop to calculate Interseries Correlations ----
for (i in seq_along(rwl_files)) {

  rwl_path <- rwl_files[i]
  rwl_name <- tools::file_path_sans_ext(basename(rwl_path))

  setTxtProgressBar(pb, i)

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

  mean_interseries_cor <- tryCatch(
    {
      mean(interseries_cor$res.cor, na.rm = TRUE)
    },
    error = function(e) {
      message("Could not calculate Interseries Correlation mean ",
              rwl_name, ": ", e$message)
      NULL
    }
  )
  if (is.null(mean_interseries_cor)) {
    flag_list[[i]] <- tibble(
      rwl = rwl_name,
      flag = "Interseries cor failed"
    )
    next
  }

  interseries_cor_list[[i]] <- tibble(
                                      rwl = rwl_name,
                                      mean_inter_cor = mean_interseries_cor)

  flag_list[[i]] <- NULL
}
close(pb)

# ---- 03 Generate the data frames ----
interseries_df <- bind_rows(interseries_cor_list)
flags <- bind_rows(flag_list)

# ---- 04 Export the data frames ----
derived_out <- file.path("02_data", "03_derived_data")

write.csv(x = interseries_df,
          file = file.path(derived_out, "mean_interseries.csv"),
          row.names = FALSE)
write.csv(x = flags,
          file = file.path(derived_out, "interseries_cor_flags.csv"),
          row.names = FALSE)
# --------------------------------------------------------------------------- *
# --------------------------------------------------------------------------- *
# --------------------------------------------------------------------------- *
