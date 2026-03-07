# -------------------- Chronology Coefficient of Variation ------------------ *
### 2025.03.20: R Script: Calculate the chronologies coefficient of variation
### from the ITRDB downloaded data. Practice of forest ecology from the course
### Advance Community Ecology.
###
###
### Elaborated by Zabdi López and Rubén D. Manzanedo, Forest Ecology Lab,
### University of Bern

## R script to:
## - 01 Calculate Chronology Coefficient of Variation
## - 02 Export data frames

# ---- Load libraries ---- *
library(dplyr)

# ---- 01 Prepare data ----
derived_in <- file.path("02_data", "03_derived_data")

rwi_calc <- read.csv(file.path(derived_in, "rwi_calculations.csv"))

# ---- 02 Calculate Chronology CV ----
chron_cv <- rwi_calc %>%
  group_by(rwl) %>%
  summarise(sd = sd(rwi, na.rm = TRUE),
            cv = sd(rwi, na.rm = TRUE) / mean(rwi, na.rm = TRUE),
            .groups = "drop")

# ---- 03 Export data frame ----
write_csv(x = chron_cv, file = file.path(derived_in, "chronology_cv.csv"))
# --------------------------------------------------------------------------- *
# --------------------------------------------------------------------------- *
# --------------------------------------------------------------------------- *