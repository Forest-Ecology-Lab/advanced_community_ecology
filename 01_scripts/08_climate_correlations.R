# ----------------------------- Climate Correlations ------------------------ *
### 2025.03.20: R Script: Extracting climate variables used in the
### practice of Forest Ecology from the course Advance Community Ecology.
###
###
### Elaborated by Zabdi López and Rubén D. Manzanedo, Forest Ecology Lab,
### University of Bern

## R script to download:
## - 01 Prepare the data
## - 02 Calculate Climate Correlations
## - 03 Export maps

# ---- Load libraries ----
library(dplyr)
library(tibble)
library(purrr)
library(stringr)
library(tidyverse)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(treeclim)

# --------------------------------------------------------------------------- *

# ---- 01 Load Data ----
# ----Map of Europe ---- *
## Load the map
europe_map <- ne_countries(scale = "large",
                           continent = "Europe",
                           returnclass = "sf")

## Set Europe limits
europe_limits <- list(xlim = c(-10, 35), ylim = c(35, 75))

# ---- Load ITRDB filtered Metadata ---- *
derived_in <- file.path("02_data", "03_derived_data")

# Load the metadata previously created
metadata <- read.csv(file.path(derived_in, "metadata_age.csv")) %>%
  tibble()

# ---- Load the chronologies ----- *
itrdb_rwi <- readRDS(file.path(derived_in, "chron_list.rds"))

# ---- Load the Climate Data ----- *
itrdb_pre <- read.csv(file.path(derived_in, "precipitation_itrdb.csv"),
                      check.names = FALSE)
itrdb_tmp <- read.csv(file.path(derived_in, "temperature_itrdb.csv"),
                      check.names = FALSE)
itrdb_tmn <- read.csv(file.path(derived_in, "mintemp_itrdb.csv"),
                      check.names = FALSE)
itrdb_tmx <- read.csv(file.path(derived_in, "maxtemp_itrdb.csv"),
                      check.names = FALSE)

## Convert to long format

### Precipitation
pre_l <- itrdb_pre %>%
  pivot_longer(-rwl,
               names_to = "date",
               values_to = "pre") %>%
  separate(date, into = c("year", "month"), sep = "_") %>%
  mutate(year = as.integer(year))

### Temperature
tmp_l <- itrdb_tmp %>%
  pivot_longer(-rwl,
               names_to = "date",
               values_to = "tmp") %>%
  separate(date, into = c("year", "month"), sep = "_") %>%
  mutate(year = as.integer(year))

### Minimum Temperature
tmn_l <- itrdb_tmn %>%
  pivot_longer(-rwl,
               names_to = "date",
               values_to = "tmn") %>%
  separate(date, into = c("year", "month"), sep = "_") %>%
  mutate(year = as.integer(year))

### Maximum Temperature
tmx_l <- itrdb_tmx %>%
  pivot_longer(-rwl,
               names_to = "date",
               values_to = "tmx") %>%
  separate(date, into = c("year", "month"), sep = "_") %>%
  mutate(year = as.integer(year))

### Join all together
clim_l <-  list(pre_l, tmp_l, tmn_l, tmx_l) %>%
  reduce(left_join, by = c("rwl", "year", "month"))

# --------------------------------------------------------------------------- *

# ---- 02 Loop of climate correlations through sites ----
# This loop will take the prepare the climate in the correct format, then select
# the chronology data. Once it does that, it will run the dcc() which will 
# calculate the correlations.

# Set variables for loop
site_names <- names(itrdb_rwi)
dcc_results <- list()

# Loop through sites and variables
pb <- txtProgressBar(min = 0, max = length(site_names), style = 3)

# Start of the loop
for (i in seq_along(site_names)) {
  
  # Get the site names
  site <- site_names[i]
  
  # Select the site rwi
  chrono_site <- itrdb_rwi[[site]] %>%
    select(year, std) %>%
    as.data.frame() %>%
    column_to_rownames("year")
  
  # Calculate correlations
  dcc_results[[site]] <- tryCatch({
    
    # Set the climate data in the correct format
    climate_site <- clim_l %>%
      filter(rwl == site) %>%
      select(-rwl) %>%
      mutate(month = match(month, month.abb)) %>%
      as.data.frame()
    
    # Calculate the climate correlations for all variables and save it on the
    # list.
    dcc(
      chrono = chrono_site,
      climate = climate_site,
      selection = -6:9,
      verbose = FALSE)
    
  }, error = function(e) {
    message("Failed site: ", site)
    NULL
  })
  setTxtProgressBar(pb, i)
}

close(pb)

# --------------------------------------------------------------------------- *

# ---- Save the data frames ---- *

# Extract the information from the list
clim_coefs <- bind_rows(
  lapply(names(dcc_results), function(site) {
    
    resp <- dcc_results[[site]]
    
    if (is.null(resp) || is.null(resp$coef)) return(NULL)
    
    resp$coef %>%
      rownames_to_column(var = "window") %>%
      mutate(rwl = site) %>%
      select(rwl, window, id, varname, month, coef,
             significant, ci_lower, ci_upper)
  })
)

# Export the data frame
derived_out <- file.path("02_data", "03_derived_data")

write.csv(x = clim_l, file = file.path(derived_out, "clim_long.csv"),
          row.names = FALSE)

clim_coefs <- clim_coefs %>% 
  mutate(month = factor(month, levels = c("Jun", "Jul", "Aug", "Sep",
                                          "Oct", "Nov", "Dec",
                                          "JAN", "FEB", "MAR", "APR",
                                          "MAY", "JUN", "JUL", "AUG",
                                          "SEP")))

saveRDS(object = clim_coefs,
        file = file.path(derived_in, "clim_coefs.rds"))

# ---- 03 Plot the data ----
clim_plot <- clim_coefs %>%
  select(rwl, month, varname, coef, significant,ci_lower,ci_upper) %>%
  left_join(metadata, by = "rwl")


# mean line and proportion of significant sites
clim_mean <- clim_plot %>%
  group_by(varname, month) %>%
  summarise(
    mean_coef = mean(coef, na.rm = TRUE),
    prop_sig = mean(significant, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(mark_sig = prop_sig >= 0.5)

# REMEMBER TO UPDATE THE LABELS CORRECTLY
# This is an easy plot you can do
plot(dcc_results$swed307)

# You can plot the site correlations as this:
site_plot <- clim_plot %>%
  filter(rwl == "swed315")

# Do the plot with lines
ggplot(site_plot, aes(x = month, y = coef, colour = varname)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4) +
  geom_errorbar(
    aes(
      ymin = ci_lower,
      ymax = ci_upper,
      linetype = significant
    ),
    width = 0.15,
    linewidth = 0.7
  ) +
  geom_point(size = 2.8) +
  facet_wrap(~ varname) +
  scale_linetype_manual(
    name = "Significant",
    values = c("FALSE" = "dashed", "TRUE" = "solid"),
    labels = c("FALSE" = "No", "TRUE" = "Yes")) +
  labs(
    x = "Months",
    y = "Coefficients"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

# Now with bars
ggplot(site_plot, aes(x = month, y = coef, fill = varname)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    colour = "grey20"
  ) +
  
  geom_errorbar(
    aes(
      ymin = ci_lower,
      ymax = ci_upper,
      linetype = significant,
      group = varname
    ),
    position = position_dodge(width = 0.8),
    width = 0.2,
    linewidth = 0.6
  ) +
  
  scale_fill_manual(
    values = c(
      "pre" = "#4575b4",  # blue (precip)
      "tmn" = "#fdae61",  # light orange (min temp)
      "tmp" = "#f46d43",  # orange-red (mean temp)
      "tmx" = "#d73027"   # red (max temp)
    ),
    name = "Variable",
    labels = c(
      "pre" = "Precipitation",
      "tmn" = "Min Temp",
      "tmp" = "Mean Temp",
      "tmx" = "Max Temp"
    )
  ) +
  
  scale_linetype_manual(
    name = "Significant",
    values = c("FALSE" = "dashed", "TRUE" = "solid"),
    labels = c("FALSE" = "No", "TRUE" = "Yes")
  ) +
  
  labs(
    x = "Months",
    y = "Coefficients"
  ) +
  
  theme_minimal() +
  facet_wrap(~varname) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# A plot with all the sites on the background and the mean on top of it
cc_plot <- ggplot() +
  geom_line(
    data = clim_plot,
    aes(x = month, y = coef, group = rwl),
    colour = "grey70",
    alpha = 0.08,
    linewidth = 0.3
  ) +
  geom_line(
    data = clim_mean,
    aes(x = month, y = mean_coef, group = 1),
    colour = "black",
    linewidth = 1
  ) +
  geom_point(
    data = clim_mean,
    aes(x = month, y = mean_coef),
    colour = "black",
    size = 1.8
  ) +
  geom_point(
    data = clim_mean %>% filter(mark_sig),
    aes(x = month, y = mean_coef),
    colour = "red",
    size = 3
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ varname, scales = "free_y") +
  labs(
    x = "Months",
    y = "Climate correlation",
    title = "Climate correlation signature across sites for PISY"
  ) +
  theme_minimal(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90))

###-------------------------------------------------------------------------- *

# 04 Export the plots ----
maps_out <- file.path("03_output", "02_maps")
plots_out <- file.path("03_output", "01_figures")

ggsave(plot = cc_plot, filename = "climate_var_correlations.png",
       path = plots_out, device = "png", dpi = 300,
       width = 4500, height = 4500, units = "px")
