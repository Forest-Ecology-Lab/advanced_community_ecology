# ----------------------------- Climate Correlations ------------------------ *
### 2025.03.20: R Script: Calculate the SPEI Correlations per site for the
### practice of Forest Ecology from the course Advance Community Ecology.
###
###
### Elaborated by Zabdi López and Rubén D. Manzanedo, Forest Ecology Lab,
### University of Bern

## R script to download:
## - 01 Prepare the data
## - 02 Calculate Climate Correlations
## - 03 Plot SPEI Correlations:
## ---- 01: Mean SPEI Correlation
## ---- 02: Strongest SPEI Month Correlation
## ---- 03: Monthly SPEI Correlation
## ---- 04: SPEI BoxPlot
## ---- 05: Significant Sites by Month
## - 03 Export maps and plots

# ---- Load the Libraries ---- *
library(dplyr)
library(tidyverse)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(treeclim)

# --------------------------------------------------------------------------- *
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

# Remember to load the metadata previously created.
metadata <- read.csv(file.path(derived_in, "metadata_age.csv")) %>%
  tibble()

# ---- Load the chronologies ----- *
itrdb_rwi <- readRDS(file.path(derived_in, "chron_list.rds"))

# ---- Load SPEI Data ---- *

itrdb_spei <- read.csv(file.path(derived_in, "spei6month_itrdb.csv"),
                       check.names = FALSE)

## Convert to long format
spei_l <- itrdb_spei %>%
  pivot_longer(-rwl,
               names_to = "date",
               values_to = "spei") %>%
  separate(date, into = c("year", "month"), sep = "_") %>%
  mutate(year = as.integer(year))

# --------------------------------------------------------------------------- *
# --------------------------------------------------------------------------- *

# ---- 02 Loop of climate correlations through sites ----

# Set variables for loop
site_names <- names(itrdb_rwi)

spei_list <- vector("list", length(site_names))

# Loop through sites and variables
pb <- txtProgressBar(min = 0, max = length(site_names), style = 3)

for (i in seq_along(site_names)) {

  # This will get the site that the iteration will be working with
  site <- site_names[i]

  # This is where the results from the iteration will be stored
  out <- tryCatch({

    # Extract the chronology from the list previously created and put it in the
    # correct format.
    chrono_site <- itrdb_rwi[[site]] %>%
      select(year, std) %>%
      as.data.frame() %>%
      column_to_rownames("year")

    # Extract the SPEI data for that specific site
    spei_site <- spei_l %>%
      filter(rwl == site) %>%
      select(year, month, spei) %>%
      mutate(month = match(month, month.abb)) %>%
      as.data.frame()

    # Calculates the corelation of the site SPEI and Chronology
    resp <- dcc(
                chrono = chrono_site,
                climate = spei_site,
                selection = -6:9,
                verbose = FALSE)


    # Saves it in a data frame with the correlation information.
    data.frame(rwl = site,
               clim_var = "spei",
               window = rownames(resp$coef),
               coef = resp$coef$coef,
               significant = as.logical(resp$coef$significant),
               month = resp$coef$month,
               stringsAsFactors = FALSE)

    # This will kick in if there is an error. If there is an error it will
    # store it with NAs so we can know which site gave error
  }, error = function(e) {

    data.frame(rwl = site,
      clim_var = "spei",
      window = NA_character_,
      coef = NA_real_,
      significant = NA,
      month = NA,
      stringsAsFactors = FALSE
    )
  })

  # At the end the data frame will be stored inside the list.
  spei_list[[i]] <- out
  setTxtProgressBar(pb, i)

} # end of loop

# Close the progress bar
close(pb)

# Creates a data frame using all the data calculated in the loop.
spei_coef_df <- bind_rows(spei_list)

# Clean the data and set up the levels
spei_coefs <- spei_coef_df %>%
  mutate(
    month_plot = ifelse(
      grepl("prev", window),
      paste0("p", month),
      month
    ),
    month_plot = factor(
      month_plot,
      levels = c("pJun", "pJul", "pAug", "pSep", "pOct", "pNov", "pDec",
                 "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP")
    )
  ) %>%
  filter(!is.na(coef))

# Export the SPEI Coefficients in case you need them later

saveRDS(object = spei_coefs,
        file = file.path(derived_in, "spei_coefs.rds"))

# --------------------------------------------------------------------------- *
# --------------------------------------------------------------------------- *

# ---- 03 Plot the data ----
spei_plot <- spei_coefs %>%
  select(rwl, month_plot, clim_var, coef, significant) %>%
  left_join(metadata, by = "rwl") %>%
  filter(!is.na(latitude), !is.na(longitude))

# === 01 Overall Mean Map === *

## Calculate the mean SPEI coefficients per month and site
spei_mean_map <- spei_plot %>%
  group_by(rwl, latitude, longitude) %>%
  summarise(
    mean_coef = mean(coef, na.rm = TRUE),
    n_months = sum(!is.na(coef)),
    .groups = "drop"
  )

## Plot it
overall_mean_map <- ggplot() +
  geom_sf(data = europe_map, fill = "grey95",
          colour = "grey70", linewidth = 0.2) +
  geom_point(
    data = spei_mean_map,
    aes(x = longitude, y = latitude, colour = mean_coef),
    size = 2.2,
    alpha = 0.9
  ) +
  scale_colour_gradient2(
    low = "brown",
    mid = "white",
    high = "darkgreen",
    midpoint = 0
  ) +
  coord_sf(xlim = europe_limits$xlim,
           ylim = europe_limits$ylim,
           expand = FALSE) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Overall mean SPEI-growth correlation",
    subtitle = "Mean coefficient across all months",
    x = "Longitude",
    y = "Latitude",
    colour = "Mean\ncorrelation"
  )

## Check it
overall_mean_map

# -------- 02 Strongest Month Per Site ------- *

## Calculate the strongest month per site
spei_strongest <- spei_plot %>%
  filter(clim_var == "spei") %>%
  group_by(rwl) %>%
  slice_max(order_by = abs(coef), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(cor_sign = ifelse(coef >= 0, "Positive", "Negative"))

## Set a manual color gradient
month_cols <- c(
  "#3B4CC0",  # pJun
  "#4A6FE3",  # pJul
  "#6A8BE8",  # pAug
  "#8FA4F4",  # pSep
  "#B2B7F8",  # pOct
  "#C8C8F9",  # pNov
  "#DCDCFB",  # pDec
  "#D0F0FF",  # Jan
  "#B5E3FF",  # Feb
  "#88D498",  # Mar
  "#6BCB77",  # Apr
  "#4CAF50",  # May
  "#FFB347",  # Jun
  "#FF8C42",  # Jul
  "#FF6B6B",  # Aug
  "#D7263D"   # Sep
)

## Plot it
strongest_month_map <- ggplot() +
  geom_sf(data = europe_map, fill = "grey95",
          colour = "grey70", linewidth = 0.2) +
  geom_point(
    data = spei_strongest,
    aes(
      x = longitude,
      y = latitude,
      fill = month_plot,
      shape = cor_sign
    ),
    colour = "black",
    size = 3,
    stroke = 0.3
  ) +
  scale_shape_manual(
    values = c(
      "Positive" = 24,
      "Negative" = 25
    )
  ) +
  scale_fill_manual(
    values = month_cols,
    name = "Month"
  ) +
  coord_sf(xlim = europe_limits$xlim,
           ylim = europe_limits$ylim,
           expand = FALSE) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21,
        colour = "black",
        size = 4
      ),
      ncol = 2
    ),
    shape = guide_legend(
      override.aes = list(
        fill = "white",
        colour = "grey40",
        size = 4
      )
    )
  ) +
  labs(
    title = "Month of strongest SPEI-growth correlation by site",
    subtitle = "Triangle up = positive correlation, triangle down = negative",
    x = "Longitude",
    y = "Latitude",
    shape = "Correlation"
  )

## Check it
strongest_month_map

# ------- 03 Faceted Monthly Maps -------- *
spei_map <- ggplot() +
  geom_sf(data = europe_map, fill = "grey95",
          colour = "grey70", linewidth = 0.2) +
  geom_point(
             data = spei_plot %>% filter(!significant),
             aes(x = longitude, y = latitude, colour = coef),
             shape = 1,
             size = 1.8,
             alpha = 0.8) +

  geom_point(
             data = spei_plot %>% filter(significant),
             aes(x = longitude, y = latitude, colour = coef),
             shape = 19,
             size = 2.2,
             alpha = 0.9) +

  scale_colour_gradient2(
                         low = "brown",
                         mid = "white",
                         high = "darkgreen",
                         midpoint = 0) +

  coord_sf(xlim = europe_limits$xlim,
           ylim = europe_limits$ylim,
           expand = FALSE) +

  facet_wrap(~ month_plot, ncol = 4) +

  labs(
       ## Update the tittle accordingly to your research question!
       title = "Spatial pattern of SPEI-growth correlations for PISY",
       subtitle = "Filled points are significant",
       x = "Longitude",
       y = "Latitude",
       colour = "coef") +

  theme_minimal() +
  theme(
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey90", colour = NA),
        strip.text = element_text(face = "bold"),
        legend.position = "right")

## Check it
spei_map

# -------- 04 Boxplot by Month -------- *

boxplot_month <- ggplot(spei_plot, aes(x = month_plot, y = coef)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_boxplot(outlier.alpha = 0.35) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Distribution of SPEI-growth correlations by month",
    x = "Month",
    y = "Correlation coefficient"
  )

## Check it
boxplot_month

# -------- 05 Significant sites by month -------- *

## Calculate the amount of signigicant sites per month
sig_sites_month_df <- spei_plot %>%
  group_by(month_plot) %>%
  summarise(
    n_sig = sum(significant, na.rm = TRUE),
    n_total = sum(!is.na(significant)),
    prop_sig = n_sig / n_total,
    .groups = "drop"
  )

## Plot it
sig_sites_month_plot <- ggplot(sig_sites_month_df,
                               aes(x = month_plot, y = n_sig)) +
  geom_col() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Number of significant SPEI-growth correlations by month",
    x = "Month",
    y = "Number of significant sites"
  )

## Check it
sig_sites_month_plot

# -------- 04 Export the maps and figures --------
maps_out <- file.path("03_output", "02_maps")
plots_out <- file.path("03_output", "01_figures")

ggsave(plot = overall_mean_map, filename = "mean_spei_map.png",
       path = maps_out, device = "png", dpi = 300,
       width = 4500, height = 4500, units = "px")

ggsave(plot = strongest_month_map, filename = "strong_month_spei_map.png",
       path = maps_out, device = "png", dpi = 300,
       width = 4500, height = 4500, units = "px")

ggsave(plot = spei_map, filename = "monthly_spei_map.png",
       path = maps_out, device = "png", dpi = 300,
       width = 4500, height = 4500, units = "px")

ggsave(plot = boxplot_month, filename = "spei_boxplot.png",
       path = plots_out, device = "png", dpi = 300,
       width = 1920, height = 1080, units = "px")

ggsave(plot = sig_sites_month_plot, filename = "sig_sites_spei.png",
       path = plots_out, device = "png", dpi = 300,
       width = 1920, height = 1080, units = "px")
# --------------------------------------------------------------------------- *
# --------------------------------------------------------------------------- *
# --------------------------------------------------------------------------- *
