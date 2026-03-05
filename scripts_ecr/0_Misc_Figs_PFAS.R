library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)

setwd("~/Desktop/EagleStats/")


df <- read.csv("./input/2023_2024_2025_AllData_clean.csv")
key <- read.csv("./input/PFAS_subtypes.csv")

# --- 1. Load Wisconsin outline ---
wisconsin <- ne_states(country = "United States of America", returnclass = "sf") %>%
  dplyr::filter(name == "Wisconsin")

# --- 2. Make sure df_sub is an sf object (points) ---
df_points <- df_sub %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326)  # WGS84

# --- 3. Plot ---
ggplot() +
  geom_sf(data = wisconsin, fill = "gray95", color = "black") +
  geom_sf(data = df_points, aes(color = TOTAL_PFAS), size = 3, alpha = 0.8) +
  scale_color_viridis(option = "plasma", trans = "log1p", name = "Total PFAS") +
  labs(
    title = "Eagle Sample Locations in Wisconsin Colored by Total PFAS"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

