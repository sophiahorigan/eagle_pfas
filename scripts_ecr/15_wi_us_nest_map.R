##########################################
# Map: Wisconsin nest locations with US inset
# - Point color: year
# - Point size: number of nestlings per nest (nest_no x year)
##########################################

rm(list = ls())

library(ggplot2)
library(dplyr)
library(jsonlite)
library(cowplot)

input_file <- "./input/2023_2024_2025_AllData_clean.csv"
geo_file <- "./input/us-states.geojson"
county_geo_file <- "./input/us-counties-fips.geojson"
out_dir <- "./output/misc_figs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

df <- read.csv(input_file, stringsAsFactors = FALSE)
if (!file.exists(geo_file)) stop("Missing ", geo_file)
if (!file.exists(county_geo_file)) stop("Missing ", county_geo_file)

# Correct obvious sign errors: Wisconsin longitudes should be west (negative)
if (any(df$long > 0, na.rm = TRUE)) {
  df$long[df$long > 0] <- -df$long[df$long > 0]
}

# One point per nest-year, colored by mean TOTAL_PFAS
nest_pts <- df %>%
  filter(!is.na(nest_no), !is.na(year), !is.na(lat), !is.na(long), !is.na(TOTAL_PFAS)) %>%
  group_by(year, nest_no) %>%
  summarize(
    lat = mean(lat, na.rm = TRUE),
    long = mean(long, na.rm = TRUE),
    total_pfas = mean(TOTAL_PFAS, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(year = factor(year))

if (nrow(nest_pts) == 0) stop("No non-missing lat/long nest-year rows found.")

# Deterministic jitter to reduce overplotting of overlapping nests
set.seed(20260303)
nest_pts <- nest_pts %>%
  mutate(
    long_jit = long + runif(n(), min = -0.12, max = 0.12),
    lat_jit = lat + runif(n(), min = -0.08, max = 0.08)
  )

# Parse a GeoJSON Polygon/MultiPolygon into long/lat rows
coords_to_df <- function(coords, geometry_type, region_name, group_start = 1L, exterior_only = FALSE) {
  ring_to_mat <- function(ring) {
    do.call(
      rbind,
      lapply(ring, function(pt) c(as.numeric(pt[[1]]), as.numeric(pt[[2]])))
    )
  }

  out <- list()
  g <- group_start
  if (geometry_type == "Polygon") {
    rings <- if (exterior_only) coords[1] else coords
    for (ring in rings) {
      mat <- ring_to_mat(ring)
      out[[length(out) + 1L]] <- data.frame(
        long = mat[, 1],
        lat = mat[, 2],
        group = g,
        region = region_name
      )
      g <- g + 1L
    }
  } else if (geometry_type == "MultiPolygon") {
    for (poly in coords) {
      rings <- if (exterior_only) poly[1] else poly
      for (ring in rings) {
        mat <- ring_to_mat(ring)
        out[[length(out) + 1L]] <- data.frame(
          long = mat[, 1],
          lat = mat[, 2],
          group = g,
          region = region_name
        )
        g <- g + 1L
      }
    }
  }
  list(df = bind_rows(out), next_group = g)
}

lines_to_df <- function(coords, geometry_type, name, group_start = 1L) {
  path_to_mat <- function(path) {
    do.call(
      rbind,
      lapply(path, function(pt) c(as.numeric(pt[[1]]), as.numeric(pt[[2]])))
    )
  }

  out <- list()
  g <- group_start
  if (geometry_type == "LineString") {
    mat <- path_to_mat(coords)
    out[[1L]] <- data.frame(long = mat[, 1], lat = mat[, 2], group = g, name = name)
    g <- g + 1L
  } else if (geometry_type == "MultiLineString") {
    for (ln in coords) {
      mat <- path_to_mat(ln)
      out[[length(out) + 1L]] <- data.frame(long = mat[, 1], lat = mat[, 2], group = g, name = name)
      g <- g + 1L
    }
  }
  list(df = bind_rows(out), next_group = g)
}

gj <- fromJSON(geo_file, simplifyVector = FALSE)
features <- gj$features

poly_list <- list()
group_id <- 1L
for (ft in features) {
  nm <- ft$properties$name
  gt <- ft$geometry$type
  cd <- ft$geometry$coordinates
  parsed <- coords_to_df(cd, gt, nm, group_id)
  if (nrow(parsed$df) > 0) poly_list[[length(poly_list) + 1L]] <- parsed$df
  group_id <- parsed$next_group
}

us_poly <- bind_rows(poly_list) %>%
  filter(!region %in% c("Alaska", "Hawaii", "Puerto Rico"))

wi_poly <- us_poly %>% filter(region == "Wisconsin")
if (nrow(wi_poly) == 0) stop("Wisconsin polygon not found in GeoJSON.")

# Wisconsin county boundaries (FIPS prefix 55)
gj_county <- fromJSON(county_geo_file, simplifyVector = FALSE)
county_features <- gj_county$features

wi_county_list <- list()
group_id_county <- 1L
for (ft in county_features) {
  geoid <- ft$properties$GEO_ID
  if (is.null(geoid)) next
  if (!startsWith(geoid, "0500000US55")) next

  gt <- ft$geometry$type
  cd <- ft$geometry$coordinates
  parsed <- coords_to_df(cd, gt, region_name = "WI_county", group_start = group_id_county, exterior_only = TRUE)
  if (nrow(parsed$df) > 0) wi_county_list[[length(wi_county_list) + 1L]] <- parsed$df
  group_id_county <- parsed$next_group
}

wi_counties <- bind_rows(wi_county_list)
if (nrow(wi_counties) == 0) stop("No Wisconsin county polygons found in county GeoJSON.")

# Optional hydrological overlays near Wisconsin
plot_bbox <- list(xmin = -93.8, xmax = -85.8, ymin = 42.0, ymax = 47.8)
north_bbox <- list(xmin = -92.35, xmax = -90.2, ymin = 46.55, ymax = 47.2)


us_inset <- ggplot() +
  geom_polygon(
    data = us_poly,
    aes(x = long, y = lat, group = group),
    fill = "gray92", color = "gray55", linewidth = 0.1
  ) +
  geom_polygon(
    data = wi_poly,
    aes(x = long, y = lat, group = group),
    fill = "black", color = "black", linewidth = 0.25
  ) +
  coord_quickmap(xlim = c(-125, -66), ylim = c(24, 50), expand = FALSE) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Main Wisconsin map
p_main <- ggplot() +
  geom_polygon(data = wi_poly,
               aes(x = long, y = lat, group = group),
               fill = "gray97", color = "black", linewidth = 0.35) +
  geom_polygon(data = wi_counties,
               aes(x = long, y = lat, group = group),
               fill = NA, color = "gray80", linewidth = 0.14) +
  geom_point(data = nest_pts,
             aes(x = long_jit, y = lat_jit, fill = total_pfas),
             shape = 21, color = "white", stroke = 0.3, size = 4.0, alpha = 0.95) +
  annotate(
    "rect",
    xmin = north_bbox$xmin, xmax = north_bbox$xmax,
    ymin = north_bbox$ymin, ymax = north_bbox$ymax,
    fill = NA, color = "black", linewidth = 0.45
  ) +
  annotate(
    "segment",
    x = north_bbox$xmax, y = north_bbox$ymin + 0.05,
    xend = -88.75, yend = 46.95,
    linewidth = 0.45, color = "black",
    arrow = grid::arrow(type = "closed", length = grid::unit(0.12, "in"))
  ) +
  scale_fill_viridis_c(
    name = "Total PFAS",
    option = "C",
    direction = 1
  ) +
  guides(
    fill = guide_colorbar(
      direction = "vertical",
      barheight = grid::unit(70, "pt"),
      barwidth = grid::unit(16, "pt")
    )
  ) +
  coord_quickmap(xlim = c(plot_bbox$xmin, plot_bbox$xmax), ylim = c(plot_bbox$ymin, plot_bbox$ymax), expand = FALSE) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold"),
    legend.position = c(0.08, 0.10),
    legend.text = element_text(size= 12),
    legend.title = element_text(size = 12),
    legend.justification = c(0, 0),
    legend.direction = "vertical",
    #legend.background = element_rect(fill = scales::alpha("white", 0.90), color = "gray70", linewidth = 0.2),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(8, 8, 8, 8)
  )

print(p_main)

# Zoomed inset for northern Wisconsin nests (Ashland/Bayfield area)
p_north <- ggplot() +
  geom_polygon(data = wi_poly,
               aes(x = long, y = lat, group = group),
               fill = "gray97", color = "black", linewidth = 0.3) +
  geom_polygon(data = wi_counties,
               aes(x = long, y = lat, group = group),
               fill = NA, color = "gray80", linewidth = 0.12) +
  geom_point(data = nest_pts,
             aes(x = long_jit, y = lat_jit, fill = total_pfas),
             shape = 21, color = "white", stroke = 0.25, size = 3.3, alpha = 0.95) +
  scale_fill_viridis_c(option = "C", direction = 1, guide = "none") +
  coord_quickmap(
    xlim = c(north_bbox$xmin, north_bbox$xmax),
    ylim = c(north_bbox$ymin, north_bbox$ymax),
    expand = FALSE
  ) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.8),
    plot.background = element_rect(fill = "white", color = NA)
  )

p <- ggdraw() +
  draw_plot(p_main, x = 0.00, y = 0.00, width = 0.68, height = 1.00) +
  draw_plot(p_north, x = 0.48, y = 0.62, width = 0.38, height = 0.37) +
  draw_plot(us_inset, x = 0.6, y = 0.22, width = 0.21, height = 0.18) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

print(p)

ggsave(
  filename = file.path(out_dir, "WI_nest_map_total_pfas_points_jitter_zoom8.png"),
  plot = p, width = 9.5, height = 7.2, dpi = 600, bg = "white"
)

ggsave(
  filename = file.path(out_dir, "WI_nest_map_total_pfas_points_jitter_zoom8.tiff"),
  plot = p, width = 9.5, height = 7.2, dpi = 600, compression = "lzw", bg = "white"
)

cat("Done: Wisconsin nest map with US inset\n")
