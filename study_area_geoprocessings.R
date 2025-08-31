# ======================================================================
# Study Area Maps & Soil Landscapes — Colombia (R)
# ----------------------------------------------------------------------
# Purpose
#   Prepare soil‑landscape polygons, clip them to study areas, compute
#   per‑department landscape areas, and produce multi‑panel maps:
#     • 9 departmental soil‑landscape maps with stations & contours
#     • Microlocation (study area extent) and macrolocation (regional context)
#   Finally, estimate urban‑area composition by landscape type.
#
# Key assumptions & notes
#   - All geometries are transformed to the CRS of `rect_study_areas`.
#   - Validity: optional fix with st_make_valid() (guarded by a check).
#   - Area calculations use EPSG:9377 (MAGNA‑SIRGAS / Origen Nacional) in m².
#   - Landscape labels are harmonized (e.g., long “Tierras misceláneas...” label).
#   - Department labels corrected via `study_areas_corrected` for plotting.
#   - Station points are read from a CSV with columns X/Y in WGS84.
#   - The color palette is fixed and mapped to unique `PAISAJE` values.
# ======================================================================

require(pacman)
p_load(sf,
       geodata,
       tidyverse,
       ggspatial,
       ggrepel,
       ggforce,
       patchwork,
       cowplot,
       broom,
       gt)

# ----------------------------------------------------------------------
# Read inputs: soils and study areas (polygons)
# ----------------------------------------------------------------------
rect_soils <- read_sf("study_area_sig/clasificacion_suelo.gpkg")
rect_study_areas <- read_sf("data/boundaries/boundaries_study_area_col_v2.gpkg")

# Harmonize CRS to study areas
rect_soils <- rect_soils %>% st_transform(st_crs(rect_study_areas))

# Geometry validity check (run fix only if invalids are present)
if (FALSE %in% unique(st_is_valid(rect_soils))) {
  rect_soils <- st_make_valid(rect_soils)
}

# Clip soils to study areas extent (intersection)
rect_soils_cropped <- st_intersection(rect_soils, rect_study_areas)

# Persist cropped soils as GeoPackage (overwrite)
st_write(
  obj = rect_soils_cropped,
  dsn = "study_area_sig/clasificacion_suelo_recortado.gpkg",
  layer = "rect_soils",
  driver = "GPKG",
  delete_dsn = TRUE
)

# ======================================================================
# MAIN MAPS — departmental panels with soils, contours, and stations
# ======================================================================
# Reload study areas and prepare labels + bounding boxes for panel extents
rect_study_areas <- read_sf("data/boundaries/boundaries_study_area_col_v2.gpkg")
rect_study_areas <- rect_study_areas %>%
  select(-grupo) %>%
  mutate(
    study_areas_corrected = c(
      "ANTIOQUIA",
      "CALDAS",
      "BOGOTÁ",
      "BOYACÁ",
      "SANTANDER",
      "NORTE DE SANTANDER",
      "CESAR",
      "BOLÍVAR",
      "MAGDALENA"
    )
  )

# Precompute per‑feature bounding boxes (for coord_sf limits)
rect_study_areas <- rect_study_areas %>%
  mutate(bounding_box = map(geom, st_bbox)) %>%
  mutate(
    xmin = map_dbl(bounding_box, ~ .x[["xmin"]]),
    xmax = map_dbl(bounding_box, ~ .x[["xmax"]]),
    ymin = map_dbl(bounding_box, ~ .x[["ymin"]]),
    ymax = map_dbl(bounding_box, ~ .x[["ymax"]])
  ) %>%
  select(-bounding_box)

# ----------------------------------------------------------------------
# Base layers: world/Colombia boundaries and topographic contours
# ----------------------------------------------------------------------
# If needed once: rect_world <- geodata::world(...)
rect_world <- read_rds("../E-SIG/SHP-Limites administrativos COL - GADM/gadm/gadm36_adm0_r5_pk.rds")
rect_world <- rect_world %>%
  terra::unwrap() %>% 
  st_as_sf() %>%
  st_transform(crs = st_crs(rect_study_areas))

# World centroids for labeling (filtered to the map window)
rect_world_centroids <- rect_world %>%
  st_centroid(of_largest_polygon = TRUE) %>%
  cbind(st_coordinates(.)) %>%
  filter(X >= -85, X <= -60, Y >= -5, Y <= 15)

# Colombia admin levels (country, dept, municipality)
rect_colombia_lv0 <- read_sf("../E-SIG/SHP-Limites administrativos COL - GADM/gadm41_COL_0.shp")
rect_colombia_lv0 <- rect_colombia_lv0 %>% st_transform(crs = st_crs(rect_study_areas))

rect_colombia_lv1 <- read_sf("../E-SIG/SHP-Limites administrativos COL - GADM/gadm41_COL_1.shp")
rect_colombia_lv1 <- rect_colombia_lv1 %>% st_transform(crs = st_crs(rect_study_areas))

rect_colombia_lv2 <- read_sf("../E-SIG/SHP-Limites administrativos COL - GADM/gadm41_COL_2.shp")
rect_colombia_lv2 <- rect_colombia_lv2 %>% st_transform(crs = st_crs(rect_study_areas))

# Topographic contours (Cartografía IGAC)
line_topography <- read_sf("../E-SIG/Cartografía Colombia Agustin Codazzi/Curva_Nivel.shp")
line_topography <- line_topography %>% st_transform(crs = st_crs(rect_study_areas))

# ----------------------------------------------------------------------
# Read cropped soils; harmonize labels; build dissolved landscapes
# ----------------------------------------------------------------------
rect_soils_cropped <- read_sf("study_area_sig/clasificacion_suelo_recortado.gpkg")
rect_soils_cropped <- rect_soils_cropped %>% st_transform(crs = st_crs(rect_study_areas))
rect_soils_cropped <- rect_soils_cropped %>%
  mutate(
    PAISAJE = case_when(
      PAISAJE == "Tierras misceláneas, con pendientes mayores del 75%, relieve muy escarpado" ~ "Tierras misceláneas\npendiente mayor al 75%",
      TRUE ~ PAISAJE
    )
  )
# Dissolve by landscape type
data_placeholder <- NULL
rect_soils_landscape <- rect_soils_cropped %>%
  group_by(PAISAJE) %>%
  summarise(geom = st_union(geom), .groups = "drop")

# ----------------------------------------------------------------------
# Areas per department × landscape (top 3 by area each dept)
# ----------------------------------------------------------------------
rect_soils_per_department <- rect_soils_cropped %>%
  group_by(department, PAISAJE) %>%
  summarise(geom = st_union(geom), .groups = "drop")

rect_soils_per_department %>%
  st_transform(crs = 9377) %>%
  # area in km2
  mutate(area = st_area(geom) / 10 ^ 6) %>%
  st_drop_geometry() %>%
  group_by(department) %>%
  slice_max(order_by = area,
            n = 3,
            with_ties = FALSE) %>%
  inner_join((
    rect_study_areas %>%
      select(department, geom) %>%
      st_transform(crs = 9377) %>%
      mutate(area_dep = st_area(geom) / 10 ^ 6) %>%
      st_drop_geometry()
  )) %>%
  arrange(department, - area) %>% 
  mutate(prop = 100 * area / area_dep) %>%
  mutate(across(c(contains("area"), prop), ~ as.numeric(.x))) %>%
  write.table(
    "data/processed/statistic_summaries/study_area_soil_landscapes.txt",
    row.names = FALSE,
    sep = ";",
    dec = ","
  )

# Urban centers (MGN 2024)
rect_urban_centers <- read_sf("../E-SIG/MGN2024_URB_ZONA_URB_ZONA_URBANA.shp")
rect_urban_centers <- read_sf("../E-SIG/MGN2024_URB_ZONA_URBANA/MGN_URB_ZONA_URBANA.shp")
rect_urban_centers <- rect_urban_centers %>% st_transform(crs = st_crs(rect_study_areas))

# Monitoring stations (point layer from coordinates table)
point_stations <- read.table(
  "data/processed/SISAIRE/data_sisaire_coordinates.txt",
  sep = ";",
  dec = ",",
  header = TRUE
)
point_stations <- point_stations %>%
  st_as_sf(coords = c("X", "Y"), crs = 4326) %>%
  st_transform(crs = st_crs(rect_study_areas))

# ----------------------------------------------------------------------
# Plot tools: palette & panel titles
# ----------------------------------------------------------------------
palette <- tibble(
  paisajes = unique(rect_soils_landscape$PAISAJE),
  color = c(
    "#FF6F61",
    "#6B5B95",
    "#9B1B30",
    "#77212E",
    "#F5D6C6",
    "#FA9A85",
    "#00539C",
    "#CE5B78",
    "#935529",
    "#009B77",
    "#2A4B7C",
    "#577284",
    "#F96714",
    "#264E36",
    "#2A293E",
    "#797B3A",
    "#DD4132",
    "#C62168",
    "#5A3E36"
  )
)

data_plot_titles <- tibble(
  department = unique(rect_study_areas$department),
  plot_title = c(
    "ANTIOQUIA",
    "CALDAS",
    "BOGOTÁ",
    "BOYACÁ",
    "SANTANDER",
    "NORTE DE SANTANDER",
    "CESAR",
    "BOLÍVAR",
    "MAGDALENA"
  )
)

# ----------------------------------------------------------------------
# Function: departmental soil‑landscape plot
# ----------------------------------------------------------------------
generate_plot_paisajes <- function(department_selected) {
  # Panel extent from precomputed bbox columns
  department_limits <- rect_study_areas %>%
    filter(department == department_selected) %>%
    st_drop_geometry() %>%
    select(contains(c("x", "y")))
  
  # Title for the panel
  data_plot_title <- data_plot_titles %>%
    filter(department == department_selected) %>%
    pull(plot_title)
  
  # Construct layered map for one department
  graph_study_area_soils <- ggplot() +
    geom_sf(data = rect_colombia_lv0, linewidth = 1) +
    
    geom_sf(data = rect_study_areas,
            fill = "transparent",
            linewidth = .5) +
    
    # Landcape type
    geom_sf(data = rect_soils_landscape,
            aes(fill = PAISAJE),
            color = "gray10") +
    
    # Topography
    geom_sf(
      data = line_topography,
      aes(linetype = "solid"),
      color = "black",
      alpha = .5
    ) +
    
    # Stations
    geom_sf(
      data = point_stations,
      aes(shape = ""),
      color = "black",
      fill = "purple",
      stroke = .8,
      size = 3
    ) +
    
    scale_shape_manual(name = "AQ Monitoring\nStation", values = 21) +
    scale_linetype_manual(values = "solid", labels = "") +
    scale_fill_manual(
      breaks = palette$paisajes,
      values = palette$color,
      labels = landscape_labels_en <- c(
        "High Plateau",
        "Structural\nHigh Plateau",
        "Hill",
        "Water Body",
        "Mining Pits",
        "Low Hill",
        "Eroded\nMiscellaneous",
        "Rocky\nMiscellaneous",
        "Mountain",
        "Erosional\nStructural Mountains",
        "Swamps",
        "Piedmont",
        "Plain",
        "Alluvial Plain",
        "Miscellaneous Lands\nSlopes over 75%",
        "Valley",
        "Alluvial Valley",
        "Mountain Slope",
        "Urban Area"
      )
    ) +
    
    labs(title = data_plot_title,
         fill = "Landscape",
         linetype = "500m Elevation contours") +
    
    coord_sf(
      xlim = c(department_limits$xmin, department_limits$xmax),
      ylim = c(department_limits$ymin , department_limits$ymax),
      expand = 0
    ) +
    
    ggspatial::annotation_scale(
      location = "br",
      bar_cols = c("black", "white"),
      width_hint = .25,
      text_cex = 1,
      text_face = "bold",
      text_family = "serif"
    ) +
    
    theme_bw() +
    
    theme(
      text = element_text(family = "serif"),
      axis.text.x = element_text(size = 8,
                                 angle = 40,
                                 vjust = 0.5),
      panel.background = element_rect(fill = "#74beea"),
      panel.grid = element_line(colour = "#e0e0e0"),
      plot.title = element_text(hjust = .5, face = "bold"),
      legend.text = element_text(hjust = .5, size = 10),
      legend.position = "bottom",
      legend.title.position = "top",
      legend.key = element_rect(fill = "transparent", color = NA),
      aspect.ratio = 1
    ) +
    guides(linetype = guide_legend(override.aes = list(linewidth = 1.5)),
           shape = guide_legend(legend.text = NA))
  
  return(graph_study_area_soils)
}

# Build all nine departmental panels
graph_soils_antioquia <- generate_plot_paisajes("ANTIOQUIA")
graph_soils_bogota <- generate_plot_paisajes("BOGOTA_DC")
graph_soils_bolivar <- generate_plot_paisajes("BOLIVAR")
graph_soils_boyaca <- generate_plot_paisajes("BOYACA")
graph_soils_caldas <- generate_plot_paisajes("CALDAS")
graph_soils_cesar <- generate_plot_paisajes("CESAR")
graph_soils_magdalena <- generate_plot_paisajes("MAGDALENA")
graph_soils_norte_de_santander <- generate_plot_paisajes("NORTE_DE_SANTANDER")
graph_soils_santander <- generate_plot_paisajes("SANTANDER")

# Arrange panels with shared legend via patchwork
graphs_soils <- (
  (
    graph_soils_antioquia | graph_soils_bogota | graph_soils_bolivar
  ) /
    (graph_soils_boyaca |
       graph_soils_caldas | graph_soils_cesar) /
    (
      graph_soils_magdalena |
        graph_soils_norte_de_santander | graph_soils_santander
    )
) + plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.title.position = "top",
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.text = element_text(hjust = .5),
    plot.margin = margin(5, 10, 5, 10)
  )

# Save multi‑panel figure in three formats
ggsave(
  "graphs/hourly/potential/study_area_soils.pdf",
  graphs_soils,
  width = 14,
  height = 12,
  dpi = 500
)
ggsave(
  "graphs/hourly/potential/study_area_soils.tiff",
  graphs_soils,
  width = 14,
  height = 12,
  dpi = 500
)
ggsave(
  "graphs/hourly/potential/study_area_soils.jpg",
  graphs_soils,
  width = 14,
  height = 12,
  dpi = 500
)

# ======================================================================
# Microlocation map (zoom into study areas with stations & urban areas)
# ======================================================================
rect_study_areas_centroids <- rect_study_areas %>%
  st_centroid(of_largest_polygon = TRUE) %>%
  cbind(st_coordinates(.)) %>%
  mutate(
    Y_corrected = case_when(
      department == "ANTIOQUIA" ~ Y + .6,
      department == "CALDAS" ~ Y + .3,
      department == "BOGOTA_DC" ~ Y + .6,
      department == "BOYACA" ~ Y + .4,
      department == "SANTANDER" ~ Y + .4,
      department == "NORTE_DE_SANTANDER" ~ Y + .4,
      department == "CESAR" ~ Y + .7,
      department == "BOLIVAR" ~ Y + .7,
      TRUE ~ Y + .5
    )
  )

# Microlocation composition
graph_study_area_microlocation <- ggplot() +
  geom_sf(data = rect_colombia_lv0, linewidth = 1) +
  
  geom_sf(data = rect_study_areas,
          aes(color = ""),
          fill = "transparent",
          linewidth = 2) +
  
  geom_sf(data = rect_urban_centers, linewidth = .5, aes(fill = "gray30")) +
  
  geom_sf(
    data = point_stations,
    aes(shape = ""),
    color = "black",
    fill = "purple",
    size = 3
  ) +
  
  geom_text(
    data = rect_study_areas_centroids,
    aes(x = X, y = Y_corrected, label = study_areas_corrected),
    family = "serif",
    fontface = "bold.italic"
  ) +
  scale_shape_manual(name = "AQ Monitoring\nstation", values = 21) +
  scale_fill_manual(name = "Urban areas",
                    values = "gray30",
                    labels = "") +
  scale_color_manual(name = "Study areas", values = "#d6604d") +
  coord_sf(xlim = c(-76.1, -70), ylim = c(4.38 , 11.35)) +
  ggspatial::annotation_north_arrow(
    location = "tl",
    which_north = "true",
    style = north_arrow_nautical(text_family = "serif", text_face = "bold"),
    pad_x = unit(0.2, "cm"),
    pad_y = unit(0.2, "cm"),
    width = unit(3, "cm"),
    height = unit(3, "cm")
  ) +
  ggspatial::annotation_scale(
    location = "br",
    bar_cols = c("black", "white"),
    width_hint = .5,
    text_cex = 1.2,
    text_face = "bold",
    text_family = "serif"
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "serif"),
    axis.text = element_text(size = 15),
    axis.title = element_blank(),
    panel.background = element_rect(fill = "#74beea"),
    panel.grid = element_line(color = "#e0e0e0"),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.title.position = "left",
    legend.text = element_text(hjust = .5),
    legend.key.width = unit(1, "cm"),
    legend.key = element_rect(fill = "transparent", color = NA),
    plot.title = element_text(hjust = .5),
    plot.margin = margin(2, 5, 0, 0),
    panel.spacing = unit(0, "pt")
  )

# Save microlocation
ggsave(
  "graphs/study_area_microlocation.pdf",
  graph_study_area_microlocation,
  width = 6,
  height = 8
)

# ======================================================================
# Macrolocation map (regional context, bounding rectangle of study area)
# ======================================================================
graph_study_area_macrolocation <- ggplot() +
  
  geom_sf(data = rect_world, linewidth = 1) +
  
  geom_sf(data = rect_colombia_lv0, linewidth = 1) +
  
  geom_rect(
    aes(
      xmin = -76.05,
      xmax = -72.22,
      ymin = 4.38,
      ymax = 11.35
    ),
    fill = "transparent",
    color = "#d6604d",
    linewidth = 1.5
  ) +
  
  geom_text_repel(
    aes(x = -82, y = 2, label = "PACIFIC OCEAN"),
    angle = 60,
    family = "serif",
    fontface = "bold.italic"
  ) +
  
  geom_text_repel(
    aes(x = -76, y = 12.5, label = "CARIBBEAN SEA"),
    angle = 15,
    family = "serif",
    fontface = "bold.italic"
  ) +
  
  scale_fill_manual(name = "Urban areas",
                    values = "gray30",
                    labels = "") +
  scale_y_continuous(breaks = seq(-4, 12, 4)) +
  coord_sf(xlim = c(-85, -60),
           ylim = c(-5 , 15),
           expand = 0) +
  theme_bw() +
  theme(
    text = element_text(family = "serif"),
    axis.text = element_text(color = "black",
                             size = 8),
    axis.title = element_blank(),
    panel.background = element_rect(fill = "#74beea", color = "black"),
    legend.position = "bottom",
    legend.title.position = "top",
    panel.spacing = unit(0, "pt"),
    panel.border = element_rect(fill = "transparent", linewidth = 1),
    panel.grid = element_line(color = "#e0e0e0"),
    plot.title = element_text(hjust = .5),
    plot.margin = margin(0, 5, 0, 0),
    plot.background = element_blank()
  )

# Save macrolocation
ggsave(
  "graphs/study_area_macrolocation.pdf",
  graph_study_area_macrolocation,
  width = 10,
  height = 10
)

# ======================================================================
# Combine micro + macro panels into one canvas
# ======================================================================
 micro_macro_stack <- ggdraw() +
  draw_plot(graph_study_area_microlocation) +
  draw_plot(
    graph_study_area_macrolocation,
    x = 0.58,
    y = 0.55,
    width = 0.4,
    height = 0.4
  )

ggsave(
  "graphs/hourly/potential/study_area.pdf",
  micro_macro_stack,
  width = 10,
  height = 12,
  dpi = 500
)
ggsave(
  "graphs/hourly/potential/study_area.tiff",
  micro_macro_stack,
  width = 10,
  height = 12,
  dpi = 500
)
ggsave(
  "graphs/hourly/potential/study_area.jpg",
  micro_macro_stack,
  width = 10,
  height = 12,
  dpi = 500
)

# ======================================================================
# Urban‑area composition by landscape
# ======================================================================
# Intersect urban polygons with study areas and with departmental landscapes,
# then compute areas and proportions by landscape within the urban footprint.
rect_urban_centers_cropped <- rect_urban_centers %>% 
  select(geometry) %>% 
  st_intersection(rect_study_areas)

sf_use_s2(FALSE)
rect_urban_landscapes <- st_intersection(rect_soils_per_department, rect_urban_centers_cropped)
sf_use_s2(TRUE)

data_landscapes_area <- rect_urban_landscapes %>% group_by(department, PAISAJE) %>% 
  # Plain crs to calculate the areas (Magna SIRGAS - Origen Nacional, Colombia)
  st_transform(crs = 9377) %>% 
  reframe(area = sum(st_area(geom)) / 10^6) %>% 
  arrange(department, - area) %>% 
  group_by(department) %>% 
  mutate(total_urb_area = sum(area),
         prop_urb_area = 100 * area / total_urb_area)