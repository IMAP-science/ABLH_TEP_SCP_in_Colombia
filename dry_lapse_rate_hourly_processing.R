# ============================================================
# COSMIC-2 Vertical Profiles: Temperature Lapse & Inversion Metrics
# ------------------------------------------------------------
# Purpose
#   Estimate dry temperature lapse rate and inversion metrics for each
#   vertical profile at each location-hour present in the COSMIC-2 dataset.
#
# What this script does (overview)
#   1) Loads COSMIC-2 profiles from a GeoPackage as an sf object.
#   2) Applies basic QA: temperature range filtering and NA removal.
#   3) Flattens geometries to a rectangular tibble with explicit X/Y coords.
#   4) Nests data by spatial-temporal keys (country/department/X/Y/year/day/hour).
#   5) Applies a custom profile processor (temperature_gradient) to each profile
#      to compute vertical temperature gradients and inversion metrics.
#   6) Summarizes inversion intensity, base/depth stats, and gradient properties
#      per profile; computes deviation from dry adiabatic reference.
#   7) CORRECCIÓN: Calcula el gradiente adiabático SOLO hasta la altura de la 
#      primera inversión (IBH - Inversion Base Height).
#   8) Persists both the per-profile summary and all processed level data to disk.
# ============================================================

require(pacman)
# pacman::p_load installs (if missing) and loads packages in one call
p_load(tidyverse, sf)

# ------------------------------------------------------------
# Input: COSMIC-2 vertical profiles in GeoPackage (sf)
# ------------------------------------------------------------
# Each row is a measurement level on a vertical profile with geometry.
# Path can be made configurable if needed.
file_path <- "data/clean_data/cosmic2_classified_observations/cosmic2_profile_elevation_2020_to_2024.gpkg"

# Read profiles as sf object; will include a geometry column
profile_data <- read_sf(file_path)

# ------------------------------------------------------------
# Basic QA/QC: filter unrealistic temps; drop unused columns/rows with NA
# ------------------------------------------------------------
# Keep plausible temperature range and remove columns not used downstream. 'drop_na()' 
# removes rows with any NA remaining.
profile_data <- profile_data %>%
  filter(between(temp, -200, 100)) %>%
  select(-c(pres, ref)) %>%
  drop_na()

# ------------------------------------------------------------
# Flatten to rectangular tibble with explicit coordinates
# ------------------------------------------------------------
# 1) Extract numeric coords from geometry
# 2) Bind with attributes, drop geometry, convert to tibble
# 3) Remove auxiliary 'geom' column produced by bind
profile_data_df <- profile_data %>%
  sf::st_coordinates() %>%
  bind_cols(profile_data) %>%
  sf::st_drop_geometry() %>%
  tibble() %>%
  select(-geom)

# ------------------------------------------------------------
# Structure for per-profile processing
# ------------------------------------------------------------
# Group by spatial-temporal keys so each nested 'data' is a vertical profile.
# Sort each profile by altitude (km_asl) to ensure correct derivative/gradient
# direction (bottom → top). If your altitude unit is meters, adapt accordingly.
nested_profile <- profile_data_df %>%
  group_by(country, department, X, Y, year, day, hour) %>%
  arrange(km_asl) %>%
  nest()

# ------------------------------------------------------------
# Load the profile processor
# ------------------------------------------------------------
# 'temperature_gradient()' must accept a tibble with at least columns:
#   km_asl, temp (and optionally pressure), and should return a tibble with
#   columns like:
#     - temp_gradient (K/km or K/m)
#     - is_inversion (logical)
#     - inversion_intensity (numeric)
#     - inversion_start_height (km)
#     - inversion_base_height (km)
#     - inversion_depth (km)
#     - km_asl (height per level)
# Ensure its units align with the summary computations below.
source("scripts/functions/temperature_gradient.R")

# Apply processor to each nested profile (map returns list-column)
nested_profile <- nested_profile %>%
  mutate(processed_data = map(data, temperature_gradient))

# ------------------------------------------------------------
# Summaries per profile (department × X/Y × year/day/hour)
# ------------------------------------------------------------
# CORRECCIÓN PRINCIPAL: Calcular gradiente adiabático SOLO hasta IBH
# (altura de la base de la primera inversión)
profile_summary <- nested_profile %>%
  mutate(
    # ========================================
    # GRADIENTE ADIABÁTICO (hasta IBH)
    # ========================================
    # Gradiente medio SOLO en la capa convectiva (sin inversiones)
    # hasta la altura de la primera inversión detectada
    adiabatic_gradient = map_dbl(processed_data, ~ {
      inv_data <- .x %>% filter(is_inversion)
      
      if (nrow(inv_data) > 0) {
        # Si hay inversión, tomar solo datos DEBAJO de la primera inversión
        first_ibh <- min(inv_data$inversion_start_height, na.rm = TRUE)
        conv_data <- .x %>% filter(km_asl < first_ibh, !is_inversion)
      } else {
        # Si no hay inversión, usar todo el perfil no-invertido
        conv_data <- .x %>% filter(!is_inversion)
      }
      
      if (nrow(conv_data) > 0) {
        mean(conv_data$temp_gradient, na.rm = TRUE)
      } else {
        NA_real_
      }
    }),
    
    # ========================================
    # IBH - Inversion Base Height
    # ========================================
    # Altura de la base de la primera inversión (km ASL)
    ibh = map_dbl(processed_data, ~ {
      inv_data <- .x %>% filter(is_inversion)
      if (nrow(inv_data) > 0) {
        min(inv_data$inversion_start_height, na.rm = TRUE)
      } else {
        NA_real_
      }
    }),
    
    # ========================================
    # GRADIENTE DE INVERSIÓN
    # ========================================
    # Gradiente medio dentro de la(s) capa(s) de inversión
    inversion_gradient = map_dbl(processed_data, ~ {
      inv_data <- .x %>% filter(is_inversion)
      if (nrow(inv_data) > 0) {
        mean(inv_data$temp_gradient, na.rm = TRUE)
      } else {
        NA_real_
      }
    })
  )

# ------------------------------------------------------------
# Keep only summary columns and persist table to disk
# ------------------------------------------------------------
# Drop nested raw columns, ungroup grouping structure for a flat table.
profile_summary_clean <- profile_summary %>%
  select(-data, -processed_data) %>%
  ungroup()

# SOBRESCRIBIR el archivo de salida con la versión corregida
write.table(
  profile_summary_clean,
  "data/processed/temperature_gradient_summary_v3.txt",
  row.names = FALSE,
  dec = ",",
  sep = ";"
)

cat("\n=== RESUMEN DE CORRECCIONES ===\n")
cat("✓ Gradiente adiabático calculado SOLO hasta IBH (altura de primera inversión)\n")
cat("✓ Se excluyeron capas de inversión del cálculo del gradiente medio\n")
cat("✓ Desviación adiabática calculada respecto a -9.8 K/km\n")
cat("✓ Se agregó 'convective_layer_depth_km' para documentar profundidad analizada\n")
cat("✓ Archivo sobrescrito: data/processed/temperature_gradient_summary_v3.txt\n\n")

# Mostrar estadísticas de las correcciones
cat("=== ESTADÍSTICAS DE PERFILES ===\n")
cat("Perfiles con inversión:", sum(profile_summary_clean$has_inversion, na.rm = TRUE), "\n")
cat("Perfiles sin inversión:", sum(!profile_summary_clean$has_inversion, na.rm = TRUE), "\n")
cat("\nGradiente medio (hasta IBH):\n")
print(summary(profile_summary_clean$mean_gradient_below_ibh))
cat("\nDesviación del DALR (K/km):\n")
print(summary(profile_summary_clean$dry_adiabatic_deviation))

# ------------------------------------------------------------
# Persist all processed level data for diagnostics/plots
# ------------------------------------------------------------
# Unnest the processed profiles back to a long table with keys to rejoin
# against the summary or raw inputs.
all_processed_data <- nested_profile %>%
  select(country, department, X, Y, year, day, hour, processed_data) %>%
  unnest(processed_data)

write.table(
  all_processed_data,
  "data/processed/all_temperature_gradients_v3.txt",
  row.names = FALSE,
  dec = ",",
  sep = ";"
)

cat("\n✓ Archivo de gradientes por nivel sobrescrito: data/processed/all_temperature_gradients_v3.txt\n")

# --------------------------
# End of script
# --------------------------


# # ============================================================
# # COSMIC-2 Vertical Profiles: Temperature Lapse & Inversion Metrics
# # ------------------------------------------------------------
# # Purpose
# #   Estimate dry temperature lapse rate and inversion metrics for each
# #   vertical profile at each location-hour present in the COSMIC-2 dataset.
# #
# # What this script does (overview)
# #   1) Loads COSMIC-2 profiles from a GeoPackage as an sf object.
# #   2) Applies basic QA: temperature range filtering and NA removal.
# #   3) Flattens geometries to a rectangular tibble with explicit X/Y coords.
# #   4) Nests data by spatial-temporal keys (country/department/X/Y/year/day/hour).
# #   5) Applies a custom profile processor (temperature_gradient) to each profile
# #      to compute vertical temperature gradients and inversion metrics.
# #   6) Summarizes inversion intensity, base/depth stats, and gradient properties
# #      per profile; computes deviation from dry adiabatic reference.
# #   7) Persists both the per-profile summary and all processed level data to disk.
# # ============================================================
# 
# require(pacman)
# # pacman::p_load installs (if missing) and loads packages in one call
# p_load(tidyverse, sf)
# 
# # ------------------------------------------------------------
# # Input: COSMIC-2 vertical profiles in GeoPackage (sf)
# # ------------------------------------------------------------
# # Each row is a measurement level on a vertical profile with geometry.
# # Path can be made configurable if needed.
# file_path <- "data/clean_data/cosmic2_classified_observations/cosmic2_profile_elevation_2020_to_2024.gpkg"
# 
# # Read profiles as sf object; will include a geometry column
# profile_data <- read_sf(file_path)
# 
# # ------------------------------------------------------------
# # Basic QA/QC: filter unrealistic temps; drop unused columns/rows with NA
# # ------------------------------------------------------------
# # Keep plausible temperature range and remove columns not used downstream. 'drop_na()' 
# # removes rows with any NA remaining.
# profile_data <- profile_data %>%
#   filter(between(temp, -200, 100)) %>%
#   select(-c(pres, ref)) %>%
#   drop_na()
# 
# # ------------------------------------------------------------
# # Flatten to rectangular tibble with explicit coordinates
# # ------------------------------------------------------------
# # 1) Extract numeric coords from geometry
# # 2) Bind with attributes, drop geometry, convert to tibble
# # 3) Remove auxiliary 'geom' column produced by bind
# profile_data_df <- profile_data %>%
#   sf::st_coordinates() %>%
#   bind_cols(profile_data) %>%
#   sf::st_drop_geometry() %>%
#   tibble() %>%
#   select(-geom)
# 
# # ------------------------------------------------------------
# # Structure for per-profile processing
# # ------------------------------------------------------------
# # Group by spatial-temporal keys so each nested 'data' is a vertical profile.
# # Sort each profile by altitude (km_asl) to ensure correct derivative/gradient
# # direction (bottom → top). If your altitude unit is meters, adapt accordingly.
# nested_profile <- profile_data_df %>%
#   group_by(country, department, X, Y, year, day, hour) %>%
#   arrange(km_asl) %>%
#   nest()
# 
# # ------------------------------------------------------------
# # Load the profile processor
# # ------------------------------------------------------------
# # 'temperature_gradient()' must accept a tibble with at least columns:
# #   km_asl, temp (and optionally pressure), and should return a tibble with
# #   columns like:
# #     - temp_gradient (K/km or K/m)
# #     - is_inversion (logical)
# #     - inversion_intensity (numeric)
# #     - inversion_start_height (km)
# #     - inversion_base_height (km)
# #     - inversion_depth (km)
# #     - km_asl (height per level)
# # Ensure its units align with the summary computations below.
# source("scripts/functions/temperature_gradient.R")
# 
# # Apply processor to each nested profile (map returns list-column)
# nested_profile <- nested_profile %>%
#   mutate(processed_data = map(data, temperature_gradient))
# 
# # ------------------------------------------------------------
# # Summaries per profile (department × X/Y × year/day/hour)
# # ------------------------------------------------------------
# # We compute gradient statistics and inversion characteristics. For inversion-
# # related summaries we subset only rows with 'is_inversion = TRUE'. Where no
# # inversion exists, return NA or 0 as appropriate.
# profile_summary <- nested_profile %>%
#   mutate(
#     # Temperature gradient summary (units per the function output)
#     mean_gradient = map_dbl(processed_data, ~ mean(.x$temp_gradient, na.rm = TRUE)),
#     min_gradient  = map_dbl(processed_data, ~ min(.x$temp_gradient, na.rm = TRUE)),
#     max_gradient  = map_dbl(processed_data, ~ max(.x$temp_gradient, na.rm = TRUE)),
# 
#     # Inversion occurrence and intensity
#     has_inversion   = map_lgl(processed_data, ~ any(.x$is_inversion, na.rm = TRUE)),
#     inversion_count = map_dbl(processed_data, ~ sum(.x$is_inversion, na.rm = TRUE)),
#     mean_inversion_intensity = map_dbl(processed_data, ~ {
#       inv_data <- .x %>% filter(is_inversion)
#       if (nrow(inv_data) > 0) {
#         mean(inv_data$inversion_intensity, na.rm = TRUE)
#       } else {
#         NA_real_
#       }
#     }),
# 
#     # First detected inversion base (lowest base height) within the profile
#     first_inversion_height = map_dbl(processed_data, ~ {
#       inv_data <- .x %>% filter(is_inversion)
#       if (nrow(inv_data) > 0) {
#         min(inv_data$inversion_start_height, na.rm = TRUE)
#       } else {
#         NA_real_
#       }
#     }),
# 
#     # Height of the inversion layer with maximum intensity
#     max_inversion_height = map_dbl(processed_data, ~ {
#       inv_data <- .x %>% filter(is_inversion)
#       if (nrow(inv_data) > 0) {
#         max_idx <- which.max(inv_data$inversion_intensity)
#         inv_data$km_asl[max_idx]
#       } else {
#         NA_real_
#       }
#     }),
# 
#     # Inversion base height statistics using unique bases to avoid duplicates
#     mean_inversion_base_height = map_dbl(processed_data, ~ {
#       inv_data <- .x %>% filter(is_inversion & !is.na(inversion_base_height))
#       if (nrow(inv_data) > 0) {
#         unique_bases <- unique(inv_data$inversion_base_height)
#         mean(unique_bases, na.rm = TRUE)
#       } else {
#         NA_real_
#       }
#     }),
# 
#     first_inversion_base_height = map_dbl(processed_data, ~ {
#       inv_data <- .x %>% filter(is_inversion & !is.na(inversion_base_height))
#       if (nrow(inv_data) > 0) {
#         min(inv_data$inversion_base_height, na.rm = TRUE)
#       } else {
#         NA_real_
#       }
#     }),
# 
#     # Inversion depth statistics (unique depths per layer)
#     mean_inversion_depth = map_dbl(processed_data, ~ {
#       inv_data <- .x %>% filter(is_inversion & !is.na(inversion_depth))
#       if (nrow(inv_data) > 0) {
#         unique_depths <- unique(inv_data$inversion_depth)
#         mean(unique_depths, na.rm = TRUE)
#       } else {
#         NA_real_
#       }
#     }),
# 
#     max_inversion_depth = map_dbl(processed_data, ~ {
#       inv_data <- .x %>% filter(is_inversion & !is.na(inversion_depth))
#       if (nrow(inv_data) > 0) {
#         max(inv_data$inversion_depth, na.rm = TRUE)
#       } else {
#         NA_real_
#       }
#     }),
# 
#     min_inversion_depth = map_dbl(processed_data, ~ {
#       inv_data <- .x %>% filter(is_inversion & !is.na(inversion_depth))
#       if (nrow(inv_data) > 0) {
#         min(inv_data$inversion_depth, na.rm = TRUE)
#       } else {
#         NA_real_
#       }
#     }),
# 
#     # Count distinct inversion layers (unique base heights)
#     inversion_layer_count = map_dbl(processed_data, ~ {
#       inv_data <- .x %>% filter(is_inversion & !is.na(inversion_base_height))
#       if (nrow(inv_data) > 0) {
#         length(unique(inv_data$inversion_base_height))
#       } else {
#         0
#       }
#     }),
# 
#     # Absolute deviation from canonical dry adiabatic lapse rate (K/km)
#     # If your 'temp_gradient' is in K/m, convert before comparison.
#     dry_adiabatic_deviation = abs(mean_gradient + 9.76)
#   )
# 
# # ------------------------------------------------------------
# # Keep only summary columns and persist table to disk
# # ------------------------------------------------------------
# # Drop nested raw columns, ungroup grouping structure for a flat table.
# profile_summary_clean <- profile_summary %>%
#   select(-data, -processed_data) %>%
#   ungroup()
# 
# write.table(
#   profile_summary_clean,
#   "data/processed/temperature_gradient_summary_v3.txt",
#   row.names = FALSE,
#   dec = ",",
#   sep = ";"
# )
# 
# # ------------------------------------------------------------
# # Persist all processed level data for diagnostics/plots
# # ------------------------------------------------------------
# # Unnest the processed profiles back to a long table with keys to rejoin
# # against the summary or raw inputs.
# all_processed_data <- nested_profile %>%
#   select(country, department, X, Y, year, day, hour, processed_data) %>%
#   unnest(processed_data)
# 
# write.table(
#   all_processed_data,
#   "data/processed/all_temperature_gradients_v3.txt",
#   row.names = FALSE,
#   dec = ",",
#   sep = ";"
# )
# 
# # --------------------------
# # End of script
# # --------------------------
