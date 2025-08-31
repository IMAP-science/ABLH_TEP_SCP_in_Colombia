# ======================================================================
# Temperature gradient & thermal inversion detection for one vertical profile
# ----------------------------------------------------------------------
# Purpose
#   Given a single vertical temperature profile (km above sea level + temp °C),
#   compute finite‑difference temperature gradients (°C/km) between successive
#   levels and flag inversion layers (positive dT/dz). Additionally, summarize
#   inversion intensity, start/base heights, and depth per continuous inversion.
#
# Inputs (expected columns)
#   - profile_data: data.frame / tibble with, at minimum:
#       * km_asl  : numeric, altitude in kilometers above sea level
#       * temp    : numeric, temperature in °C
#     Extra columns are allowed and preserved via dplyr pipelines.
#
# Outputs (tibble per input row minus the last, which lacks a "lead")
#   Columns added:
#     * temp_next             : lead(temp)
#     * km_asl_next           : lead(km_asl)
#     * delta_temp            : temp_next - temp (°C)
#     * delta_height          : km_asl_next - km_asl (km)
#     * temp_gradient         : delta_temp / delta_height (°C/km)
#     * is_inversion          : TRUE if temp_gradient > 0
#     * inversion_intensity   : temp_gradient when is_inversion else 0
#     * inversion_start_height: km_asl where an inversion segment starts
#     * inversion_base_height : base of each continuous inversion layer (km)
#     * inversion_depth       : layer thickness (top - base) in km
#   Notes:
#     - The final row (no lead) is dropped; NA gradients removed.
#     - If the input has < 2 rows or lacks required columns, graceful handling.
#
# Conventions & assumptions
#   - Positive dT/dz denotes inversion (warming with height).
#   - Gradients are °C/km (not K/m). Literature dry‑adiabatic reference ~‑9.76 °C/km.
#   - Input is sorted ascending by km_asl inside the function.
#   - Continuous inversion groups are detected by runs of TRUE in is_inversion.
#
# Edge cases handled
#   - Empty / 1‑row inputs -> return zero‑row tibble with full schema.
#   - Missing required columns -> stop with informative error.
#   - Profiles with no inversions -> return NA for base/depth fields.
# ======================================================================

temperature_gradient <- function(profile_data) {
  
  # ---- Input validation ------------------------------------------------
  # Empty or too short profile -> return empty tibble with full schema
  if (nrow(profile_data) < 2) {
    return(tibble(
      km_asl = numeric(0),
      temp = numeric(0),
      temp_next = numeric(0),
      km_asl_next = numeric(0),
      delta_temp = numeric(0),
      delta_height = numeric(0),
      temp_gradient = numeric(0),
      is_inversion = logical(0),
      inversion_intensity = numeric(0),
      inversion_start_height = numeric(0),
      inversion_base_height = numeric(0),
      inversion_depth = numeric(0)
    ))
  }
  
  # Required columns present?
  if (!all(c("km_asl", "temp") %in% names(profile_data))) {
    stop("Required columns 'km_asl' and 'temp' not found in profile_data")
  }
  
  # ---- Ordering & finite differences ----------------------------------
  # Sort by altitude to ensure monotonic vertical coordinate
  profile_sorted <- profile_data %>%
    arrange(km_asl)
  
  # Lead/lag deltas and gradient in °C/km; flag inversions and their starts
  result <- profile_sorted %>%
    mutate(
      temp_next   = lead(temp),
      km_asl_next = lead(km_asl),
      delta_temp  = temp_next - temp,              # °C
      delta_height= km_asl_next - km_asl,          # km
      temp_gradient = delta_temp / delta_height,   # °C/km
      is_inversion   = temp_gradient > 0,          # inversion if warming with height
      inversion_intensity = if_else(is_inversion, temp_gradient, 0),
      inversion_start_height = if_else(is_inversion, km_asl, NA_real_)
    ) %>%
    # Drop the last row where lead() produced NA and any NA gradients
    filter(!is.na(temp_gradient))
  
  # If nothing remains after filtering, return empty schema (rare but safe)
  if (nrow(result) == 0) {
    return(tibble(
      km_asl = numeric(0),
      temp = numeric(0),
      temp_next = numeric(0),
      km_asl_next = numeric(0),
      delta_temp = numeric(0),
      delta_height = numeric(0),
      temp_gradient = numeric(0),
      is_inversion = logical(0),
      inversion_intensity = numeric(0),
      inversion_start_height = numeric(0),
      inversion_base_height = numeric(0),
      inversion_depth = numeric(0)
    ))
  }
  
  # ---- Continuous inversion grouping ----------------------------------
  # Create run IDs for sequences of TRUE/FALSE in is_inversion (no wraparound)
  result <- result %>%
    mutate(
      inversion_group = {
        if (nrow(result) == 1) {
          1
        } else {
          cumsum(c(TRUE, diff(is_inversion) != 0))
        }
      }
    )
  
  # Subset to inversion rows and summarize base/top/depth per group
  inversion_data <- result %>%
    filter(is_inversion)
  
  if (nrow(inversion_data) > 0) {
    inversion_layers <- inversion_data %>%
      group_by(inversion_group) %>%
      summarise(
        inversion_base_height = min(km_asl),                     # km
        inversion_top_height  = max(km_asl_next, na.rm = TRUE),  # km
        inversion_depth       = inversion_top_height - inversion_base_height,
        .groups = "drop"
      )
    
    # Join layer summaries back; non‑inversion rows get NA for base/depth
    result <- result %>%
      left_join(inversion_layers, by = "inversion_group") %>%
      mutate(
        inversion_base_height = if_else(is_inversion, inversion_base_height, NA_real_),
        inversion_depth       = if_else(is_inversion, inversion_depth, NA_real_)
      )
  } else {
    # No inversions found -> add NA columns to keep a consistent schema
    result <- result %>%
      mutate(
        inversion_base_height = NA_real_,
        inversion_depth       = NA_real_
      )
  }
  
  # ---- Cleanup & return ------------------------------------------------
  # Drop helper columns and return tidy profile diagnostics
  result <- result %>%
    select(-inversion_group) %>%
    select(-any_of("inversion_top_height"))
  
  return(result)
}