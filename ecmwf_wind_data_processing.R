# =============================================================
# ERA5 wind fields -> hourly departmental summaries (R)
# -------------------------------------------------------------
# Purpose
#   Process ERA5 U/V wind components from yearly GRIBs and produce
#   per‑department hourly means of u, v, wind speed and direction.
#   Uses a regular interpolation grid + gstat predictions, assigns
#   grid cells to departments via bounding‑box fuzzy join, and writes
#   block CSVs (10‑day chunks) plus a consolidated output.
#
# Key
#   - Input GRIBs contain alternating U/V layers at hourly cadence.
#   - Layers are ordered [U1, V1, U2, V2, ...]; indices split accordingly.
#   - Output CRS is WGS84 (EPSG:4326). Study polygons are transformed to it.
#   - Interpolation: gstat with constant mean surfaces (no variogram model
#     specified). This approximates a smooth surface but is not full kriging.
#   - Department assignment uses bounding boxes for speed (approximation!).
#     For exact polygon containment, consider st_join with st_within.
#   - Memory: processed by 10‑day blocks to keep RAM bounded.
#   - Parallel: furrr::future_map via multisession; progressr enabled.
#   - File I/O: Intermediate block CSVs -> final consolidated CSV.
# =============================================================

# Process ERA5 wind fields into hourly departmental summaries
# Inputs
#   years: vector of years to process
#   extent_sf: sf polygons with a column 'department' defining areas of aggregation
#   input_dir: directory containing yearly GRIB files named 'wind_data_<year>.grib'
#   output_dir: directory to write intermediate block csv files and the final aggregation
#   x_res, y_res: target grid resolution in degrees for interpolation grid
#   n_cores: number of parallel workers used by furrr
process_daily_wind_data <- function(years = 2020:2024,
                                    extent_sf,
                                    input_dir = "data/clean_data/era5/",
                                    output_dir = "data/processed/wind/",
                                    x_res = 0.15,
                                    y_res = 0.15,
                                    n_cores = 4) {
  require(pacman)
  p_load(tidyverse,
         sf,
         terra,
         fuzzyjoin,
         lubridate,
         gstat,
         furrr,
         progressr)
  
  # Parallel and progress configuration
  # - multisession: spawns background R sessions (Windows‑friendly)
  # - handlers: enable progressr notifications (console by default)
  plan(multisession, workers = n_cores)
  handlers(global = TRUE)
  
  # Ensure output directory exists (create recursively if needed)
  if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE)
  
  # Normalize study extent to WGS84 and extract per‑polygon bounding boxes
  # Resulting extent_df is a light table for fast bbox joins
  extent_sf <- st_transform(extent_sf, crs = 4326)
  extent_df <- extent_sf %>%
    mutate(bounding_box = map(geom, st_bbox)) %>%
    mutate(
      xmin = map_dbl(bounding_box, ~ .x[["xmin"]]),
      xmax = map_dbl(bounding_box, ~ .x[["xmax"]]),
      ymin = map_dbl(bounding_box, ~ .x[["ymin"]]),
      ymax = map_dbl(bounding_box, ~ .x[["ymax"]])
    ) %>%
    select(department, xmin, xmax, ymin, ymax) %>%
    st_drop_geometry()
  
  # Iterate per year to keep resource usage predictable
  walk(years, function(year) {
    file_path <- file.path(input_dir, paste0("wind_data_", year, ".grib"))
    if (!file.exists(file_path)) {
      cat("File not found:", file_path, "\n")
      return(NULL)
    }
    
    # Load GRIB as SpatRaster; align CRS to WGS84 if needed
    r <- terra::rast(file_path)
    if (terra::crs(r) != terra::crs("EPSG:4326"))
      r <- terra::project(r, "EPSG:4326")
    
    # Separate U and V components by layer index and get timestamps
    # Assumes interleaved ordering U,V,U,V,... across time
    times <- terra::time(r)
    idx_u <- seq(1, nlyr(r), by = 2)
    idx_v <- seq(2, nlyr(r), by = 2)
    r_u <- r[[idx_u]]
    r_v <- r[[idx_v]]
    
    # Group unique days into blocks of 10 days (chunking strategy)
    days <- as.Date(times[idx_u])
    day_blocks <- split(unique(days), ceiling(as.numeric(as.factor(unique(
      days
    ))) / 10))
    
    cat("Processing year:", year, "\n")
    
    # Process each block independently
    walk(day_blocks, function(day_block) {
      cat("  Processing block:",
          min(day_block),
          "to",
          max(day_block),
          "\n")
      idx_block <- which(as.Date(times[idx_u]) %in% day_block)
      if (length(idx_block) == 0)
        return(NULL)
      
      # Slice U and V rasters for the block; collect block datetimes
      r_u_block <- r_u[[idx_block]]
      r_v_block <- r_v[[idx_block]]
      datetimes <- times[idx_u][idx_block]
      
      # Build regular interpolation grid over full raster extent
      # Note: grid resolution controlled by x_res/y_res (degrees)
      grid_df <- expand.grid(x = seq(ext(r)[1], ext(r)[2], by = x_res),
                             y = seq(ext(r)[3], ext(r)[4], by = y_res))
      
      # For each timestamp in the block:
      #  1) Extract raster cell centers + U,V values
      #  2) Fit gstat constant‑mean models for u and v (no variogram)
      #  3) Predict onto the regular grid
      wind_data <- map_dfr(seq_along(datetimes), function(i) {
        u_vals <- terra::values(r_u_block[[i]])
        v_vals <- terra::values(r_v_block[[i]])
        xy <- terra::xyFromCell(r_u_block[[i]], 1:ncell(r_u_block[[i]]))
        df <- tibble(
          x = xy[, 1],
          y = xy[, 2],
          u = u_vals,
          v = v_vals,
          datetime = datetimes[i]
        )
        df <- df %>% filter(complete.cases(.))
        if (nrow(df) < 10)
          return(NULL)
        
        # gstat setup: two formulas (u and v) on same locations; simple surface
        gdf <- gstat::gstat(formula = u ~ 1,
                            locations = ~ x + y,
                            data = df)
        gdf <- gstat::gstat(
          gdf,
          formula = v ~ 1,
          locations = ~ x + y,
          data = df
        )
        pred <- predict(gdf, newdata = grid_df)
        
        tibble(
          datetime = datetimes[i],
          x = pred$x,
          y = pred$y,
          u = pred$var1.pred,
          v = pred$var2.pred
        )
      })
      
      if (nrow(wind_data) == 0)
        return(NULL)
      
      # Assign grid points to departments via bbox fuzzy join (fast, approximate)
      # If edge precision matters, replace with st_join using polygons.
      assigned <- fuzzy_left_join(
        wind_data,
        extent_df,
        by = c(
          "x" = "xmin",
          "x" = "xmax",
          "y" = "ymin",
          "y" = "ymax"
        ),
        match_fun = list(`>=`, `<=`, `>=`, `<=`)
      ) %>% drop_na()
      
      if (!"department" %in% names(assigned))
        return(NULL)
      
      # Aggregate to department × hour; compute u,v means, speed & direction
      summarized <- assigned %>%
        group_by(department, datetime) %>%
        summarise(
          u = mean(u, na.rm = TRUE),
          v = mean(v, na.rm = TRUE),
          wind_speed = sqrt(u ^ 2 + v ^ 2),
          wind_dir = (atan2(v, u) * 180 / pi + 360) %% 360,
          .groups = "drop"
        )
      
      # Persist block CSV for incremental processing/recovery
      file_suffix <- paste0(year, "_", format(min(day_block), "%j"))
      write_csv(summarized, file.path(output_dir, paste0("wind_", file_suffix, ".csv")))
      
      # Clean up block objects to free memory
      rm(r_u_block, r_v_block, wind_data, assigned, summarized)
      gc()
    })
    
    # Year‑level cleanup
    rm(r, r_u, r_v, times, idx_u, idx_v)
    gc()
  })
  
  # Consolidate all block outputs to a single CSV (department × hour)
  block_files <- list.files(output_dir, pattern = "^wind_.*\\.csv$", full.names = TRUE)
  final_result <- map_dfr(block_files, read_csv, show_col_types = FALSE)
  write_csv(final_result,
            file.path(output_dir, "wind_hourly_by_department.csv"))
  
  cat("\nProcess completed. Results saved to:", output_dir, "\n")
  return(final_result)
}


# -------------------------------------------------------------
# EXECUTION BLOCK
# -------------------------------------------------------------
# Load dependencies
require(pacman)
p_load(tidyverse,
       sf,
       terra,
       fuzzyjoin,
       lubridate,
       gstat,
       furrr,
       progressr)

# Read study areas and transform to WGS84 (EPSG:4326)
rect_study_areas <- read_sf("data/boundaries/boundaries_study_area_col_v2.gpkg") %>%
  st_transform(crs = 4326)

# Run the processing pipeline for default years and directories
process_daily_wind_data(extent_sf = rect_study_areas)

# Post processing and reshaping of the consolidated output
w <- read_csv("data/processed/wind/wind_hourly_by_department.csv")

names(w)

# Add date & hour, pivot longer -> average duplicates (safeguard) -> wider by variable
w_v1 <- w %>% mutate(date = lubridate::date(datetime),
                     hour = lubridate::hour(datetime)) %>% select(-datetime) %>%
  pivot_longer(
    cols = -c(department, date, hour),
    names_to = "variable",
    values_to = "value"
  ) %>%
  group_by(department, date, hour, variable) %>%
  reframe(value_mean = mean(value, na.rm = TRUE)) %>%
  pivot_wider(names_from = variable, values_from = value_mean)

# Simple completeness check by department and year (row counts per year)
w_v1 %>% group_by(department, year(date)) %>% reframe(n = n()) %>%
  distinct(n)

# Persist the reshaped summary table
write.csv(w_v1, "data/processed/wind/wind_hourly_by_department_v1.csv",
          row.names = FALSE)