# ======================================================================
# COSMIC‑2 -> XGBoost vertical profiles + inversion metrics (R)
# ----------------------------------------------------------------------
# Purpose
#   1) Load filtered COSMIC‑2 temperature profiles over Colombia.
#   2) Engineer features and train an XGBoost regression model for temperature
#      as a function of (X, Y, altitude, day‑of‑year, hour‑of‑day).
#   3) Generate hourly vertical temperature profiles on department centroids.
#   4) Persist per‑department profiles and a combined RDS.
#   5) Compute vertical gradients and simple inversion diagnostics per hour.
# ======================================================================

# ---- Libraries ----
# Pacman simplifies conditional loads; packages cover spatial (sf), time (lubridate),
# fast tables (data.table), modeling (xgboost) and tidy ops (tidyverse).
# No side effects beyond loading.
require(pacman)
p_load(tidyverse, data.table, sf, lubridate, xgboost)

# ----------------------------------------------------------------------
# Discover COSMIC‑2 files and filter to Colombia subset
# ----------------------------------------------------------------------
# 'list.files(..., full.names = TRUE)' keeps absolute/relative paths.
# Filtering by pattern assumes filenames contain the country token "Colombia".
cosmic2_files <- list.files("data/clean_data/cosmic2_filtered_per_country/",
                            full.names = TRUE)

# Keep only files corresponding to Colombia (name contains "Colombia")
cosmic2_files <- cosmic2_files[grep("Colombia", cosmic2_files)]

# ----------------------------------------------------------------------
# Read all files as sf, enforce WGS84, extract coordinates, drop geometry
# ----------------------------------------------------------------------
# read_sf -> sf object with geometry; st_as_sf(crs=4326) ensures CRS defined.
# st_coordinates(.) returns matrix [X,Y]; bind_cols appends numeric lon/lat
# for modeling; st_drop_geometry() flattens to a rectangular table.
cosmic2_data <- map(.x = cosmic2_files, .f = ~ read_sf(.x)) %>%
  bind_rows() %>%
  st_as_sf(crs = 4326) %>%
  bind_cols(., st_coordinates(.)) %>% 
  st_drop_geometry()

# Convert to data.table for memory/speed efficiency (by‑reference ops)
cosmic2_data_v1 <- data.table(cosmic2_data) 

# Keep only needed fields; convert km_asl -> meters (altitude_m)
# Using scientific notation 10^3 for clarity; remove original km column.
cosmic2_data_v1 <- 
  cosmic2_data_v1[,.(year, day, hour, km_asl, temp, X, Y)][, altitude_m := km_asl * 10^3][, - "km_asl"]

# Ensure data.table by‑reference behavior
setDT(cosmic2_data_v1)

# ----------------------------------------------------------------------
# Build date and datetime (UTC) from year + day‑of‑year + hour
# ----------------------------------------------------------------------
# Day‑of‑year to Date: ymd(Y-01-01) + (doy-1); then compose POSIXct hour stamps.
cosmic2_data_v1[, date := ymd(paste0(year, "-01-01")) + days(day - 1)]
cosmic2_data_v1[
  , datetime := as.POSIXct(
    paste0(date, " ", sprintf("%02d:00:00", hour)),
    tz = "UTC"
  )
]

# Feature engineering for modeling: hour‑of‑day and day‑of‑year (cyclic proxies)
# Note: for strict cyclicity one could add sin/cos transforms; here raw indices suffice.
cosmic2_data_v1[
  , ':='(
    hour_of_day = hour(datetime),
    day_of_year = yday(datetime)
  )
]

# ----------------------------------------------------------------------
# Build department centroids as prediction grid (X/Y only)
# ----------------------------------------------------------------------
# The centroid per department acts as a representative location for vertical
# profile generation. Coordinates are extracted into plain numeric columns.
study_sf <- read_sf("data/boundaries/boundaries_study_area_col_v2.gpkg")
centroids <- study_sf %>%
  st_centroid() %>%
  transmute(
    department,
    X = st_coordinates(.)[,1],
    Y = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry() %>%
  as.data.table()

# ----------------------------------------------------------------------
# Assemble training matrix for XGBoost
# ----------------------------------------------------------------------
# Feature set: longitude (X), latitude (Y), altitude in meters, and temporal
# indices (day_of_year, hour_of_day).
features <- c("X","Y","altitude_m","day_of_year","hour_of_day")

# xgb.DMatrix requires numeric matrix (features) and numeric label (temp)
dtrain <- xgb.DMatrix(
  data  = as.matrix(cosmic2_data_v1[, ..features]),
  label = cosmic2_data_v1$temp
)

# XGBoost hyperparameters:
#  - tree_method='hist' for large datasets
#  - eta=0.05 (learning rate) with nrounds=500 (moderate depth=10)
#  - subsample/colsample control variance & speed; threads = cores-4 for margin.
params <- list(
  objective        = "reg:squarederror",
  tree_method      = "hist",
  max_depth        = 10,
  eta              = 0.05,
  subsample        = 0.8,
  colsample_bynode = 0.8,
  nthread          = parallel::detectCores() - 4
)

# Train model with a fixed seed for reproducibility
set.seed(2025)
xgb_mod <- xgb.train(
  params   = params,
  data     = dtrain,
  nrounds  = 500,
  verbose  = 1
)

# Free training matrix to reduce RAM usage (model is kept)
rm(dtrain); gc()

# ----------------------------------------------------------------------
# Define prediction grids: altitude sequence and hourly time sequence
# ----------------------------------------------------------------------
# Vertical resolution: 50 m steps across observed altitude range.
# Temporal resolution: every hour spanning min->max observation timestamps.
alt_rng <- range(cosmic2_data_v1$altitude_m, na.rm = TRUE)
z_seq   <- seq(floor(alt_rng[1]/50)*50,
               ceiling(alt_rng[2]/50)*50,
               by = 50)

t_seq <- seq(
  from = floor_date(min(cosmic2_data_v1$datetime),  "hour"),
  to   = ceiling_date(max(cosmic2_data_v1$datetime), "hour"),
  by   = "1 hour"
)

# Output directory for predicted profiles; clean previous outputs for a fresh run
# (removes old per‑department RDS files to avoid mixing runs)
out_dir <- "data/processed/TEMP/xgb_profiles"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

old_rds <- list.files(out_dir,
                      pattern    = "^vertical_profile_XGB_.*\\.rds$",
                      full.names = TRUE)
if (length(old_rds)) unlink(old_rds)

# ----------------------------------------------------------------------
# Generate hourly vertical profiles per department centroid
# ----------------------------------------------------------------------
# For each centroid: build Cartesian grid (CJ) of datetime × altitude,
# derive temporal features, predict temperature, and persist compressed RDS.
for (i in seq_len(nrow(centroids))) {
  dept <- centroids$department[i]
  X0   <- centroids$X[i]
  Y0   <- centroids$Y[i]
  
  # Build prediction grid for one department across altitude and time
  newdata <- CJ(
    datetime   = t_seq,
    altitude_m = z_seq
  )[, ':='(
    X           = X0,
    Y           = Y0,
    day_of_year = yday(datetime),
    hour_of_day = hour(datetime)
  )]
  
  # Predict temperature on the grid
  preds <- predict(
    xgb_mod,
    as.matrix(newdata[, ..features])
  )
  newdata[, temp_pred := preds]
  newdata[, department := dept]
  
  # Arrange columns for clarity (dept, time, height, predicted temp)
  setcolorder(
    newdata,
    c("department", "datetime", "altitude_m", "temp_pred")
  )
  
  # Save department level profile (compressed RDS)
  saveRDS(
    newdata,
    file     = file.path(out_dir,
                         sprintf("vertical_profile_XGB_%s.rds", dept)),
    compress = "xz"
  )
  
  message("Perfil XGB guardado para ", dept)
  rm(newdata, preds); gc()
}

# ----------------------------------------------------------------------
# Combine all department profiles into a single object for downstream use
# ----------------------------------------------------------------------
# rbindlist reads each per‑dept RDS and stacks into one long table.
rds_files <- list.files(out_dir,
                        pattern    = "^vertical_profile_XGB_.*\\.rds$",
                        full.names = TRUE)
all_profiles <- rbindlist(
  lapply(rds_files, readRDS),
  use.names = TRUE
)

saveRDS(
  all_profiles,
  file     = file.path(out_dir, "all_vertical_profiles_XGB.rds"),
  compress = "xz"
)
message("Todos los perfiles combinados en all_vertical_profiles_XGB.rds")

# ----------------------------------------------------------------------
# Compute vertical gradient and basic inversion metrics (per dept × date × hour)
# ----------------------------------------------------------------------
# Order by altitude (ascending) before finite differences; compute dT/dz.
# Inversion criterion: positive gradient (>0) between successive levels.
# Basic diagnostics include: mean gradient across column, count of inversions,
# first inversion base/top/depth, and gradient at first inversion.
all_profiles[
  , ':='(
    date = as.Date(datetime),
    hour = hour(datetime)
  )
]

inversion_metrics <- all_profiles[
  # Ensure records are ordered by altitude before computing finite differences
  order(altitude_m),
  {
    alt   <- altitude_m
    tmp   <- temp_pred
    alt_n <- shift(alt, type = "lead")
    tmp_n <- shift(tmp, type = "lead")
    grad  <- (tmp_n - tmp) / (alt_n - alt)
    
    # Identify layers with positive gradient as inversions
    inv_idx      <- which(grad > 0)
    n_inversions <- length(inv_idx)
    has_inv      <- n_inversions > 0
    
    # First inversion layer diagnostics
    base1  <- if (has_inv) alt[inv_idx[1]]      else NA_real_
    top1   <- if (has_inv) alt_n[inv_idx[1]]    else NA_real_
    depth1 <- if (has_inv) top1 - base1         else NA_real_
    grad1  <- if (has_inv) grad[inv_idx[1]]     else NA_real_
    mean_g <- mean(grad, na.rm = TRUE)
    
    .(
      mean_gradient           = mean_g,
      n_inversions            = n_inversions,
      has_inversion           = has_inv,
      base_inversion_m         = base1,
      top_inversion_m          = top1,
      depth_inversion_m        = depth1,
      grad_first_inversion     = grad1
    )
  },
  by = .(department, date, hour)
]

# Persist inversion diagnostics (compressed RDS)
saveRDS(
  inversion_metrics,
  file     = file.path(out_dir, "all_inversion_metrics_XGB.rds"),
  compress = "xz"
)
message("Métricas de inversión guardadas en all_inversion_metrics_XGB.rds")
