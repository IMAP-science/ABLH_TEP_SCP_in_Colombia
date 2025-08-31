# ======================================================================
# GLM linkage of pollutants with ABLH / Ventilation Coefficient / Inversions
# ----------------------------------------------------------------------
# Purpose
#   Ingest hourly SISAIRE pollutants, ERA5 ABLH, wind‑derived VC, and
#   inversion metrics; assemble a modeling table; fit Gaussian GLMs for
#   pollutant concentration using atmospheric drivers, stratified
#   by department × pollutant × season × stability; compute predictive metrics,
#   export coefficient tables, derive scale‑aware relative effects, and map them.
#
# Key notes
#   - Seasons differ by Caribbean vs Andean departments (see logic below).
#   - Stability proxy: "unstable" for 14–21 LT, else "stable".
#   - Gradients (mean_gradient, grad_first_inversion) are converted to °C/km.
#   - GLM family is Gaussian; coefficients are later rescaled to relative effects
#     (per 10% change of covariates using their means) for interpretability.
#   - Outputs: metrics tables (txt/html), coefficients (absolute & relative),
#     and a facetted map figure summarizing relative effect magnitudes/signs.
# ======================================================================

require(pacman)
p_load(tidyverse,
       data.table,
       lubridate,
       RColorBrewer,
       gt,
       beepr,
       mgcv,
       yardstick,
       broom,
       sf,
       ggspatial,
       ggrepel,
       patchwork)

# ----------------------------------------------------------------------
# 1) Import data: ABLH (ERA5), pollutants (SISAIRE), wind (CAMS), inversions (COSMIC2)
# ----------------------------------------------------------------------
# ERA5 ABLH (hourly lumped per department)
data_abl <- fread(
  "data/processed/ABLH/ERA5/ablh_era5_lumped.txt",
  sep = ";",
  dec = ",",
  header = TRUE
)

data_abl[, date := lubridate::date(date)]

# SISAIRE hourly pollutants by station; later aggregated by dept/hour

data_sisaire <- fread(
  "data/processed/SISAIRE/data_sisaire_hourly.txt",
  sep = ";",
  dec = ",",
  header = TRUE
)

data_sisaire[, date := lubridate::date(date)]

# Hourly wind summaries (u, v, speed, dir). Will be used to compute VC.

data_winds <- fread("data/processed/WIND/wind_hourly_by_department_v1.csv")


data_winds[, date := lubridate::date(date)]

# ----------------------------------------------------------------------
# 2) Atmospheric circulation table: add Ventilation Coefficient (VC)
# ----------------------------------------------------------------------
# VC = ABLH × wind_speed.

data_atm <- data_abl %>%
  left_join(data_winds, by = c("department", "date", "hour")) %>%
  mutate(vc = wind_speed * ablh_m) %>%
  select(department, date, hour, ablh_m, vc)

# ----------------------------------------------------------------------
# 3) Inversion metrics (from XGB profiles) & modeling table assembly
# ----------------------------------------------------------------------
# Keep gradient/base/gradient‑first inversions. Convert date to Date class.

data_inversion <- read_rds("data/processed/TEMP/all_inversion_metrics_XGB.rds") %>%
  select(department,
         date,
         hour,
         mean_gradient,
         base_inversion_m,
         grad_first_inversion) %>%
  mutate(date = lubridate::date(date))

# SISAIRE -> dept/hour means per pollutant; join with ABLH/VC and inversions.
# Long pivot over pollutants to allow per‑pollutant models.

data_all_vars <- data_sisaire %>%
  select(department, date, hour, parameter, value) %>%
  arrange(department, parameter, date, hour) %>%
  group_by(department, date, hour, parameter) %>%
  reframe(mean_value = mean(value, na.rm = TRUE)) %>%
  pivot_wider(names_from = parameter, values_from = mean_value) %>%
  mutate(date = lubridate::date(date)) %>%
  full_join(data_atm) %>%
  full_join(data_inversion) %>%
  pivot_longer(
    cols = c("O3", "PM10", "PM2.5"),
    names_to = "pollutant",
    values_to = "value"
  ) %>%
  mutate(
    department = factor(department),
    pollutant = factor(pollutant),
    season = ifelse(
      department %in% c("BOLIVAR", "CESAR", "MAGDALENA"),
      ifelse(month(date) %in% c(4:5, 9:11), "wet", "dry"),
      ifelse(month(date) %in% c(3:5, 9:11), "wet", "dry")
    )
  ) %>%
  group_by(department, pollutant) %>%
  arrange(date, hour, .by_group = TRUE) %>%
  # Short ABLH memory: include lags 1,2,3,6 h
  mutate(
    ablh_lag1  = lag(ablh_m, 1),
    ablh_lag2  = lag(ablh_m, 2),
    ablh_lag3  = lag(ablh_m, 3),
    ablh_lag6 = lag(ablh_m, 6)
  ) %>%
  ungroup() %>%
  # Discard rows with any missing predictors/response
  drop_na(
    ablh_m,
    ablh_lag1,
    ablh_lag2,
    ablh_lag3,
    ablh_lag6,
    vc,
    mean_gradient,
    base_inversion_m,
    grad_first_inversion,
    value
  )

# ----------------------------------------------------------------------
# 4) Stratification: stability window and refined season definition
# ----------------------------------------------------------------------
# Stability proxy: daytime/afternoon hours as "unstable"; otherwise "stable".

data_all_vars <- data_all_vars %>%
  mutate(
    stability = ifelse(hour %in% 14:21, "unstable", "stable") %>%
      factor(levels = c("unstable", "stable")),
    season = ifelse(
      department %in% c("ANTIOQUIA", "SANTANDER", "BOGOTA_DC", "BOYACA"),
      ifelse(month(date) %in% c(3:5, 9:11), "wet", "dry"),
      ifelse(month(date) %in% c(4:5, 9:11), "wet", "dry")
    ) %>%
      factor(levels = c("wet", "dry"))
  ) %>%
  mutate(across(.cols = contains("grad"), .fns = ~ .x * 10 ^ 3))

# Persist modeling table (for reproducibility / downstream checks)
# Gradients saved in °C/km
write.table(
  data_all_vars,
  "data/processed/all_variables.txt",
  sep = ";",
  dec = ",",
  row.names = FALSE
)

# ----------------------------------------------------------------------
# 5) Modeling: Gaussian GLM per (dept × pollutant × season × stability)
# ----------------------------------------------------------------------
# Predictors: contemporaneous ABLH, short lags (1,2,3,6h), VC, mean gradient,
# inversion base height, and first inversion gradient.
# Metrics: RMSE/MAE/R² computed in‑sample on fitted values (predict on training).

models <- data_all_vars %>%
  group_by(department, pollutant, season, stability) %>%
  nest() %>%
  mutate(
    fit = map(
      data,
      ~ glm(
        value ~
          ablh_m +
          ablh_lag1 +
          ablh_lag2 +
          ablh_lag3 +
          ablh_lag6 +
          vc +
          mean_gradient +
          base_inversion_m +
          grad_first_inversion,
        data = .x,
        family = gaussian()
      )
    ),
    preds        = map2(fit, data, ~ predict(.x, .y)),
    test_metrics = map2_dfr(
      data,
      preds,
      ~ tibble(
        RMSE    = rmse_vec(.x$value, .y),
        MAE     = mae_vec(.x$value, .y),
        R2_test = rsq_vec(.x$value, .y)
      )
    ),
    # coefficients (with 95% CI)
    coefs = map(fit, tidy, conf.int = TRUE)
  )

# Save & re‑load (ensures consistent downstream reads)
write_rds(models,
          "data/processed/statistic_summaries/hourly/models.rds")

models <- read_rds("data/processed/statistic_summaries/hourly/models.rds")

# ----------------------------------------------------------------------
# 6) Performance table (RMSE/MAE/R²)
# ----------------------------------------------------------------------
perf_tbl <- models %>%
  select(department, pollutant, season, stability, test_metrics) %>%
  unnest(test_metrics) %>%
  transmute(
    department    = department,
    season        = season,
    stability     = stability,
    pollutant     = pollutant,
    RMSE          = round(RMSE, 3),
    MAE           = round(MAE, 3),
    R2_test       = round(R2_test, 3)
  ) %>%
  arrange(department, season, stability, pollutant)

# Export pretty HTML and TXT versions
perf_tbl %>%
  gt() %>%
  tab_header("GLM Performance by Department & Pollutant") %>%
  gtsave("performance_metrics.html", path = "data/processed/statistic_summaries/hourly/")

perf_tbl %>%
  write.table(
    "data/processed/statistic_summaries/hourly/performance_metrics.txt",
    sep = ";",
    dec = ",",
    row.names = FALSE
  )

# ----------------------------------------------------------------------
# 7) Coefficient table (tidy) and relative‑effect transformation
# ----------------------------------------------------------------------
coeff_tbl <- models %>%
  select(department, pollutant, coefs) %>%
  unnest(coefs) %>%
  filter(term != "(Intercept)") %>%
  transmute(
    department,
    pollutant,
    season,
    stability,
    term,
    variable = recode(
      term,
      ablh_m             = "ABLH",
      ablh_lag1          = "ABLH[1*h~lag]",
      ablh_lag2          = "ABLH[2*h~lag]",
      ablh_lag3          = "ABLH[3*h~lag]",
      ablh_lag6          = "ABLH[6*h~lag]",
      vc                 = "VC",
      mean_gradient      = "Adiabatic~Gradient",
      base_inversion_m   = "Inversion~Base~Height",
      grad_first_inversion = "Inversion~Gradient"
    ) %>%
      factor(
        levels = c(
          "ABLH",
          "ABLH[1*h~lag]",
          "ABLH[2*h~lag]",
          "ABLH[3*h~lag]",
          "ABLH[6*h~lag]",
          "VC",
          "Adiabatic~Gradient",
          "Inversion~Base~Height",
          "Inversion~Gradient"
        )
      ),
    estimate = round(estimate, 3),
    ci_lower = round(conf.low, 3),
    ci_upper = round(conf.high, 3),
    p_value  = signif(p.value, 3)
  ) %>%
  arrange(department, pollutant, variable)

# Pretty HTML + TXT
coeff_tbl %>%
  gt() %>%
  tab_header("GLM Coefficients by Department & Pollutant") %>%
  fmt_scientific(columns = p_value, scale_by = 1) %>%
  gtsave("variable_effects.html", path = "data/processed/statistic_summaries/hourly/")

coeff_tbl %>%
  write.table(
    "data/processed/statistic_summaries/hourly/variable_effects.txt",
    sep = ";",
    dec = ",",
    row.names = FALSE
  )

# ----------------------------------------------------------------------
# 8) Relative effects (scale‑aware): per 10% increase using variable means
# ----------------------------------------------------------------------
# Build means for each numeric predictor over the whole modeling table.
means <- data_all_vars %>%
  summarise(across(
    c(
      ablh_m,
      ablh_lag1,
      ablh_lag2,
      ablh_lag3,
      ablh_lag6,
      vc,
      mean_gradient,
      base_inversion_m,
      grad_first_inversion
    ),
    ~ mean(.x, na.rm = TRUE)
  )) %>%
  pivot_longer(everything(), names_to = "var_key", values_to = "mean_val")

# Map tidy terms -> means, compute 10% change effects and 95% CI in response units
coeff_rel <- coeff_tbl %>%
  mutate(
    var_key = case_when(
      term == "ablh_m"               ~ "ablh_m",
      term == "ablh_lag1"            ~ "ablh_lag1",
      term == "ablh_lag2"            ~ "ablh_lag2",
      term == "ablh_lag3"            ~ "ablh_lag3",
      term == "ablh_lag6"            ~ "ablh_lag6",
      term == "vc"                   ~ "vc",
      term == "mean_gradient"        ~ "mean_gradient",
      term == "base_inversion_m"     ~ "base_inversion_m",
      term == "grad_first_inversion" ~ "grad_first_inversion"
    )
  ) %>%
  left_join(means, by = "var_key") %>%
  mutate(
    estimate_rel = estimate * mean_val * 0.1,
    ci_lower_rel = ci_lower * mean_val * 0.1,
    ci_upper_rel = ci_upper * mean_val * 0.1
  ) %>%
  select(-c(mean_val, estimate, ci_lower, ci_upper, var_key, variable))

# Export tidy relative‑effect table (wide per pollutant × season, by stability)
coeff_rel %>%
  mutate(summary_result = paste0(
    ifelse(
      p_value < 0.05,
      paste0(formatC(
        estimate_rel, format = "f", digits = 2
      ), "*"),
      formatC(estimate_rel, format = "f", digits = 2)
    ),
    " [",
    formatC(ci_lower_rel, format = "f", digits = 2),
    " - ",
    formatC(ci_upper_rel, format = "f", digits = 2),
    "]"
  )) %>%
  select(-c(estimate_rel, ci_lower_rel, ci_upper_rel, p_value)) %>% 
  pivot_wider(
  names_from = c(pollutant, season),
  names_sep = "-",
  values_from = summary_result
) %>%
  select(department, term, stability, contains(c("PM10", "PM2.5", "O3"))) %>%
  write.table(
    x = .,
    "data/processed/statistic_summaries/hourly/variable_effects_relative.txt",
    sep = ";",
    dec = ",",
    row.names = FALSE
  )

# ----------------------------------------------------------------------
# 9) Visualization: facetted map of relative effects (sign/magnitude)
# ----------------------------------------------------------------------
# Base layers and centroids with label nudges for display
rect_study_areas <- sf::read_sf("data/boundaries/boundaries_study_area_col_v2.gpkg") |>
  dplyr::select(-dplyr::any_of("grupo"))

rect_colombia_lv0 <- sf::read_sf("../E-SIG/SHP-Limites administrativos COL - GADM/gadm41_COL_0.shp") |>
  sf::st_transform(crs = sf::st_crs(rect_study_areas))

rect_urban_centers <- sf::read_sf("../E-SIG/MGN2024_URB_ZONA_URBANA/MGN_URB_ZONA_URBANA.shp") |>
  sf::st_transform(crs = sf::st_crs(rect_study_areas))

# Centroids (nudge labels vertically per department to avoid overlaps)
rect_study_areas_centroids <- rect_study_areas %>%
  sf::st_centroid(of_largest_polygon = TRUE) %>%
  cbind(sf::st_coordinates(.)) %>%
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
  ) %>% 
  dplyr::mutate(
    department_fixed = dplyr::recode(department,
                                                 ANTIOQUIA = "Antioquia",
                                                 BOGOTA_DC = "Bogotá",
                                                 BOLIVAR   = "Bolívar",
                                                 BOYACA    = "Boyacá",
                                                 CALDAS    = "Caldas",
                                                 CESAR     = "Cesar",
                                                 MAGDALENA = "Magdalena",
                                                 NORTE_DE_SANTANDER = "Norte de Santander",
                                                 SANTANDER = "Santander"
                                                 
  ))

# Using relative effects for mapping: sign via fill, magnitude via size.
coeff_map <- coeff_rel %>%
  dplyr::filter(term %in% c("ablh_m","vc","mean_gradient","base_inversion_m","grad_first_inversion")) %>%
  dplyr::transmute(
    department, pollutant, season, stability, term,
    rel_effect = estimate_rel,
    p_value    = p_value
  ) %>%
  dplyr::mutate(
    variable = dplyr::recode(term,
                             ablh_m               = "ABLH",
                             vc                   = "VC",
                             mean_gradient        = "Adiabatic Gradient",
                             base_inversion_m     = "Inversion Base Height",
                             grad_first_inversion = "Inversion Gradient"
    ),
    variable  = factor(variable,
                       levels = c("ABLH","VC","Adiabatic Gradient",
                                  "Inversion Base Height","Inversion Gradient")),
    season    = factor(season, levels = c("wet","dry")),
    stability = factor(stability, levels = c("unstable","stable")),
    col_key   = paste0(pollutant, " – ", season)
  ) %>%
  dplyr::left_join(
    rect_study_areas_centroids %>% sf::st_drop_geometry(),
    by = "department"
  ) %>%
  dplyr::select(- term)

# Map: each facet shows variable (rows) × pollutant–season (cols)
# Point size encodes |relative effect| (µg/m³ per 10% increase in driver).
# Fill color encodes sign (blue positive, red negative). Only p<0.05 shown.

graph_model_coefficients <- ggplot() +
  geom_sf(data = rect_colombia_lv0, linewidth = 1) +
  
  geom_text(
    data = rect_study_areas_centroids %>% 
      st_drop_geometry(),
    aes(x = X, y = Y_corrected, label = department_fixed),
    family = "serif",
    fontface = "bold.italic"
  ) +
  
  geom_point(data = coeff_map %>%
               mutate(slope_sign = if_else(rel_effect >= 0, "positive", "negative")) %>%
               filter(p_value < .05),
             aes(x = X, y = Y, 
                 shape  = stability,
                 fill   = slope_sign,
                 size   = abs(rel_effect)),
             position   = position_dodge(width = .7),
             alpha = .6) +
  
  scale_shape_manual(name = "",
                     breaks = c("stable", "unstable"),
                     values = c(21, 23),
                     labels = c("Stable", "Unstable")) +
  
  scale_size_area(
    max_size = 16,
    limits = c(0, NA),
    name = expression("|Rel. effect| ("* µ *"g/m"^3*" per 10% increase)")
  ) +
  
  scale_fill_manual(
    values = c(positive = "blue", negative = "red"),
    labels = c(positive = "Positive", negative = "Negative"),
    name = ""
  ) +
  
  scale_y_continuous(breaks = seq(5, 11, 2)) +
  scale_x_continuous(breaks = seq(-76, -72, 2)) +
  
  coord_sf(xlim = c(-76.5, -71.5), ylim = c(4.35 , 11.4)) +
  
  facet_grid(rows = vars(variable),
             cols = vars(col_key),
             scales = "fixed") +
  
  theme_bw() +
  theme(
    text = element_text(family = "serif"),
    axis.text = element_text(size = 15),
    axis.title = element_blank(),
    panel.background = element_rect(fill = "#74beea"),
    panel.grid = element_line(color = "#e0e0e0"),
    panel.spacing.x = unit(2, "pt"),
    panel.spacing.y = unit(2, "pt"),
    legend.position = "bottom",
    legend.title = element_text(size = 15, hjust = .5),
    legend.title.position = "top",
    legend.text = element_text(hjust = .5, size = 15),
    legend.key.width = unit(1, "cm"),
    legend.key = element_rect(fill = "transparent", color = NA),
    plot.title = element_text(hjust = .5),
    plot.margin = margin(2, 5, 0, 0),
    panel.spacing = unit(0, "pt"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15,
                              color = "black")
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 21,
                                                 size = 8,
                                                 colour = "black")),
         shape = guide_legend(override.aes = list(size = 8,
                                                  stroke = 2,
                                                  colour = "black")),
         size = guide_legend(nrow = 1))

# Save multi‑format map exports
ggsave(
  "graphs/hourly/potential/models_map.jpg",
  graph_model_coefficients,
  width = 15,
  height = 15,
  dpi = 500
)
ggsave(
  "graphs/hourly/potential/models_map.pdf",
  graph_model_coefficients,
  width = 15,
  height = 15,
  dpi = 500
)
ggsave(
  "graphs/hourly/potential/models_map.tiff",
  graph_model_coefficients,
  width = 15,
  height = 15,
  dpi = 500
)