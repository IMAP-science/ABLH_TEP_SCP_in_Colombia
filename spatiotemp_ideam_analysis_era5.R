# ======================================================================
# IDEAM / SISAIRE: Hourly Pollutant Analysis & Atmosphere Coupling (R)
# ----------------------------------------------------------------------
# Purpose
#   Analyze hourly SISAIRE pollutants reported by IDEAM across departments,
#   compute robust statistics, exceedance rates, normality/stationarity tests,
#   and build a consolidated summary together with atmospheric variables
#   (ABLH, Ventilation Coefficient VC, inversion metrics).
#
# What this script does (overview)
#   1) Loads hourly SISAIRE, ABLH (ERA5), wind summaries, and inversion metrics.
#   2) Computes per‑department/parameter/year descriptive statistics (mean, SD,
#      min/Q1/median/Q3/max, kurtosis, skewness, CV) and persists results.
#   3) Calculates daily exceedance rates against simple thresholds and summarizes
#      by station and department.
#   4) Tests normality (Lilliefors/KS) and stationarity (ADF) per station×parameter.
#   5) Builds an atmospheric dataset (ABLH, VC, gradients) and harmonizes units.
#   6) Produces consolidated daily summaries (means/min/max and 98th percentile)
#      for pollutants + atmospheric variables and writes them to disk.
#   7) Visualizes monthly distributions of extreme observations (> 98th pct).
# ======================================================================

# THIS SCRIPT HAS THE PURPOSE OF ANALYZING THE DATA FOR ALL POLLUTANTS REPORTED BY IDEAM

require(pacman)
p_load(
  tidyverse,
  data.table,
  readxl,
  lubridate,
  broom,
  moments,
  tseries,
  boot,
  nortest,
  stringi,
  ggridges,
  gt,
  gtsummary,
  RColorBrewer
)


# ----------------------------------------------------------------------
# Load inputs
# ----------------------------------------------------------------------
# SISAIRE hourly pollutants (wide coverage across stations/departments)
data_sisaire <- fread(
  "data/processed/SISAIRE/data_sisaire_hourly.txt",
  sep = ";",
  dec = ",",
  header = TRUE
) %>% mutate(date = ymd(date))

# ERA5 ABLH (hourly, lumped/aggregated)
data_abl <- fread(
  "data/processed/ABLH/ERA5/ablh_era5_lumped.txt",
  sep = ";",
  dec = ",",
  header = TRUE
)

# Hourly wind summaries by department (u, v, speed, direction)
data_wind <- fread("data/processed/WIND/wind_hourly_by_department_v1.csv")

# Inversion metrics from XGB vertical profiles (precomputed)
# Keep relevant variables for coupling with SISAIRE/ABLH

data_inversion <- read_rds("data/processed/TEMP/all_inversion_metrics_XGB.rds")
data_inversion <- data_inversion %>%
  select(department,
         date,
         hour,
         mean_gradient,
         base_inversion_m,
         grad_first_inversion)

# ----------------------------------------------------------------------
# Descriptive statistics per department × parameter × year (robust summaries)
# ----------------------------------------------------------------------
# Metrics include: mean, sd, min, quartiles, median, max, kurtosis, skewness,
# and coefficient of variation (cv = sd/mean). Persist to disk.
{
  data_sisaire %>%
    mutate(month = month(date)) %>%
    group_by(department, parameter, year) %>%
    reframe(
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      min = min(value, na.rm = TRUE),
      q1 = quantile(value, .25, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      q3 = quantile(value, .75, na.rm = TRUE),
      max = max(value, na.rm = TRUE),
      kurtosis = kurtosis(value, na.rm = TRUE),
      skewness = skewness(value, na.rm = TRUE)
    ) %>%
    mutate(cv = sd / mean) %>%
    write.table(
      .,
      "data/processed/statistic_summaries/hourly/sisaire_metrics_dep_par_year.txt",
      sep = ";",
      dec = ",",
      row.names = FALSE
    )
  }


# ----------------------------------------------------------------------
# Exceedances: daily means vs. simple per‑parameter thresholds
# ----------------------------------------------------------------------
# 1) Compute daily mean by station.
# 2) Mark daily exceedance (value > threshold) using per‑pollutant cutoffs.
# 3) Summarize exceedance rates per station, then aggregate by department.

data_excedences <- data_sisaire %>%
  group_by(department, parameter, station, year, date) %>%
  reframe(value = mean(value, na.rm = TRUE)) %>%
  mutate(
    threshold = case_when(parameter == "PM10" ~ 75L, parameter == "PM2.5" ~ 37L, TRUE ~ 100L),
    exceeds = ifelse(value > threshold, TRUE, FALSE)
  ) %>%
  group_by(department, parameter, station, year) %>%
  reframe(rate = sum(exceeds) / n()) %>%
  group_by(department, parameter) %>%
  reframe(
    median_rate = median(rate),
    q1 = quantile(rate, .25),
    q3 = quantile(rate, .75),
    p98 = quantile(rate, .98),
    mean_rate = mean(rate),
    sd_rate = sd(rate),
    n_stations = n_distinct(station)
  )

# Persist exceedance summaries
data_excedences %>%
  write.table(
    "data/processed/statistic_summaries/hourly/sisaire_stations_and_excedences.txt",
    sep = ";",
    dec = ",",
    row.names = FALSE
  )
  

# ----------------------------------------------------------------------
# Normality test (Lilliefors / Kolmogorov–Smirnov)
# ----------------------------------------------------------------------
# Null hypothesis: data follows a normal distribution with unknown mean/SD.
# We list any station×parameter combos with p>=0.05 (rare per prior findings).

data_sisaire %>%
  group_by(department, station, parameter) %>%
  summarise(ks_p_value = lillie.test(value)$p.value, .groups = "drop") %>%
  filter(ks_p_value >= 0.05)

# There is not a single law-comparable station-parameter combination with normality
# on its data.

# ----------------------------------------------------------------------
# Stationarity test (ADF with k=0)
# ----------------------------------------------------------------------
# Null hypothesis: unit root (non‑stationary). We list groups with p>=0.05.

data_sisaire %>%
  nest(data = -c(department, station, parameter)) %>%
  mutate(df_test = map(.x = data, .f = ~ adf.test(.x$value, k = 0)),
         df_p_value = map_dbl(df_test, ~ .x[["p.value"]])) %>%
  filter(df_p_value >= 0.05)

# Same with this
# SUMMARY TABLE OF THE STATIONARY BEHAVIOR OF THE DATA

# ----------------------------------------------------------------------
# Couple atmospheric variables and build consolidated daily summaries
# ----------------------------------------------------------------------
# 1) Combine ABLH + wind to derive VC (Ventilation Coefficient = ABLH × wind_speed).
# 2) Join inversion metrics (convert gradients to °C/km).
# 3) Summarize daily min/max/mean and then aggregate year‑round statistics.

data_atmosphere <- data_abl %>%
  full_join(data_wind) %>%
  mutate(vc = ablh_m * wind_speed) %>%
  select(-c(u, v, wind_dir, wind_speed)) %>%
  mutate(date = lubridate::date(date)) %>%
  full_join(data_inversion %>% mutate(date = lubridate::date(date)) %>% 
              mutate(across(.cols = contains("grad"),
                            .fns = ~ .x * 10^3))) %>%
  pivot_longer(
    cols = -c(department, date, hour),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  group_by(department, parameter, date) %>%
  reframe(
    daily_min = min(value),
    daily_max = max(value),
    daily_mean = mean(value)
  ) %>%
  group_by(department, parameter) %>%
  reframe(
    value_mean = mean(daily_mean, na.rm = TRUE),
    value_sd = sd(daily_mean, na.rm = TRUE),
    max_mean = mean(daily_max, na.rm = TRUE),
    max_sd = sd(daily_max, na.rm = TRUE),
    min_mean = mean(daily_min, na.rm = TRUE),
    min_sd = sd(daily_min, na.rm = TRUE),
    p98 = quantile(daily_mean, .98 , na.rm = TRUE)
  )

# Pollutant daily summaries aggregated similarly

data_summary <- data_sisaire %>%
  select(department, station, parameter, date, hour, value) %>%
  group_by(department, station, parameter, date) %>%
  reframe(
    daily_min = min(value),
    daily_max = max(value),
    daily_mean = mean(value)
  ) %>%
  group_by(department, parameter) %>%
  reframe(
    value_mean = mean(daily_mean, na.rm = TRUE),
    value_sd = sd(daily_mean, na.rm = TRUE),
    max_mean = mean(daily_max, na.rm = TRUE),
    max_sd = sd(daily_max, na.rm = TRUE),
    min_mean = mean(daily_min, na.rm = TRUE),
    min_sd = sd(daily_min, na.rm = TRUE),
    p98 = quantile(daily_mean, .98 , na.rm = TRUE)
  ) %>%
  bind_rows(data_atmosphere) %>%
  arrange(department, parameter)

# Persist consolidated summary
write.table(
  data_summary,
  "data/processed/statistic_summaries/hourly/data_summary.txt",
  sep = ";",
  dec = ",",
  row.names = FALSE
)

# ----------------------------------------------------------------------
# Tidy the summary for presentation: combine mean±sd strings & reshape
# ----------------------------------------------------------------------
# Creates columns "mean_and_sd", "mean_max_and_sd", "mean_min_and_sd" and keeps p98.
# Then pivots to a wide table keyed by parameter for easier inspection.

data_summary %>% 
  mutate(mean_and_sd = case_when(parameter %in% c("grad_first_inversion", "mean_gradient") ~ 
                                   paste0(round(value_mean, 2), " - ", round(value_sd, 2)),
                                 TRUE ~ paste0(round(value_mean, 2), " - ", round(value_sd, 2))),
         mean_max_and_sd = case_when(parameter %in% c("grad_first_inversion", "mean_gradient") ~ 
                                       paste0(round(max_mean, 2), " - ", round(max_sd, 2)),
                                     TRUE ~ paste0(round(max_mean, 2), " - ", round(max_sd, 2))),
         mean_min_and_sd = case_when(parameter %in% c("grad_first_inversion", "mean_gradient") ~ 
                                       paste0(round(min_mean, 2), " - ", round(min_sd, 2)),
                                     TRUE ~ paste0(round(min_mean, 2), " - ", round(min_sd, 2))),
         p98 = as.character(round(p98,2))) %>% 
  select(-c(value_mean, value_sd,
            max_mean, max_sd,
            min_mean, min_sd)) %>% 
  pivot_longer(cols = - c(department, parameter),
               names_to = "names",
               values_to = "values") %>% 
  pivot_wider(names_from = parameter,
              values_from = values) %>% 
  write.table(
    .,
    "data/processed/statistic_summaries/hourly/data_summary_v1.txt",
    sep = ";",
    dec = ",",
    row.names = FALSE
  )


# ----------------------------------------------------------------------
# Extreme values behavior (monthly): > 98th percentile counts per month
# ----------------------------------------------------------------------
# For each department × parameter × year: compute the 98th percentile threshold,
# flag hourly extremes, then count extremes by month and visualize distributions
# with ridgeline densities.

department_labels <- c(
  "ANTIOQUIA" = "Antioquia",
  "BOGOTA_DC" = "Bogotá",
  "BOLIVAR" = "Bolívar",
  "BOYACA" = "Boyacá",
  "CALDAS" = "Caldas",
  "CESAR" = "Cesar",
  "MAGDALENA" = "Magdalena",
  "NORTE_DE_SANTANDER" = "Norte de Santander",
  "SANTANDER" = "Santander"
)

# STUDYING MAXIMUM VALUES BEHAVIOR
graph_extreme_values <- data_sisaire %>%
  select(department, parameter, date, hour, value) %>%
  mutate(month = month(date), year = year(date)) %>%
  group_by(department, parameter, year) %>%
  mutate(
    threshold = quantile(value, .98),
    is_rare = (value > threshold),
    n_data = n()
  ) %>%
  group_by(department, parameter, year, month) %>%
  reframe(rare_proportion = sum(is_rare)) %>%
  unique() %>%
  mutate(department_label = recode(department, !!!department_labels)) %>%
  mutate(department_label = factor(department_label, levels = department_labels)) %>%
  ggplot(aes(x = rare_proportion)) +
  geom_density_ridges(
    aes(y = factor(month), fill = department_label),
    alpha = .5,
    scale = 1.05,
    show.legend = TRUE
  ) +
  labs(x = "N. observations over 98th percentile", y = "Month") +
  # scale_x_log10() +
  scale_fill_manual(values = RColorBrewer::brewer.pal(9, "Paired"), drop = FALSE) +
  coord_cartesian(expand = 0) +
  facet_wrap(~ parameter,
             labeller = labeller(parameter = as_labeller(
               c(PM10 = "PM[10]", PM2.5 = "PM[2.5]", O3 = "O[3]"), label_parsed
             )),
             scales = "free_x") +
  theme_bw() +
  theme(
    text = element_text(family = "serif"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    strip.background = element_blank()
  )

# Save ridgeline figure in multiple formats at high resolution
ggsave(
  "graphs/hourly/potential/monthly_extreme_values.jpg",
  plot = graph_extreme_values,
  width = 8,
  height = 6,
  dpi = 500
)
ggsave(
  "graphs/hourly/potential/monthly_extreme_values.tiff",
  plot = graph_extreme_values,
  width = 8,
  height = 6,
  dpi = 500
)
