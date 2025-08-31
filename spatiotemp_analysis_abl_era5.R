# ==========================================================
# ABLH–Air Quality Coupling: Hourly & Diurnal Analyses (R)
# ----------------------------------------------------------
# Purpose: Explore relationships between boundary layer height (ABLH),
#          wind speed (for ventilation coefficient, VC), inversion metrics,
#          and air quality (SISAIRE: PM10, PM2.5, O3) at hourly/daily scales
#          across multiple Colombian departments.
#
# What this script does (high level):
#   1) Loads packages and raw/processed data sets (ABLH, SISAIRE, wind,
#      inversion/temperature profiles).
#   2) Checks (weak/strict) stationarity of ABLH via ADF tests overall and
#      per-day. Assesses normality with Lilliefors (Kolmogorov–Smirnov) tests.
#   3) Visualizes hourly distributions and diurnal profiles by department.
#   4) Builds a time-block grouping of hours via k-means on daily ABLH vectors
#      and compares groups (Wilcoxon/Hodges–Lehmann) given non-normality.
#   5) Computes Kendall rank correlations between pollutants and atmospheric
#      variables (ABLH, VC, ABLH lags, inversion metrics) per station.
#   6) Summarizes and plots correlation directions/significance proportions.
#   7) Produces seasonal/within-year mean cycles and density/diurnal profiles.
#
# Notes:
#   - This script assumes the presence of preprocessed files under data/processed/.
#   - Random elements: k-means uses multiple starts; set a seed externally if you
#     need fully reproducible clusters (e.g., set.seed(123)).
#   - Parallelization uses furrr w/ multisession; adjust workers to your machine.
# ==========================================================

require(pacman)
# pacman::p_load will install (if missing) and load packages in one call
p_load(
  tidyverse,      # data manipulation & plotting
  data.table,     # dealing with large datasets
  lubridate,      # date/time helpers
  ggridges,       # ridgeline density plots
  scales,         # scales for axes and labels
  tseries,        # time-series tests (ADF)
  cluster,        # clustering utils (silhouette)
  nortest,        # normality tests (Lilliefors)
  broom,          # tidier test outputs
  rsample,        
  viridis,        
  yardstick,      
  RColorBrewer,   # color palettes for ggplot
  gt, gtsummary,  
  beepr,          # audio beep to signal long tasks complete
  parallel,       # base parallelization
  furrr           # future + purrr for parallel map
)

# ---------------------------
# 1) IMPORT ALL DATA FILES
# ---------------------------

# ERA5 ABLH (hourly), semicolon-separated, comma decimal
# Ensure date parses correctly; if it includes time, we coerce to Date
data_abl <- fread(
  "data/processed/ABLH/ERA5/ablh_era5_lumped.txt",
  sep = ";",
  dec = ",",
  header = TRUE
)

# Coerce to Date (drop time component if present)
data_abl[, date := lubridate::date(date)]

# SISAIRE air quality data (hourly)
data_sisaire <- fread(
  "data/processed/SISAIRE/data_sisaire_hourly.txt",
  sep = ";",
  dec = ",",
  header = TRUE
)

data_sisaire[, date := lubridate::date(date)]

# Copernicus Atmospheric Data Service (wind) - already CSV
data_winds <- fread("data/processed/WIND/wind_hourly_by_department_v1.csv")

data_winds[, date := lubridate::date(date)]

# Temperature/inversion metrics (precomputed RDS)
# all_inversion_metrics_XGB.rds: inversion base height, gradient, etc.
# all_vertical_profiles_XGB.rds: vertical profiles (not used directly below)
data_inversion   <- read_rds("data/processed/TEMP/all_inversion_metrics_XGB.rds")
data_temperature <- read_rds("data/processed/TEMP/all_vertical_profiles_XGB.rds")

# ---------------------------------------------------
# 2) FIRST ANALYSIS OF ABLH: STATIONARITY & VARIATION
# ---------------------------------------------------
# Question: How does ABLH behave across departments in time?
# We test for (weak) stationarity via ADF (null: unit root/non-stationary).
# adf.test(..., alternative = "stationary"): small p-value => reject non-stationarity.

data_abl %>%
  arrange(department, date, hour) %>%
  nest(data = -department) %>%
  mutate(adf_p = map_dbl(data, ~ adf.test(
    .x$ablh_m, alternative = "stationary", k = 0
  )[["p.value"]]))

# adf_p < 0.01 across all depts, then ABLH’s overall hourly series are stationary
# over the 5 years.

# Next, check stationarity on DAILY segments (per department, per date).
# This assesses intra-day stationarity of the hourly profile for each day.

data_abl %>%
  arrange(department, date, hour) %>%
  nest(data = -c(department, date)) %>%
  mutate(is_stationary = map_lgl(data, ~ {
    p_value <- adf.test(.x$ablh_m, alternative = "stationary", k = 0)[["p.value"]]
    stationary_95 <- p_value < 0.05
    return(stationary_95)
  })) %>%
  # Approximate denominator by year length for each date (leap-year aware)
  mutate(n_days = ifelse(year(date) %% 4 == 0, 366, 365)) %>%
  group_by(department) %>%
  reframe(n_stationary = round(100 * sum(is_stationary) / n_days, 2))

# Daily ABLH profiles tend to be non-stationary (changing statistics within a day),
# while the long-run series across many days/years are be stationary by ADF.

# ------------------------------
# 3) HOURLY DISTRIBUTIONS BY DEP
# ------------------------------
# Boxplots of ABLH by hour, faceted by department.
unique_departments <- unique(data_abl$department)
unique_departments_names <- c(
  "Antioquia",
  "Bogotá",
  "Bolívar",
  "Boyacá",
  "Caldas",
  "Cesar",
  "Magdalena",
  "Norte de Santander",
  "Santander"
)
names(unique_departments_names) <- unique_departments

graph_boxplot <- data_abl %>%
  ggplot() +
  geom_boxplot(
    aes(
      x = hour,
      y = ablh_m,
      group = interaction(hour, department)
    ),
    outlier.size = 1,
    outlier.alpha = .3
  ) +
  labs(x = "Hour", y = "ABLH (m)") +
  scale_x_continuous(breaks = seq(0, 23, 3)) +
  facet_wrap(~ department, labeller = labeller(department = unique_departments_names)) +
  theme_minimal() +
  theme(
    text = element_text(family = "serif"),
    axis.ticks = element_line(),
    panel.background = element_rect(fill = NULL),
    panel.grid = element_blank(),
    panel.spacing.y = unit(0, "pt"),
    panel.spacing.x = unit(4, "pt")
  )

# Export plots in multiple formats (publication friendly)
ggsave(
  "graphs/hourly/potential/boxplots_daily_profile.jpg",
  graph_boxplot,
  width = 8,
  height = 6,
  dpi = 500
)
ggsave(
  "graphs/hourly/potential/boxplots_daily_profile.tiff",
  graph_boxplot,
  width = 8,
  height = 6,
  dpi = 500
)
ggsave(
  "graphs/hourly/potential/boxplots_daily_profile.pdf",
  graph_boxplot,
  width = 8,
  height = 6,
  dpi = 500
)

# --------------------------------------------------
# 4) HOUR BLOCKING VIA K-MEANS ON DAILY ABLH VECTORS
# --------------------------------------------------
# Build a 24-dim vector per day (one feature per hour), scale columns, and
# search k in [2..8] by silhouette to find coherent hour groups.
# Note: We cluster hours by using the cluster assignments of centers (hour2cluster).

scaled_mat <- data_abl %>%
  pivot_wider(
    id_cols     = c(date, department),
    names_from  = hour,
    values_from = ablh_m
  ) %>%
  arrange(date) %>%
  select(-c(department, date)) %>%
  as.matrix() %>% scale()

sil_df <- map_df(2:8, function(k) {
  km <- kmeans(
    scaled_mat,
    centers = k,
    nstart = 25,
    iter.max = 100,
    algorithm = "Lloyd"
  )
  ss <- silhouette(km$cluster, dist(scaled_mat))
  tibble(k       = k, avg_sil = mean(ss[, "sil_width"]))
})

best_k  <- sil_df$k[which.max(sil_df$avg_sil)]
final_km <- kmeans(scaled_mat, centers = best_k, nstart = 50)

sil_df %>%
  ggplot(aes(best_k, avg_sil)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = 2:8) +
  labs(x     = "Number of clusters (k)", y     = "Average silhouette width", title = "Silhouette Method with Increased iter.max & Lloyd Algorithm")

# Map cluster centers back to hours: contiguous runs form hour blocks
centers <- final_km$centers
hour2cluster <- apply(centers, 2, which.max)
hour_blocks <- split(0:23, cumsum(c(1, diff(hour2cluster) != 0)))
names(hour_blocks) <- paste0("Cluster_", unique(hour2cluster))
hour_blocks

# If best_k == 2, you often get two broad periods (e.g., pre-/post-noon)
# Create a binary daytime flag accordingly (domain-driven choice here)
data_abl_v1 <- data_abl %>%
  mutate(daytime = ifelse(hour  %in%  c(0:13, 22:23), "morning_late_night", "afternoon"))

# Given non-normality (see next section), compare ABLH across 2 groups using
# Wilcoxon rank-sum (Mann–Whitney) and report HL median difference with CI.
wilcox.test(formula = ablh_m ~ daytime,
            data = data_abl_v1,
            conf.int = TRUE)

# ------------------------------
# 5) NORMALITY CHECKS (Lilliefors)
# ------------------------------
# We test distributions for ABLH, wind_speed, and SISAIRE values.
# Non-normality motivates the choice of non-parametric statistics later.

data_abl %>%
  nest(data = -department) %>%
  mutate(ks_normality = map_dbl(data, ~ lillie.test(.x$ablh_m)[["p.value"]]))

	data_winds %>%
  nest(data = -department) %>%
  mutate(ks_normality = map_dbl(data, ~ lillie.test(.x$wind_speed)[["p.value"]]))

	data_sisaire %>%
  nest(data = -c(department, parameter)) %>%
  mutate(ks_normality = map_dbl(data, ~ lillie.test(.x$value)[["p.value"]]))
# Interpretation: p-values typically small -> reject normality.

# -------------------------------------------
# 6) BUILD MASTER DATA SET & DERIVED FEATURES
# -------------------------------------------
# Join ABLH + wind + SISAIRE + inversion metrics hourly and compute
# ventilation coefficient VC = wind_speed * ablh_m.

data <- data_abl %>%
  left_join(data_winds, by = c("department", "date", "hour")) %>%
  left_join(
    data_sisaire %>% pivot_wider(names_from = parameter, values_from = value),
    by = c("department", "date", "hour")
  ) %>%
  left_join(data_inversion, by = c("department", "date", "hour"))

	data <- data %>% mutate(vc = wind_speed * ablh_m)

# Non-parametric plan moving forward due to non-normality
n_stations <- data_sisaire %>% distinct(department, station) %>%
  group_by(department) %>% reframe(n_station = n())

# -------------------------------------------------
# 7) Kendall Rank Correlations per Station/Variable
# -------------------------------------------------
# We compute Kendall's tau between pollutants and atmospheric variables
# (including ABLH lags and inversion metrics) per department × station.

setDT(data)

data_dt <- data[, .(
  department,
  station,
  date,
  hour,
  ablh_m,
  PM10,
  PM2.5,
  O3,
  vc,
  mean_gradient,
  base_inversion_m,
  grad_first_inversion
)]

setkey(data_dt, department, station, date, hour)

# Create ABLH lagged features (1, 2, 3, 6 hours) within each department×station
data_dt[, c("ablh_lag1", "ablh_lag2", "ablh_lag3", "ablh_lag6") := .(shift(ablh_m, 1L),
                                                                     shift(ablh_m, 2L),
                                                                     shift(ablh_m, 3L),
                                                                     shift(ablh_m, 6L)), by = .(department, station)]

# Reshape (long) pollutants
data_dt <- data.table::melt(
  data_dt,
  id.vars      = c(
    "department",
    "station",
    "date",
    "hour",
    "ablh_m",
    "ablh_lag1",
    "ablh_lag2",
    "ablh_lag3",
    "ablh_lag6",
    "vc",
    "mean_gradient",
    "base_inversion_m",
    "grad_first_inversion"
  ),
  measure.vars = c("PM10", "PM2.5", "O3"),
  variable.name = "pollutant",
  value.name   = "pollutant_value"
)

# Reshape (long) atmospheric variables
data_dt <- data.table::melt(
  data_dt,
  id.vars      = c(
    "department",
    "station",
    "date",
    "hour",
    "pollutant",
    "pollutant_value"
  ),
  measure.vars = c(
    "ablh_m",
    "ablh_lag1",
    "ablh_lag2",
    "ablh_lag3",
    "ablh_lag6",
    "vc",
    "mean_gradient",
    "base_inversion_m",
    "grad_first_inversion"
  ),
  variable.name = "atmospheric_variable",
  value.name   = "atmospheric_value"
)

# Drop rows with missing values produced by lags/joins
data_dt <- drop_na(data_dt)

# Parallel setup for future_map calls - adjust workers to your CPU
plan(multisession, workers = 8)

# Compute Kendall’s tau and p-value for each group
# (skip tiny groups to avoid test failures)

data_kendall_correlations <- data_dt %>%
  group_by(department, station, pollutant, atmospheric_variable) %>%
  nest() %>%
  mutate(# future_map over each nested tibble
    result = future_map(data, \(df) {
      # Small groups provide unstable estimates; skip if < 3
      if (nrow(df) < 3)
        return(tibble(kendall_corr = NA_real_, kendall_p = NA_real_))
      t <- cor.test(df$pollutant_value, df$atmospheric_value, method = "kendall")
      tibble(kendall_corr = unname(t$estimate),
             kendall_p = t$p.value)
    })) %>%
  unnest(result) %>%
  select(-data) %>%
  ungroup()

# Audible cue when the long-running step finishes
beepr::beep(sound = 3)

# Persist raw correlation results for downstream summaries

data_kendall_correlations %>%
  write.table(
    "data/processed/statistic_summaries/hourly/sisaire_ablh_vc_kendall_correlations_raw.txt",
    sep = ";",
    dec = ",",
    row.names = FALSE
  )

# Re-load from disk (ensures a clean, portable step boundary)

data_kendall_correlations <- fread(
  "data/processed/statistic_summaries/hourly/sisaire_ablh_vc_kendall_correlations_raw.txt",
  sep = ";",
  dec = ",",
  header = TRUE
)

# --------------------------------------
# 8) SUMMARIZE CORRELATION DIRECTIONS
# --------------------------------------
# Convert Kendall tau sign to direction (direct/inverse) and flag significance.
# Then aggregate to absolute and relative frequencies by department×pollutant×var.

data_kendall_correlations_summarized <- data_kendall_correlations %>%
  mutate(
    direction = ifelse(kendall_corr > 0, "direct", "inverse"),
    sig = case_when(kendall_p < .05 ~ "sig_95", TRUE ~ "no_sig")
  ) %>%
  select(department, pollutant, atmospheric_variable, direction, sig) %>%
  group_by(department, pollutant, atmospheric_variable, direction, sig) %>%
  reframe(n = n())

write.table(
  data_kendall_correlations_summarized,
  "data/processed/statistic_summaries/hourly/sisaire_ablh_vc_kendall_correlations_absolute_freq.txt",
  sep = ";",
  dec = ",",
  row.names = FALSE
)

# Relative frequencies

data_kendall_correlations_summarized_v2 <- data_kendall_correlations_summarized %>%
  group_by(department, pollutant, atmospheric_variable) %>%
  reframe(n_involved_stations = sum(n)) %>%
  inner_join(data_kendall_correlations_summarized) %>%
  mutate(proportion = n / n_involved_stations) %>%
  select(-c(n, n_involved_stations))

write.table(
  data_kendall_correlations_summarized_v2,
  "data/processed/statistic_summaries/hourly/sisaire_ablh_vc_kendall_correlations_relative_freq.txt",
  sep = ";",
  dec = ",",
  row.names = FALSE
)

# Re-load

data_kendall_correlations_summarized_v2 <- fread(
  "data/processed/statistic_summaries/hourly/sisaire_ablh_vc_kendall_correlations_relative_freq.txt",
  sep = ";",
  dec = ",",
  header = TRUE
)

# ------------------------------
# 9) PLOT: CORRELATION PROPORTIONS
# ------------------------------
# Stacked proportions by significance (95%) × direction per department,
# faceted by pollutant and atmospheric variable (with parsed math labels).

graph_station_correlation_props <- data_kendall_correlations_summarized_v2 %>%
  mutate(
    department = factor(department, levels = rev(unique_departments)),
    atmospheric_variable_rec = recode(
      atmospheric_variable,
      ablh_m               = "ABLH",
      ablh_lag1            = "ABLH[1*h~lag]",
      ablh_lag2            = "ABLH[2*h~lag]",
      ablh_lag3            = "ABLH[3*h~lag]",
      ablh_lag6            = "ABLH[6*h~lag]",
      vc                   = "VC",
      mean_gradient        = "AG",
      base_inversion_m     = "IBH",
      grad_first_inversion = "IG"
    ),
    pollutant_rec = recode(
      pollutant,
      O3 = "O[3]",
      PM10 = "PM[10]",
      PM2.5 = "PM[2.5]"
    ),
    direction = factor(direction, levels = c("inverse", "direct")),
    sig = factor(sig, levels = c("no_sig", "sig_95")),
    sig_dir = paste(sig, direction, sep = "-")
  ) %>%
  unique() %>%
  ggplot() +
  geom_col(
    aes(
      y = department,
      x = proportion,
      fill = sig_dir,
      group = interaction(sig, direction)
    ),
    position = position_stack(),
    color = "black",
    width = 1,
    linewidth = .5
  ) +
  facet_grid(
    cols = vars(pollutant_rec),
    rows = vars(atmospheric_variable_rec),
    scales = "free_y",
    labeller = label_parsed
  ) +
  scale_y_discrete(breaks = unique_departments, labels = unique_departments_names) +
  scale_x_continuous(
    breaks = seq(0, 1, .25),
    labels = function(.x) { paste0(.x * 100, "%") }
  ) +
  scale_fill_manual(
    breaks = c(
      "no_sig-inverse",
      "sig_95-inverse",
      "no_sig-direct",
      "sig_95-direct"
    ),
    labels = c("Inverse - NS", "Inverse - Sig", "Direct - NS", "Direct - Sig"),
    values = RColorBrewer::brewer.pal(6, "Paired")[c(5:6, 1:2)]
  ) +
  theme_bw() +
  theme(
    text               = element_text(family = "serif", color = "black"),
    axis.text.x        = element_text(size = 9, angle = 30, vjust = .8, color = "black"),
    axis.text.y        = element_text(size = 10, color = "black"),
    axis.title.x       = element_blank(),
    axis.title.y       = element_blank(),
    strip.text         = element_text(face = "bold", size = 12, color = "black"),
    strip.background   = element_blank(),
    axis.title         = element_text(size = 14),
    legend.position    = "bottom",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 11),
    panel.grid         = element_blank()
  )

# Save plot (multi-format)
ggsave(
  "graphs/hourly/potential/kendall_correlations.jpg",
  plot = graph_station_correlation_props,
  width = 10,
  height = 12,
  dpi = 500
)
ggsave(
  "graphs/hourly/potential/kendall_correlations.tiff",
  plot = graph_station_correlation_props,
  width = 10,
  height = 12,
  dpi = 500
)
ggsave(
  "graphs/hourly/potential/kendall_correlations.pdf",
  plot = graph_station_correlation_props,
  width = 10,
  height = 12,
  dpi = 500
)

# ---------------------------------------------------
# 10) WITHIN-YEAR MEAN CYCLES (SEASONAL SIGNATURES)
# ---------------------------------------------------
# We compute department-level daily means and then average them by day-of-year
# (yday) to reveal typical seasonal patterns for ABLH/VC and SISAIRE metrics.

# SISAIRE daily summaries -> yday means

data_sisaire_v1 <- data_sisaire %>%
  group_by(department, parameter, date) %>%
  reframe(
    daily_min  = min(value, na.rm = TRUE),
    daily_max  = max(value, na.rm = TRUE),
    daily_mean = mean(value, na.rm = TRUE)
  ) %>%
  mutate(y_day = yday(date)) %>%
  group_by(department, parameter, y_day) %>%
  reframe(
    mean_sisaire_min   = min(daily_min,  na.rm = TRUE),
    mean_sisaire_max   = max(daily_max,  na.rm = TRUE),
    mean_sisaire_value = mean(daily_mean, na.rm = TRUE)
  )

# ABLH/VC daily summaries -> yday means

data_abl_v2 <- data_abl %>%
  left_join(data_winds, by = c("department", "date", "hour")) %>%
  mutate(vc = wind_speed * ablh_m) %>%
  select(department, date, hour, ablh_m, vc) %>%
  pivot_longer(
    cols = -c(department, date, hour),
    names_to = "atm_parameter",
    values_to = "atm_value"
  ) %>%
  group_by(department, date, atm_parameter) %>%
  reframe(
    daily_min  = min(atm_value, na.rm = TRUE),
    daily_max  = max(atm_value, na.rm = TRUE),
    daily_mean = mean(atm_value, na.rm = TRUE)
  ) %>%
  mutate(y_day = yday(date)) %>%
  group_by(department, atm_parameter, y_day) %>%
  reframe(
    mean_atm_min   = min(daily_min,  na.rm = TRUE),
    mean_atm_max   = max(daily_max,  na.rm = TRUE),
    mean_atm_value = mean(daily_mean, na.rm = TRUE)
  )

# Combined seasonal plot (two y-axes with scaling factor for pollutants)

graph_ts <- ggplot() +
  # ABLH/VC mean seasonal curves
  geom_line(data = data_abl_v2,
            aes(
              x = y_day,
              y = mean_atm_value,
              color = atm_parameter,
              group = interaction(department, atm_parameter)
            )) +
  # SISAIRE seasonal curves scaled to align with ABLH axis
  geom_line(data = data_sisaire_v1,
            aes(
              x = y_day,
              y = mean_sisaire_value * 30,
              color = parameter,
              group = interaction(department, parameter)
            )) +
  labs(x = "Day of the year") +
  coord_cartesian(expand = 0) +
  scale_y_continuous(
    name = expression("ABLH (m) and VC (" * m ^ 2 / s * ")"),
    sec.axis = sec_axis(~ . / 30, name = expression("Concentration (" *
                                                      mu * "g/" * m ^ 3 * ")"))
  ) +
  scale_color_manual(
    breaks = c("ablh_m", "vc", "PM10", "PM2.5", "O3"),
    labels = c("ABLH", "VC", expression(PM[10]), expression(PM[2.5]), expression(O[3])),
    values = c(brewer.pal(5, "Paired"))
  ) +
  facet_grid(
    department ~ parameter,
    scales = "free_y",
    labeller = labeller(
      department = unique_departments_names,
      parameter = as_labeller(c(
        PM10 = "PM[10]", PM2.5 = "PM[2.5]", O3 = "O[3]"
      ), label_parsed)
    )
  ) +
  theme_bw() +
  theme(
    text              = element_text(family = "serif"),
    legend.background = element_blank(),
    legend.position   = "bottom",
    legend.title      = element_blank(),
    legend.text       = element_text(size = 15),
    panel.spacing     = unit(3, "pt"),
    strip.background  = element_blank(),
    strip.text        = element_text(size = 12),
    axis.text         = element_text(size = 12),
    axis.title        = element_text(size = 15),
    axis.line.y.right = element_line()
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.2)))

# Save seasonal plot
ggsave(
  "graphs/hourly/potential/ts_ablh_vc_sisaire.jpg",
  plot = graph_ts,
  width = 10,
  height = 14,
  dpi = 500
)
ggsave(
  "graphs/hourly/potential/ts_ablh_vc_sisaire.tiff",
  plot = graph_ts,
  width = 10,
  height = 14,
  dpi = 500
)
ggsave(
  "graphs/hourly/potential/ts_ablh_vc_sisaire.pdf",
  plot = graph_ts,
  width = 10,
  height = 14,
  dpi = 500
)

# ----------------------------------------------------
# 11) RIDGELINE DENSITY (HOURLY) WITH DUAL X-AXES
# ----------------------------------------------------
# Build a unified long data for atmospheric (ABLH/VC) and pollutants, and then
# produce ridgeline densities per hour. A scaling factor aligns units visually.

data_abl_v3 <- data_abl %>%
  left_join(data_winds, by = c("department", "date", "hour")) %>%
  mutate(vc = wind_speed * ablh_m) %>%
  select(department, date, hour, ablh_m, vc) %>%
  pivot_longer(
    cols = -c(department, date, hour),
    names_to = "parameter",
    values_to = "mean_value"
  )

# SISAIRE as-is for joining
data_sisaire_v2 <- data_sisaire %>%
  select(department, date, hour, parameter, mean_value = value)

# Combine and create ordered hour factor
data_v2 <- rbind(data_abl_v3, data_sisaire_v2) %>%
  mutate(hour_day = factor(as.character(hour), 0:23))

# Scale only ABLH/VC to share axis with pollutants in the plot
scale_ab      <- 80
inv_scale_ab  <- function(x) { x * scale_ab }

data_v2_scaled <- data_v2 %>%
  mutate(plot_value = if_else(
      parameter %in% c("ablh_m", "vc"),
      mean_value / scale_ab,
      mean_value
    ))

# Ridgeline density plot per hour, faceted by department

graph_hourly_density <- ggplot(data_v2_scaled, aes(x = plot_value, y = hour_day, fill = parameter)) +
  geom_density_ridges(
    alpha = 0.6,
    scale = 1.5,
    rel_min_height = 0.01,
    linewidth = .5
  ) +
  scale_x_continuous(
    limits = c(0, 100),
    name    = expression("Concentration (" * mu * "g/" * m ^ 3 * ")"),
    sec.axis = sec_axis(~ inv_scale_ab(.), name = "ABLH or VC (m or m·m/s)")
  ) +
  scale_y_discrete(breaks = seq(0, 23, 3)) +
  scale_fill_manual(
    breaks = c("ablh_m", "vc", "PM10", "PM2.5", "O3"),
    labels = c("ABLH", "VC", expression(PM[10]), expression(PM[2.5]), expression(O[3])),
    values = brewer.pal(5, "Paired")
  ) +
  facet_wrap(~ department, labeller = labeller(department = unique_departments_names)) +
  theme_bw(base_family = "serif") +
  theme(
    legend.position        = "inside",
    legend.position.inside = c(.95, .43),
    legend.title           = element_blank(),
    legend.background      = element_rect(fill = "gray90", color = "black"),
    panel.spacing          = unit(3, "pt"),
    strip.background       = element_blank(),
    axis.text.y            = element_text(size = 10),
    axis.line.x.top        = element_line(),
    strip.text             = element_text(size = 12)
  )

# Save ridgeline plot
ggsave(
  "graphs/hourly/potential/hourly_density_dualaxis.jpg",
  graph_hourly_density,
  width  = 10,
  height = 10,
  dpi = 500
)
ggsave(
  "graphs/hourly/potential/hourly_density_dualaxis.tiff",
  graph_hourly_density,
  width  = 10,
  height = 10,
  dpi = 500
)
ggsave(
  "graphs/hourly/potential/hourly_density_dualaxis.pdf",
  graph_hourly_density,
  width  = 10,
  height = 10,
  dpi = 500
)

# -------------------------------------------
# 12) DIURNAL PROFILES (QUANTILES & SMOOTHED)
# -------------------------------------------
# Compute hourly summaries across dates per department for pollutants (median,
# IQR) and combine with atmospheric variables (ABLH/VC) for dual-axis plots.

# Pollutants: hourly central tendency and spread

data_sisaire_v3 <- data_sisaire_v2 %>%
  group_by(department, hour, parameter) %>%
  reframe(
    value        = mean(mean_value, na.rm = TRUE),
    value_median = median(mean_value, na.rm = TRUE),
    value_q1     = quantile(mean_value, .25, na.rm = TRUE),
    value_q3     = quantile(mean_value, .75, na.rm = TRUE),
    value_sd     = sd(mean_value, na.rm = TRUE)
  )

# Diurnal ribbon (IQR) + median lines for ABLH & VC, and scaled pollutants

graph_diurnal_profile <- data_abl_v3 %>%
  rename(atm_variable = "parameter", atm_value = "mean_value") %>%
  group_by(department, hour, atm_variable) %>%
  reframe(
    atm_value_mean   = mean(atm_value, na.rm = TRUE),
    atm_value_sd     = sd(atm_value, na.rm = TRUE),
    atm_value_median = median(atm_value, na.rm = TRUE),
    atm_value_q1     = quantile(atm_value, .25, na.rm = TRUE),
    atm_value_q3     = quantile(atm_value, .75, na.rm = TRUE)
  ) %>%
  ggplot(aes(x = hour)) +
  geom_ribbon(
    aes(ymin = atm_value_q1, ymax = atm_value_q3, fill = atm_variable),
    alpha = .5, show.legend = FALSE
  ) +
  geom_line(aes(y = atm_value_median, color = atm_variable), linewidth = 1.2) +
  geom_ribbon(
    data = data_sisaire_v3,
    aes(ymin = value_q1 * 20, ymax = value_q3 * 20, fill = parameter),
    alpha = .5, show.legend = FALSE
  ) +
  geom_line(data = data_sisaire_v3,
            aes(y = value_median * 20, color = parameter),
            linewidth = 1.2) +
  labs(x = "Hour of the day") +
  scale_y_continuous(
    name = expression("ABLH (m) and VC (" * m ^ 2 / s * ")"),
    sec.axis = sec_axis(~ . / 20, name = expression("SCP Concentration (" *
                                                      mu * "g/" * m ^ 3 * ")"))
  ) +
  scale_color_manual(
    breaks = c("ablh_m", "vc", "PM10", "PM2.5", "O3"),
    labels = c("ABLH", "VC", expression(PM[10]), expression(PM[2.5]), expression(O[3])),
    values = RColorBrewer::brewer.pal(5, "Paired")
  ) +
  scale_fill_manual(
    breaks = c("ablh_m", "vc", "PM10", "PM2.5", "O3"),
    labels = c("ABLH", "VC", expression(PM[10]), expression(PM[2.5]), expression(O[3])),
    values = RColorBrewer::brewer.pal(5, "Paired")
  ) +
  facet_grid(
    department ~ parameter,
    labeller = labeller(
      department = unique_departments_names,
      parameter = as_labeller(c(
        PM10 = "PM[10]", PM2.5 = "PM[2.5]", O3 = "O[3]"
      ), label_parsed)
    ),
    scales = "free_y"
  ) +
  theme_bw() +
  theme(
    text              = element_text(family = "serif"),
    axis.title        = element_text(size = 15),
    axis.line.y.right = element_line(),
    axis.text         = element_text(size = 12, color = "black"),
    strip.background  = element_blank(),
    strip.text        = element_text(size = 12),
    legend.position   = "bottom",
    legend.title      = element_blank(),
    legend.text       = element_text(size = 15)
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 2)))

# Save diurnal quantile plot
ggsave(
  "graphs/hourly/potential/ablh_pollutants_diurnal_profile.jpg",
  graph_diurnal_profile,
  width = 10,
  height = 14,
  dpi = 500
)
ggsave(
  "graphs/hourly/potential/ablh_pollutants_diurnal_profile.pdf",
  graph_diurnal_profile,
  width = 10,
  height = 14,
  dpi = 500
)
ggsave(
  "graphs/hourly/potential/ablh_pollutants_diurnal_profile.tiff",
  graph_diurnal_profile,
  width = 10,
  height = 14,
  dpi = 500
)

# Smoothed diurnal profiles (geom_smooth) for an overall trend view

data_sisaire_v4 <- data_sisaire_v2 %>%
  group_by(department, date, hour, parameter) %>%
  reframe(
    value   = mean(mean_value, na.rm = TRUE),
    value_sd = sd(mean_value, na.rm = TRUE)
  )

graph_diurnal_profile_smoothed <- data_abl_v3 %>%
  rename(atm_variable = "parameter", atm_value = "mean_value") %>%
  ggplot() +
  geom_smooth(aes(x = hour, y = atm_value, color = atm_variable, fill = atm_variable),
              se = FALSE) +
  geom_smooth(data = data_sisaire_v4,
              aes(x = hour, y = value * 20, color = parameter, fill = parameter),
              se = FALSE) +
  labs(x = "Hour of the day") +
  scale_y_continuous(
    name = expression("ABLH (m) and VC (" * m ^ 2 / s * ")"),
    sec.axis = sec_axis(~ . / 20, name = expression("Concentration (" *
                                                      mu * "g/" * m ^ 3 * ")"))
  ) +
  scale_color_brewer(
    breaks = c("ablh_m", "vc", "PM10", "PM2.5", "O3"),
    labels = c("ABLH", "VC", expression(PM[10]), expression(PM[2.5]), expression(O[3])),
    palette = "Paired"
  ) +
  scale_fill_brewer(
    breaks = c("ablh_m", "vc", "PM10", "PM2.5", "O3"),
    labels = c("ABLH", "VC", expression(PM[10]), expression(PM[2.5]), expression(O[3])),
    palette = "Paired"
  ) +
  facet_grid(
    parameter ~ department,
    labeller = labeller(
      department = unique_departments_names,
      parameter = as_labeller(c(
        PM10 = "PM[10]", PM2.5 = "PM[2.5]", O3 = "O[3]"
      ), label_parsed)
    ),
    scales = "free_y"
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "serif"),
    axis.line.y.right = element_line(),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()
  )

# Save smoothed diurnal plot
ggsave(
  "graphs/hourly/potential/ablh_pollutants_diurnal_profile_smoothed.jpg",
  graph_diurnal_profile_smoothed,
  width = 15,
  height = 7,
  dpi = 500
)
ggsave(
  "graphs/hourly/potential/ablh_pollutants_diurnal_profile_smoothed.tiff",
  graph_diurnal_profile_smoothed,
  width = 15,
  height = 7,
  dpi = 500
)
ggsave(
  "graphs/hourly/potential/ablh_pollutants_diurnal_profile_smoothed.pdf",
  graph_diurnal_profile_smoothed,
  width = 15,
  height = 7,
  dpi = 500
)

# ---------------------------------------------------
# End of script
# ---------------------------------------------------