# =============================================================
# IDEAM Station Catalog -> Cleaned Coordinate Reference (R)
# -------------------------------------------------------------
# Purpose
#   Build a harmonized catalog of station coordinates for SISAIRE pollutants
#   (2019–present) from the raw IDEAM Excel sheet. The output table serves as
#   a stable reference to geolocate stations in downstream analyses.
#
# What this script does (overview)
#   1) Loads the raw Excel catalog and keeps rows from 2019 onward.
#   2) Filters variables to SISAIRE pollutants only (NO2, O3, PM10, PM2.5, CO, SO2).
#   3) Drops metadata columns not required for the coordinate catalog.
#   4) Restricts to a target set of departments.
#   5) Normalizes text fields (station/department/city) to a consistent, ASCII,
#      upper‑case, underscore‑separated format to avoid duplicates by spelling.
#   6) Removes “(S)” suffixed stations and keeps the most recent year per station
#      as the coordinate anchor (assumes station location is stable in latest year).
#   7) Inspects duplicates and then resolves multi‑coordinate entries by averaging
#      lat/lon within env_authority×department×city×station groups.
#   8) Systematically disambiguates generic station names (e.g., SENA, BOMBEROS,
#      HOSPITAL) using department context; removes specific stations when needed.
#   9) Writes the cleaned, harmonized coordinate table to disk as a delimited file.
# =============================================================

# Load required packages
require(pacman)
p_load(tidyverse, data.table, readxl, janitor, stringi)

# Read raw station catalog from Excel
data_station <- read_excel(
  "data/raw_ideam/estaciones_2011_2023_ideam_v0.xlsx",
  sheet = "Hoja1",
  col_names = TRUE
)

# Initial filtering
# Keep years from 2019 onward
# Keep only SISAIRE pollutant variables
# Drop columns not needed for the coordinate catalog
# Restrict to the selected set of departments
data_station_v0 <- data_station %>%
  filter(year >= 2019) %>%
  filter(variable %in% c("NO2", "O3", "PM10", "PM2.5", "CO", "SO2")) %>%
  select(-c(
    time_exposition,
    units,
    station_id,
    metropolitan_area,
    variable
  )) %>%
  filter(
    department %in% c(
      "Antioquia",
      "Bogotá DC",
      "Bolívar",
      "Boyacá",
      "Caldas",
      "Cesar",
      "Cundinamarca",
      "Magdalena",
      "Norte de Santander",
      "Santander"
    )
  )

# Utility for consistent station and locality text normalization
# Removes accents, trims spaces, removes dots, and uses underscores as separators
text_cleanning <- function(string) {
  string_clean <- string %>%
    stri_trans_general("Latin-ASCII") %>%  # accents to ASCII
    toupper() %>%                          # upper case
    gsub("\\s+$", "", .) %>%               # trim right
    gsub("^\\s+", "", .) %>%               # trim left
    gsub("\\.", "", .) %>%                 # remove dots
    gsub("\\s+", "_", .) %>%               # spaces to underscores
    gsub("\\-", "_", .) %>%                # dashes to underscores
    gsub("\\_+", "_", .)                   # collapse multiple underscores
  return(string_clean)
}

data_station_v1 <- data_station_v0 %>%
  # Remove stations tagged with suffix S which represent special or secondary entities
  filter(!(str_detect(station, "\\s\\(S\\).*$"))) %>%
  mutate(
    station = text_cleanning(station),
    department = text_cleanning(department),
    city = text_cleanning(city)
  ) %>%
  # For each station keep the most recent year entry to anchor coordinates
  group_by(department, city, station) %>%
  mutate(last_year = max(year)) %>%
  filter(year == last_year) %>%
  select(-last_year) %>%
  unique()

# Identify stations with more than one coordinate entry
# The join displays duplicates to be inspected manually
data_station_v1 %>%
  group_by(department, city, station) %>%
  reframe(count = n()) %>%
  filter(count > 1) %>%
  inner_join(data_station_v1)

# Resolve duplicated coordinates and harmonize specific station names
data_station_v2 <- data_station_v1 %>%
  group_by(env_authority, department, city, station) %>%
  # Average latitude and longitude when multiple values exist
  mutate(latitude = mean(latitude),
         longitude = mean(longitude)) %>%
  select(-year) %>%
  unique() %>%
  ungroup() %>%
  # Disambiguate SENA stations by department
  mutate(station = ifelse(
    station == "SENA",
    ifelse(department == "CUNDINAMARCA", "MOSQUERA_SENA", "SOGAMOSO_SENA"),
    station
  )) %>%
  # Disambiguate BOMBEROS stations by department
  mutate(station = ifelse(
    station == "BOMBEROS",
    ifelse(
      department == "BOYACA",
      "NOBSA_BOMBEROS",
      "VALLEDUPAR_BOMBEROS"
    ),
    station
  )) %>%
  # Remove stations that should not be part of the final catalog
  filter(
    !(
      station == "ALCALDIA_DE_GUARNE_2022" |
        station == "CA_BELLO" | station == "CA_GOMEZ_PLATA" |
        station == "CA_SABANETA"
    )
  ) %>%
  # Disambiguate generic HOSPITAL station names using department
  mutate(station = ifelse(
    station == "HOSPITAL",
    case_when(
      department == "BOYACA" ~ "HOSPITAL_SOGAMOSO",
      department == "CUNDINAMARCA" ~ "HOSPITAL_SOACHA",
      TRUE ~ "HOSPITAL_LA_ESTRELLA"
    ),
    station
  ))

# Persist the cleaned and harmonized station coordinate catalog
write.table(
  data_station_v2,
  "data/clean_data/ideam/station_coordinates.txt",
  sep = ";",
  dec = ",",
  row.names =