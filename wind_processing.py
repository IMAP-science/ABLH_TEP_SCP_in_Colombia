"""
ERA5 Wind Data Extraction Pipeline - FAST YEARLY VERSION
--------------------------------------------------------
- Optimizado para:
    * Reusar datos ya descargados (archivos anuales o consolidados 2020–2024)
    * Descargar ÚNICAMENTE años faltantes, 1 petición por año y área
    * Consolidar en un archivo por área: era5_<Area>_2020-2024.nc
- Control de descargas:
    * ALLOW_DOWNLOADS = False -> NO descarga nada, solo usa datos existentes.
    * ALLOW_DOWNLOADS = True  -> descarga años faltantes.
"""

import warnings
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Set
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp

import cdsapi
import geopandas as gpd
import xarray as xr
import pandas as pd
import numpy as np

warnings.filterwarnings("ignore")

# ------------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------------
BASE_DIR = Path(r"D:\Luismi\I-ABL")
BOUNDARIES_FILE = BASE_DIR / "data" / "boundaries" / "boundaries_study_area_col_v2.gpkg"
RAW_WIND_DIR = BASE_DIR / "data" / "clean_data" / "wind"
PROC_WIND_DIR = BASE_DIR / "data" / "processed" / "wind"

RAW_WIND_DIR.mkdir(parents=True, exist_ok=True)
PROC_WIND_DIR.mkdir(parents=True, exist_ok=True)

N_CORES = 8
YEARS = ["2020", "2021", "2022", "2023", "2024"]

# Si está en False, NO se descarga nada aunque falten años
ALLOW_DOWNLOADS = False

# Nombre de la columna en el GPKG con el nombre del departamento
DEPT_FIELD = "department"


# ------------------------------------------------------------------
# HELPER FUNCTIONS
# ------------------------------------------------------------------

def validate_nc_file(nc_file: Path) -> bool:
    """Validate NetCDF file exists and contains required ERA5 wind variables."""
    if not nc_file.exists():
        return False
    try:
        with xr.open_dataset(nc_file) as ds:
            return "u10" in ds.variables and "v10" in ds.variables
    except Exception:
        return False


def get_time_coordinate(ds: xr.Dataset) -> str:
    """Get the name of the time coordinate."""
    return "valid_time" if "valid_time" in ds.coords else "time"


def check_csv_exists(area_name: str, year: str) -> bool:
    """Check if processed CSV already exists and is valid."""
    output_file = PROC_WIND_DIR / f"{area_name}_wind_{year}.csv"
    if not output_file.exists():
        return False
    try:
        df = pd.read_csv(output_file, nrows=5)
        return "datetime" in df.columns and len(df) > 0
    except Exception:
        return False


# ------------------------------------------------------------------
# SMART INVENTORY
# ------------------------------------------------------------------

def find_all_nc_files(area_name: str) -> Dict[str, str]:
    """
    Find ALL NetCDF files for an area, including:
    - era5_<Area>_2020-2024.nc (consolidated)
    - era5_<Area>_2020.nc, ..., era5_<Area>_2024.nc (individual years)
    """
    found_files: Dict[str, str] = {}

    # Check for consolidated file first
    consolidated = RAW_WIND_DIR / f"era5_{area_name}_2020-2024.nc"
    if validate_nc_file(consolidated):
        found_files["consolidated"] = str(consolidated)
        return found_files

    # Check for individual year files
    for year in YEARS:
        year_file = RAW_WIND_DIR / f"era5_{area_name}_{year}.nc"
        if validate_nc_file(year_file):
            found_files[year] = str(year_file)

    return found_files


def identify_missing_data(area_name: str) -> Tuple[str, Optional[str], Set[str]]:
    """
    Identify what data is missing for an area.

    Returns:
        status: 'ready', 'needs_consolidation', or 'needs_download'
        consolidated_file_path: path if ready, otherwise None
        missing_years: set of years that need downloading (if any)
    """
    files = find_all_nc_files(area_name)

    # Best case: have consolidated file
    if "consolidated" in files:
        return "ready", files["consolidated"], set()

    # Have some (or all) individual years
    years_found = {y for y in files.keys() if y in YEARS}
    missing_years = set(YEARS) - years_found

    if not missing_years and years_found:
        # all years exist individually, only consolidation is missing
        return "needs_consolidation", None, set()

    # Some years missing
    return "needs_download", None, missing_years


# ------------------------------------------------------------------
# FAST CONSOLIDATION
# ------------------------------------------------------------------

def consolidate_years_fast(area_name: str, year_files: Dict[str, str]) -> Optional[str]:
    """Quickly consolidate year files without unnecessary operations."""
    full_file = RAW_WIND_DIR / f"era5_{area_name}_2020-2024.nc"

    if full_file.exists() and validate_nc_file(full_file):
        return str(full_file)

    try:
        print(f"  📦 [{area_name}] Consolidando años disponibles...")

        datasets = [
            xr.open_dataset(year_files[year])
            for year in sorted(YEARS)
            if year in year_files
        ]
        if not datasets:
            print(f"  ✗ [{area_name}] No hay datasets anuales para consolidar.")
            return None

        time_coord = get_time_coordinate(datasets[0])
        combined = xr.concat(datasets, dim=time_coord)
        combined.to_netcdf(full_file)

        for ds in datasets:
            ds.close()
        combined.close()

        print(f"  ✓ [{area_name}] Consolidado en {full_file.name}")
        return str(full_file)

    except Exception as e:
        print(f"  ✗ [{area_name}] Error consolidando: {e}")
        return None


# ------------------------------------------------------------------
# YEARLY DOWNLOAD (OPTIONAL)
# ------------------------------------------------------------------

def download_year_single(client: cdsapi.Client, area_info: Dict, year: str) -> bool:
    """
    Download one full year for a given area in a single CDS request (only if ALLOW_DOWNLOADS=True).

    - For 2024, only months up to current month are requested.
    - For 2020–2023, all 12 months are included.
    """
    area_name = area_info["name"]
    year_file = RAW_WIND_DIR / f"era5_{area_name}_{year}.nc"

    if validate_nc_file(year_file):
        print(f"  ✓ [{area_name} {year}] Año ya descargado (archivo válido).")
        return True

    # If downloads are not allowed, skip
    if not ALLOW_DOWNLOADS:
        print(f"  ✗ [{area_name} {year}] Falta año y ALLOW_DOWNLOADS=False. No se descarga.")
        return False

    # Define months to request
    if year == "2024":
        current_month = datetime.now().month
        months = [f"{m:02d}" for m in range(1, current_month + 1)]
    else:
        months = [f"{m:02d}" for m in range(1, 13)]

    print(f"  ⬇ [{area_name} {year}] Descargando año completo (1 petición)...")

    try:
        client.retrieve(
            "reanalysis-era5-single-levels",
            {
                "product_type": "reanalysis",
                "variable": ["10m_u_component_of_wind", "10m_v_component_of_wind"],
                "year": year,
                "month": months,
                "day": [f"{d:02d}" for d in range(1, 32)],
                "time": [f"{h:02d}:00" for h in range(24)],
                "area": area_info["bbox"],
                "format": "netcdf",
            },
            str(year_file),
        )
        print(f"  ✓ [{area_name} {year}] Año completo descargado.")
        return validate_nc_file(year_file)
    except Exception as e:
        print(f"  ✗ [{area_name} {year}] Error en descarga anual: {e}")
        return False


# ------------------------------------------------------------------
# SMART ORCHESTRATION
# ------------------------------------------------------------------

def prepare_area_data(area_info: Dict) -> Tuple[str, bool, Optional[str], Set[str]]:
    """
    Prepare data for an area (consolidate or download as needed).

    Returns:
        (area_name, success, consolidated_file_path, missing_years)
    """
    area_name = area_info["name"]
    status, consolidated_file, missing_years = identify_missing_data(area_name)

    print(f"\n{'=' * 70}")
    print(f"[{area_name}]")
    print(f"{'=' * 70}")

    # Case 1: Already ready
    if status == "ready":
        print(f"✓ Datos completos: {Path(consolidated_file).name}")
        return area_name, True, consolidated_file, set()

    # Case 2: Need consolidation only
    if status == "needs_consolidation":
        print("⚡ Consolidando archivos existentes...")
        files = find_all_nc_files(area_name)
        consolidated = consolidate_years_fast(area_name, files)
        if consolidated:
            return area_name, True, consolidated, set()
        return area_name, False, None, set()

    # Case 3: Need to download
    print(f"⬇ Años faltantes: {sorted(missing_years)}")

    if not ALLOW_DOWNLOADS:
        print(f"✗ [{area_name}] ALLOW_DOWNLOADS=False → no se descargan años faltantes.")
        return area_name, False, None, missing_years

    try:
        client = cdsapi.Client()
        files = find_all_nc_files(area_name)

        # Download missing years
        for year in sorted(missing_years):
            if not download_year_single(client, area_info, year):
                return area_name, False, None, missing_years

            year_file = RAW_WIND_DIR / f"era5_{area_name}_{year}.nc"
            files[year] = str(year_file)

        # Consolidate all available year files into a single 2020–2024 file
        consolidated = consolidate_years_fast(area_name, files)
        if consolidated:
            print(f"✓ [{area_name}] Completo")
            return area_name, True, consolidated, set()

        return area_name, False, None, missing_years

    except Exception as e:
        print(f"✗ [{area_name}] Error: {e}")
        return area_name, False, None, missing_years


# ------------------------------------------------------------------
# FAST PROCESSING (Skip what's already done)
# ------------------------------------------------------------------

def process_year_worker(args: Tuple[str, str, str]) -> Dict:
    """Process one year (skip if CSV exists)."""
    area_name, year, nc_file = args
    output_file = PROC_WIND_DIR / f"{area_name}_wind_{year}.csv"

    # Skip if already processed
    if check_csv_exists(area_name, year):
        return {
            "area": area_name,
            "year": year,
            "success": True,
            "file": str(output_file),
            "cached": True,
        }

    try:
        # Calculate expected hours
        if year == "2024":
            expected_hours = pd.date_range(f"{year}-01-01", datetime.now(), freq="H").shape[0]
        else:
            days = 366 if int(year) % 4 == 0 else 365
            expected_hours = days * 24

        # Load and process
        with xr.open_dataset(nc_file) as ds:
            time_coord = get_time_coordinate(ds)
            ds_year = ds.sel({time_coord: slice(f"{year}-01-01", f"{year}-12-31")})

            if len(ds_year[time_coord]) == 0:
                return {
                    "area": area_name,
                    "year": year,
                    "success": False,
                    "error": "No data",
                }

            u10 = ds_year["u10"].mean(dim=["latitude", "longitude"])
            v10 = ds_year["v10"].mean(dim=["latitude", "longitude"])

            wind_speed = np.sqrt(u10**2 + v10**2)
            wind_direction = (np.arctan2(-u10, -v10) * 180 / np.pi) % 360

            df = pd.DataFrame(
                {
                    "datetime": pd.to_datetime(wind_speed[time_coord].values),
                    "wind_speed_ms": wind_speed.values,
                    "wind_direction_deg": wind_direction.values,
                    "u10": u10.values,
                    "v10": v10.values,
                }
            )

        df = df.sort_values("datetime").reset_index(drop=True)
        df.to_csv(output_file, index=False)

        obtained = len(df)
        percentage = (obtained / expected_hours) * 100

        print(f"  ✓ [{area_name} {year}] {obtained}/{expected_hours} ({percentage:.1f}%)")

        return {
            "area": area_name,
            "year": year,
            "success": True,
            "expected": expected_hours,
            "obtained": obtained,
            "percentage": percentage,
            "file": str(output_file),
        }

    except Exception as e:
        print(f"  ✗ [{area_name} {year}] Error: {e}")
        return {"area": area_name, "year": year, "success": False, "error": str(e)}


# ------------------------------------------------------------------
# MAIN EXECUTION
# ------------------------------------------------------------------

def main():
    """Main execution - optimized for speed, with optional downloads."""
    print("=" * 70)
    print("ERA5 WIND DATA - FAST PROCESSING MODE")
    print("=" * 70)
    print("⚡ Skipping existing data | Downloading only what's needed")
    print(f"Cores: {N_CORES}")
    print(f"ALLOW_DOWNLOADS: {ALLOW_DOWNLOADS}")

    # Step 1: Load study areas and build Area -> Department mapping
    print("\n[1/4] Cargando áreas de estudio...")
    gdf = gpd.read_file(BOUNDARIES_FILE)

    if DEPT_FIELD not in gdf.columns:
        print(f"⚠ ADVERTENCIA: la columna '{DEPT_FIELD}' no existe en el GPKG.")
        print("   Usaré 'Area_0', 'Area_1', ... como identificadores,")
        print("   y no podré asociar nombre de departamento explícito.")

    study_areas: List[Dict] = []
    area_department_map = []      # para CSV de mapeo
    global area_dept_lookup
    area_dept_lookup = {}         # dict: "Area_0" -> "Nombre Departamento"

    for idx, row in gdf.iterrows():
        bounds = row.geometry.bounds

        # mantén los nombres coherentes con tus NetCDF: Area_0, Area_1, ...
        area_name = row.get("name", f"Area_{idx}")

        if DEPT_FIELD in gdf.columns:
            dept_name = row[DEPT_FIELD]
        else:
            dept_name = None

        study_areas.append(
            {
                "id": idx,
                "name": area_name,
                "bounds": bounds,
                "bbox": [bounds[3], bounds[0], bounds[1], bounds[2]],
            }
        )

        dept_label = str(dept_name) if dept_name not in [None, ""] else area_name
        area_dept_lookup[area_name] = dept_label

        area_department_map.append(
            {"id": idx, "area_name": area_name, "department": dept_label}
        )

    print(f"✓ {len(study_areas)} áreas encontradas")

    # Guardar mapeo Area_X -> Departamento
    mapping_file = PROC_WIND_DIR / "area_department_mapping.csv"
    pd.DataFrame(area_department_map).to_csv(mapping_file, index=False)
    print(f"✓ Mapeo área-departamento guardado en: {mapping_file}")

    # Step 2: Quick inventory
    print("\n[2/4] Inventario rápido...")
    download_needed: List[Dict] = []
    ready_areas: Dict[str, str] = {}

    for area in study_areas:
        status, consolidated, missing = identify_missing_data(area["name"])
        if status == "ready":
            ready_areas[area["name"]] = consolidated
            print(f"  ✓ {area['name']}: Listo")
        elif status == "needs_consolidation":
            print(f"  ⚡ {area['name']}: Solo consolidar")
            download_needed.append(area)
        else:
            # needs_download
            if ALLOW_DOWNLOADS:
                print(f"  ⬇ {area['name']}: Falta {len(missing)} años → se intentará descargar")
                download_needed.append(area)
            else:
                print(f"  ✗ {area['name']}: Faltan {len(missing)} años y ALLOW_DOWNLOADS=False → se omite")

    # Step 3: Prepare data (consolidate + optional download)
    print("\n[3/4] Preparando datos (consolidación y descargas si se permiten)...")
    if download_needed:
        for area in download_needed:
            result = prepare_area_data(area)
            if result[1]:
                ready_areas[result[0]] = result[2]
    else:
        print("✓ No hay áreas que requieran consolidación o descarga.")

    if not ready_areas:
        print("⚠ No hay áreas listas para procesar (ready_areas vacío).")
        print("=" * 70)
        return

    # Step 4: Process in parallel
    print(f"\n[4/4] Procesando a CSV (paralelo, {N_CORES} cores)...")

    tasks: List[Tuple[str, str, str]] = []
    for area in study_areas:
        if area["name"] in ready_areas:
            for year in YEARS:
                tasks.append((area["name"], year, ready_areas[area["name"]]))

    print(f"  Total: {len(tasks)} tareas")

    results: List[Dict] = []
    if tasks:
        with ProcessPoolExecutor(max_workers=N_CORES) as executor:
            futures = {executor.submit(process_year_worker, task): task for task in tasks}
            for future in as_completed(futures):
                results.append(future.result())

    # Summary
    print("\n" + "=" * 70)
    print("RESUMEN")
    print("=" * 70)

    successful = sum(1 for r in results if r.get("success"))
    cached = sum(1 for r in results if r.get("cached", False))
    processed = successful - cached

    print(f"✓ Completadas: {successful}/{len(results)}")
    print(f"  - Ya existían: {cached}")
    print(f"  - Procesadas ahora: {processed}")

    # Save completeness report
    report_data = [
        {
            "Area": r["area"],
            "Year": r["year"],
            "Expected_Hours": r.get("expected", 0),
            "Obtained_Hours": r.get("obtained", 0),
            "Completeness_%": r.get("percentage", 0),
        }
        for r in results
        if r.get("success") and "expected" in r
    ]

    if report_data:
        df_report = pd.DataFrame(report_data)
        report_file = PROC_WIND_DIR / "completeness_report.csv"
        df_report.to_csv(report_file, index=False)
        avg = df_report["Completeness_%"].mean()
        print(f"Completitud promedio: {avg:.2f}%")
        print(f"Reporte: {report_file}")
    else:
        print("⚠ No se generó completeness_report (no hubo años con 'expected').")

    # Archivo maestro con todas las series
    print("\nGenerando archivo maestro con todas las series temporales...")

    all_rows = []
    for r in results:
        if not r.get("success"):
            continue
        area = r["area"]
        year = r["year"]
        csv_file = PROC_WIND_DIR / f"{area}_wind_{year}.csv"
        if not csv_file.exists():
            continue

        df_year = pd.read_csv(csv_file)
        df_year["area"] = area
        df_year["department"] = area_dept_lookup.get(area, area)
        df_year["year"] = year
        all_rows.append(df_year)

    if all_rows:
        df_all = pd.concat(all_rows, ignore_index=True)
        master_file = PROC_WIND_DIR / "wind_timeseries_all_departments.csv"
        df_all.to_csv(master_file, index=False)
        print(f"✓ Archivo maestro generado: {master_file}")
        print(f"  Filas totales: {len(df_all)}")
    else:
        print("⚠ No se encontraron series para construir el archivo maestro.")

    print("=" * 70)


if __name__ == "__main__":
    mp.freeze_support()
    main()
