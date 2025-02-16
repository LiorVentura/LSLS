# Green and Dense Cities: A Sustainable Path Towards a Nature-Positive Future

This repository contains the R scripts and relevant functions used in the manuscript **"Green and dense cities: A sustainable path towards nature-positive future"** (in review, *Nature Cities*). The code implements urban growth scenarios, simulates land-cover changes at 100 m resolution, and builds species distribution models (SDMs) to evaluate impacts on bird diversity. It also estimates human access to nature under various growth and restoration scenarios.

## Table of Contents
1. [Overview](#overview)
2. [Repository Structure](#repository-structure)
3. [Data Requirements](#data-requirements)
4. [Scripts](#scripts)
   - [1. `scenario_fnct_.R`](#1-scenario_fnct_r)
   - [2. `scenario_run_.R`](#2-scenario_run_r)
   - [3. `SDM_fnct_.R`](#3-sdm_fnct_r)
   - [4. `SDM_maxent_.R`](#4-sdm_maxent_r)
   - [5. `SDM_maxent_ebird_.R` (Optional)](#5-sdm_maxent_ebird_r-optional)
5. [Usage Examples](#usage-examples)
6. [References](#references)
7. [Contact](#contact)

---

## Overview

Our research explores how different urban growth trajectories (sprawling vs. densifying) and restoration strategies (regional vs. local land-sparing) influence:
- **Bird diversity** (estimated via Maxent species distribution models).
- **Human access to nature** (quantified as proximity and extent of green spaces).

We combine four main urban growth scenarios:

1. **Land Sharing (LSh)**: Continued urban sprawl on metropolitan edges.  
2. **Land Sparing (LSp)**: Increasing residential building heights within existing urban areas (densification).  
3. **Full Land Sparing (LSpF)**: Regional restoration of outer urban areas combined with densification in the center.  
4. **Local Land Sparing (LSpL)**: Localized patches of restored green spaces within each urban block, combined with densification of surrounding pixels.

Using these scenarios, we iteratively simulate land-cover changes from a baseline year (2020) through time steps (2025, 2030, 2035, 2040, 2045, 2050). We then project species distributions under each scenario to evaluate potential impacts on biodiversity.

---

## Repository Structure

```
repo_root/
├── scenario_fnct_.R
├── scenario_run_.R
├── SDM_fnct_.R
├── SDM_maxent_.R
├── SDM_maxent_ebird_.R
├── 100_grid_variables/  (example input raster folder)
├── tables/             (example tabular data folder)
├── Maxent/             (example output folder for models)
├── scenarios/          (example output folder for scenario runs)
└── ...                 (other folders, e.g., for GIS layers or results)
```

---

## Data Requirements

1. **Spatial Grid (100 m resolution)**  
   - A set of raster maps at 100 m resolution representing baseline (2020) land-cover and environmental variables (building heights, building volume, green spaces, roads, etc.).  
   - Example placeholders in `./100_grid_variables/`.

2. **Pixel Table of Environment (`env_table`)**  
   - A CSV or RDS table containing one row per 100 m pixel with columns for coordinates (`X`, `Y`) and environmental variables (e.g., `bldg_vol`, `bldg_cover`, `tree_cover`, etc.).  
   - Example placeholders in `./tables/env100_raster_table_.csv`.

3. **Auxiliary Tables** (for scenario constraints)  
   - Municipal / district plans defining building height limits.  
   - Protected areas polygons or attributes.  
   - Statistical area definitions for local vs. regional restoration.  
   - Example placeholders in `./tables/`.

4. **Population/Housing Projections**  
   - An external CSV specifying how much housing volume must be added by each year-step.  
   - See `pop_projection_2050.csv` in `./tables/` as an example.

5. **Bird Observations**  
   - For the Maxent models, you need occurrence points and possibly background or pseudo-absence data.  
   - The repository uses two data sources: (1) a standardized bird survey dataset (`birdsurvey_alldata_clean.csv`), (2) an eBird dataset (`ebird_2024_thinned_final.csv`).  
   - Example shapefiles in `./survey_points_new.shp` or `./eBird_points_filtered_ITM.shp`.

*All data references shown here are placeholders illustrating required data formats. Actual data used in the study are not included in this repository for privacy/licensing reasons.*

---

## Scripts

### 1. `scenario_fnct_.R`
**Purpose**: Defines all core functions for simulating urban land-cover changes according to each scenario.  
**Key Functions**:

1. **`focal_500(rast, pad_rast)`**  
   - Applies a 500 m focal (moving window) mean to a raster variable (e.g., `tree_cover`) by extending with a "padding" raster.  
   - Used to compute 500 m neighborhood variables (e.g., `tree_cover_500`).

2. **`get_changedist(env_df)`**  
   - Calculates distances to already-changed pixels in the simulation, used as weights to cluster new developments.

3. **`get_builtdist(env_df)`**  
   - Calculates distances to built pixels to help locate near-edge expansions.

4. **`get_opentdist(env_df)`**  
   - Calculates distances to "periphery" or open areas within each district.

5. **Scenario Functions**:
   - **`LShH(env_table, vol_trgt)`**: *Land Sharing (High)*  
     Sprawling development by converting outer, low-building-volume pixels to match moderate-density pixels.
   - **`LSpH(env_table, vol_trgt, non_res_trgt=0)`**: *Land Sparing (High)*  
     Densification by increasing building height/volume in existing built pixels until a housing volume target is met.  
   - **`LSpF(env_table, vol_trgt, Spare_trgt)`**: *Full Regional Land Sparing*  
     Combination of densification in city center and restoration (i.e., clearing) in outer metropolitan areas.  
   - **`LSpL(env_table, vol_trgt, spare_trgt)`**: *Local Land Sparing*  
     Restoration patches within each block, accompanied by densification of the remaining built pixels.  

6. **`scen_sum(env_list)`**  
   - Summarizes scenario outputs (total building volume, cover of roads, vegetation, etc.) across all time steps.

7. **`scen_export(env_list_i, path)`**  
   - Exports the final environment table to disk as CSV and writes new raster layers (e.g., `tree_cover`, `tree_cover_500`, etc.) to a specified path.

---

### 2. `scenario_run_.R`
**Purpose**: Main driver script to run the scenario simulations. It imports data, calls the functions from `scenario_fnct_.R`, iterates over year-steps, and saves outputs.

**Workflow**:
1. **Data Import**  
   - Reads pixel table `env100_raster_table_.csv`, building height limits, protected areas, TAD (Transit Access District) info, etc.
2. **Data Preparation**  
   - Merges auxiliary info to create a single `env_table` with building volumes, protected statuses, etc.
   - Loads the time-series population/housing targets from `pop_projection_2050.csv`.
3. **Scenario Runs**  
   - For each scenario (e.g., LShH, LSpH, LSpF, LSpL) and for multiple iterations (`k` times), runs the relevant function (e.g., `LShH(env_table, vol_trgt)`).
   - Stores intermediate and final environment states (e.g., `env_table` in 2025, 2030, 2035...).
4. **Export Results**  
   - Calls `scen_sum(...)` to produce summary statistics.  
   - Calls `scen_export(...)` to save scenario outputs (rasters and CSV tables) to `./scenarios/.../`.

**Main sections in `scenario_run_.R`**:
- **`#### LShH - Sharing`**  
- **`#### LSpH - Sparing`**  
- **`#### LSpF - Full Sparing`**  
- **`#### LSpL - Local Sparing`**  

Each block runs the scenario over time steps, saves outputs, and logs metadata in CSV files.

---

### 3. `SDM_fnct_.R`
**Purpose**: Contains helper functions used for preparing inputs to Maxent models.

**Key Function**:
- **`get_focals(predictors, focal_list)`**  
  - For each raster in `focal_list`, computes a 5×5 moving window average (500 m if pixels are 100 m).  
  - Returns an expanded raster stack with new variables named, e.g., `tree_cover_500`.

This script is primarily sourced by the Maxent scripts to prepare or augment the predictor raster stack.

---

### 4. `SDM_maxent_.R`
**Purpose**: Runs species distribution models (SDMs) using **presence–absence** or **presence–pseudoabsence** data from a standardized **bird survey** dataset.

**Workflow**:
1. **Data Import**  
   - Loads `SDM_fnct_.R` for the `get_focals(...)` function.  
   - Reads table of bird survey points (e.g., `birdsurvey_alldata_clean.csv`) and shapefile of point coordinates.  
   - Reads base predictor raster stack from `./100_grid_variables/`.
2. **Predictor Preparation**  
   - Uses `get_focals(...)` to compute neighborhood variables (e.g., 500 m means).
3. **Cross-validation Setup**  
   - Defines number of folds (`k_fold = 5`) and repeated cross-validation runs (`k_cross`).
   - Creates presence-background or presence-absence sets.
   - Implements spatial blocking (via `blockCV` or custom code).
4. **Maxent Model Fitting**  
   - For each species, runs cross-validation (5-folds × `k_cross` repeats).  
   - Evaluates models via TSS and Boyce Index.  
   - Exports summary metrics, predictions, and model objects to `./Maxent/...`.
5. **Final Model**  
   - Fits a final Maxent model on all filtered presences plus background.  
   - Exports final predicted surfaces (`prediction.grd`), response curves, and evaluation metrics.
6. **Scenario Prediction**  
   - Optionally loops over the scenario rasters (produced by `scenario_run_.R`) to generate future distribution maps for each species under each scenario.

---

### 5. `SDM_maxent_ebird_.R` (Optional)
**Purpose**: Similar to `SDM_maxent_.R` but specifically handles **presence-only** eBird data.  
**Main differences**:
- Incorporates a **bias layer** (e.g., `bias_TAD`) to address uneven sampling effort typical in citizen-science data.  
- Thins occurrences spatially to reduce auto-correlation (via `GeoThinneR` or internal code).  
- Uses presence-background approach with weighting according to the bias layer in Maxent.

This script follows a similar cross-validation procedure (spatially blocked) and produces final Maxent models plus predictions for eBird-based species distributions.

---

## Usage Examples

Below is a rough outline of how to replicate major analysis steps. Actual workflow depends on your data paths and environment:

1. **Set up directories and data**  
   - Organize your baseline rasters in `./100_grid_variables/`, your input tables in `./tables/`, etc.
   - Ensure your CSV tables (`env100_raster_table_.csv`, `pop_projection_2050.csv`, etc.) match the format expected by the scripts.

2. **Run Urban Growth Simulations**  
   - Open `scenario_run_.R`.
   - Adjust `main_path` to your repository root or a path containing the subfolders.  
   - Confirm scenario names, iteration counts (`k`), and time steps (`years`).  
   - Run the script. It will produce new scenario folders in `./scenarios/...` with CSV tables and GeoTIFF rasters for each scenario-year-iteration.

3. **Run SDM with Survey Data**  
   - Open `SDM_maxent_.R`.  
   - Confirm the path to your bird survey dataset (`birdsurvey_alldata_clean.csv`) and the location of your predictor rasters.  
   - Adjust cross-validation settings as needed.  
   - Run the script to generate Maxent models in `./Maxent/...`.

4. **(Optional) Run SDM with eBird Data**  
   - Open `SDM_maxent_ebird_.R`.  
   - Provide paths to eBird observations, bias layer (`bias_TAD`), and predictor rasters.  
   - This script yields final eBird-based species distribution models in a similar manner.

5. **Project Models onto Scenarios**  
   - Inside `SDM_maxent_.R` (or `SDM_maxent_ebird_.R`), search for the section labeled “Scenario SDMs” or “produce predictions.”  
   - Update `sce_names` and `years` to match the scenario outputs you produced with `scenario_run_.R`.  
   - Ensure the file paths to scenario rasters are correct.  
   - Run the code. It creates future distribution maps for each species under each scenario-year-iteration (in `./Maxent/[scenario]/[year]/iter_[k]/`).

---

## References

- **Maxent**: Phillips, S. J., Anderson, R. P., & Schapire, R. E. (2006). [Maximum entropy modeling of species geographic distributions. *Ecological Modelling*, 190(3–4), 231–259.](https://www.sciencedirect.com/science/article/pii/S030438000500267X)  
- **dismo R package**: Hijmans, R. J. et al. (2020). *dismo*: Species Distribution Modeling. [CRAN Link](https://CRAN.R-project.org/package=dismo)  
- **blockCV**: Valavi, R., et al. (2019). *blockCV*: An r package for generating spatially or environmentally separated folds for k-fold cross‐validation of species distribution models. *Methods in Ecology and Evolution*, 10(2), 225–232.  
- **ecospat**: Broennimann, O., et al. (2023). *ecospat*: Spatial Ecology Miscellaneous Methods. [CRAN Link](https://CRAN.R-project.org/package=ecospat)  
- **Hirzel et al. (2006)**: [Evaluating the ability of habitat suitability models to predict species presences. *Ecological Modelling*, 199(2), 142–152.](https://www.sciencedirect.com/science/article/pii/S0304380006001148)  
- **Mann-Kendall**: Kendall, M. G. (1975). *Rank Correlation Methods*. Griffin.  
- **GeoThinneR**: For thinning presence-only data. [GitHub Link](https://github.com/darunabas/GeoThinneR)

---


Please cite the corresponding publication once it is available if you use any portion of this code in a derivative work.

---

**Note**: All code is provided “as is” with no warranties. Adjust parameters (file paths, variable names, iteration settings) according to your data and computing resources.