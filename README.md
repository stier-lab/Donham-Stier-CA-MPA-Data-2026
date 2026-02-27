# CA MPA Kelp Forest — Data Processing

This repository processes raw monitoring data from three long-term kelp forest monitoring programs (PISCO, KFM/MBON, LTER) and satellite imagery (Landsat) into harmonized, analysis-ready CSVs.

**Companion analysis repo:** [Donham-Stier-CA-MPA-2026](https://github.com/stier-lab/Donham-Stier-CA-MPA-2026) — consumes the harmonized CSVs produced here to run meta-analysis, generate figures, and produce manuscript tables.

## Quick Start

```r
# 1. Set up raw data (see Data Setup below)
# 2. Run the pipeline
source(here::here("code", "R", "run_all.R"))
# 3. Find harmonized CSVs in output/harmonized/
```

## Data Setup

Raw monitoring data (~1.3 GB) is not included in this repository. Choose one option:

### Option 1: Download from Dryad

Data is archived on Dryad (DOI: TBD). Download and extract into `data/`.

### Option 2: Symlink from Google Drive (lab members)

```bash
GDRIVE="/Users/$(whoami)/Library/CloudStorage/GoogleDrive-astier@ucsb.edu/My Drive/Stier Lab/People/Emily Donham/Projects/Kelp MPA/data"
PROJECT="$HOME/Donham-Stier-CA-MPA-Data-2026/data"

ln -s "$GDRIVE/MBON" "$PROJECT/MBON"
ln -s "$GDRIVE/PISCO" "$PROJECT/PISCO"
ln -s "$GDRIVE/LTER" "$PROJECT/LTER"
ln -s "$GDRIVE/LANDSAT" "$PROJECT/LANDSAT"
ln -s "$GDRIVE/MPA" "$PROJECT/MPA"
ln -s "$GDRIVE/ALL_sizefreq_2024.csv" "$PROJECT/ALL_sizefreq_2024.csv"
```

## Output

The pipeline produces 4 harmonized CSVs in `output/harmonized/`:

| File | Rows | Description |
|------|------|-------------|
| `harmonized_response_ratios.csv` | ~3,400 | Log response ratios (lnRR) by MPA/year/species/response |
| `harmonized_raw_responses.csv` | ~6,200 | Raw density/biomass values inside and outside MPAs |
| `harmonized_landsat_rr.csv` | ~720 | Satellite-derived kelp canopy response ratios |
| `harmonized_site_metadata.csv` | 34 | MPA site locations, establishment dates, types |

See `output/harmonized/DATA_DICTIONARY.md` for column-level documentation.

## Pipeline

```
00_libraries.R           → Load packages
00c_analysis_constants.R → Named constants
01_utils.R               → Utility functions (biomass, bootstrap, RR, validation)
02_pBACIPS_function.R    → Core statistical function
03_data_import.R         → Import raw data
04_pisco_processing.R    → PISCO data processing
05_kfm_processing.R      → KFM/NPS data processing
06_lter_processing.R     → LTER data processing
06b_landsat_processing.R → Landsat satellite data
07_combine_and_export.R  → Combine sources, export harmonized CSVs
```

## Citation

Donham, E. & Stier, A. (2026). Restoration of trophic cascades within kelp forests following the establishment of a network of marine protected areas. *Conservation Letters*.
