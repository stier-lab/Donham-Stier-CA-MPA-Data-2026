# CA MPA Kelp Forest: Data Processing

This repository processes raw monitoring data from three long-term kelp forest monitoring programs (PISCO, KFM/MBON, LTER) and satellite imagery (Landsat) into harmonized, analysis-ready CSVs.

**Companion analysis repo:** [sbc-2026-donham-kelp-mpa-cascade](https://github.com/stier-lab/sbc-2026-donham-kelp-mpa-cascade). It consumes the harmonized CSVs produced here to run meta-analysis, generate figures, and produce manuscript tables.

**GitHub status:** archived/read-only after the 2026-09-03 cleanup. This repo
is provenance and harmonization support for the current analysis repo, not the
final manuscript data/code release by itself.

## Repository Scope

This public repository is the source-data harmonization workspace. Keep raw data
out of git, keep manuscript drafts/citation PDFs/cover letters in the private
manuscript repo, and keep statistical models, figures, and result summaries in
the companion analysis repo.

## Quick Start

```r
# 1. Set up raw data (see Data Setup below)
# 2. Run the pipeline
source(here::here("code", "R", "run_all.R"))
# 3. Find harmonized CSVs in output/harmonized/
```

## Data Setup

Raw monitoring data are not tracked in this repository. To rerun harmonization,
obtain the source files from the providers or project data stewards described
in `DATA_SOURCES.md`, then arrange local, untracked inputs under:

```text
data/MBON/
data/PISCO/
data/LTER/
data/LANDSAT/
data/MPA/
data/ALL_sizefreq_2024.csv
```

The current tracked harmonized outputs are included under `output/harmonized/`.
Do not add raw monitoring extracts, lab Drive paths, or local symlinks to git.

## Output

The pipeline produces 4 harmonized CSVs in `output/harmonized/`:

| File | Rows | Description |
|------|------|-------------|
| `harmonized_response_ratios.csv` | ~3,400 | Log response ratios (lnRR) by MPA/year/species/response |
| `harmonized_raw_responses.csv` | ~6,200 | Raw density/biomass values inside and outside MPAs |
| `harmonized_landsat_rr.csv` | ~720 | Satellite-derived kelp canopy response ratios |
| `harmonized_site_metadata.csv` | 34 | MPA site locations, establishment dates, types |

## Documentation

This repo has two documentation files at different layers of the pipeline:

- **`DATA_SOURCES.md`** (in repo root): single source of truth for every INPUT data file. Lists all 15 input files with their canonical citations, DOIs, EDI package IDs, row counts, schema, and the script that consumes each one. Use this when you need to know where the raw data came from or how to cite a specific source.

- **`output/harmonized/DATA_DICTIONARY.md`**: column-level documentation for the four OUTPUT CSVs above, plus an integrated "Methodology Notes" section covering per-program transect-area conventions, the size-frequency pool composition, allometric source attributions, MPA exclusion lists, and the proportion-based lnRR rationale. Use this when reading the harmonized data or designing analyses that consume them.

## Pipeline

```
00_libraries.R           → Load packages
00c_analysis_constants.R → Named constants
01_utils.R               → Utility functions (biomass, bootstrap, RR, validation)
02_pBACIPS_function.R    → Core statistical function
03_data_import.R         → Import raw data (size-frequency pool, site list, MPA features, MBON site geolocations); see DATA_SOURCES.md §5, §6
04_pisco_processing.R    → PISCO data processing; see DATA_SOURCES.md §3
05_kfm_processing.R      → KFM/NPS data processing; see DATA_SOURCES.md §1.1, §1.2, §1.3, §4
06_lter_processing.R     → LTER data processing; see DATA_SOURCES.md §2
06b_landsat_processing.R → Landsat satellite data; see DATA_SOURCES.md §7
07_combine_and_export.R  → Combine sources, export harmonized CSVs
```

## Citation

Donham, E. & Stier, A. (2026). Partial trophic-cascade recovery and kelp resilience in a southern California marine protected-area network. Manuscript in preparation for *Journal of Applied Ecology*.
