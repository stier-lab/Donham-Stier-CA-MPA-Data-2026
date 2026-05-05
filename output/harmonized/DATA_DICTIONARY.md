# Data Dictionary: Harmonized Output

These CSVs are produced by the data-processing pipeline and consumed by the analysis repo. They are also the analysis-ready dataset published on Dryad.

**Last updated:** 2026-05-04 (post Emily code review, metadata audit, and full data-source documentation)

**See also:** `../../DATA_SOURCES.md` in the data-processing repo root, the single source of truth for every INPUT data file (15 files, with citations, DOIs, EDI package IDs, row counts, schema, and which script consumes each one). This dictionary documents the OUTPUT side. `DATA_SOURCES.md` documents the input side.

---

## harmonized_response_ratios.csv

Log response ratios comparing MPA vs reference sites for 5 focal species across 3 monitoring programs.

| Column | Type | Description |
|--------|------|-------------|
| CA_MPA_Name_Short | text | Standardized MPA name |
| year | int | Survey year |
| y | text | Full scientific name (5 focal species) |
| lnDiff | float | Log response ratio: ln(MPA/Reference), computed on time-series-max-normalized proportions |
| mpa | float | Zero-corrected proportion at MPA site (proportion of MPA-taxon time-series max) |
| reference | float | Zero-corrected proportion at reference site (same normalization) |
| Diff | float | Raw ratio (mpa / reference) |
| resp | text | Response type: "Bio" (biomass) or "Den" (density) |
| time | int | Years since MPA implementation (0 = before/implementation year) |
| type | text | MPA designation type (SMR, SMCA, Special Closure) |
| Location | text | Region code |
| Hectares | float | MPA area in hectares |
| source | text | Monitoring program: PISCO, KFM, or LTER |
| BA | text | Period: "Before" or "After" MPA establishment |

### Species (y column)
- Panulirus interruptus (California spiny lobster, predator)
- Semicossyphus pulcher (California sheephead, predator)
- Strongylocentrotus purpuratus (purple sea urchin, herbivore)
- Mesocentrotus franciscanus (red sea urchin, herbivore)
- Macrocystis pyrifera (giant kelp, producer)

## harmonized_raw_responses.csv

Raw density and biomass measurements inside and outside MPAs.

| Column | Type | Description |
|--------|------|-------------|
| CA_MPA_Name_Short | text | MPA name |
| year | int | Survey year |
| taxon_name | text | Full scientific name |
| source | text | Monitoring program: PISCO, KFM, or LTER |
| status | text | "mpa" or "reference" |
| value | float | Density (individuals/m²) or biomass (g/m²) measurement |
| resp | text | Response type: "Bio" or "Den" |
| BA | text | Period: "Before" or "After" |
| time | int | Years since MPA implementation |

## harmonized_landsat_rr.csv

Satellite-derived kelp canopy response ratios from Landsat imagery (Bell/SBC LTER kelp-watch product). This dataset bypasses the main monitoring-program pipeline and is consumed directly by the effect-size calculation.

| Column | Type | Description |
|--------|------|-------------|
| CA_MPA_Name_Short | text | MPA name |
| year | int | Survey year |
| lnDiff | float | Log response ratio (annual MAX kelp biomass MPA vs reference) |
| BA | text | Before/After MPA establishment |
| time | int | Years since MPA implementation |
| y | text | Always "Macrocystis pyrifera" |
| resp | text | Always "Bio" |
| source | text | Always "Landsat" |
| type | text | MPA designation type |
| Location | text | Region code |
| Hectares | float | MPA area in hectares |

## harmonized_site_metadata.csv

MPA site information including geographic coordinates and establishment dates.

| Column | Type | Description |
|--------|------|-------------|
| Lat | float | Latitude (decimal degrees) |
| Lon | float | Longitude (decimal degrees, negative = west) |
| Location | text | Region code |
| MPA_Start | int | Year MPA was legally established |
| Quarter | int | Quarter code |
| CA_MPA_Name_Short | text | Standardized MPA name (join key) |
| Site | text | Alternative site name |
| type | text | MPA designation: SMR, SMCA, or Special Closure |
| ChannelIsland | text | Island name, or NA for mainland sites |

---

## Methodology Notes (audit verified 2026-05-04)

### Sampling unit conventions per program

These were audited directly against the raw data files. Numbers in parentheses are observation counts.

| Program | Taxon | Transect / quadrat area | Notes |
|---|---|---|---|
| PISCO swath | Urchins, lobster, kelp | 60 m² (30 m × 2 m) | No AREA column in raw data; constant by protocol. |
| PISCO fish (UCSB / VRG) | Sheephead | 60 m² seafloor footprint | Counts recorded at vertical levels BOT (29,415) / MID (323) / CNMD (47) / CAN (7); CAN excluded by default; toggle `SHEEPHEAD_BOT_ONLY` available for BOT-only sensitivity. |
| KFM (MBON quad/swath) | Urchins | **1 m² (1982–1984)**, **2 m² (1985+)** | Quad-area changed through time. Density divides by row-level `area` not by the constant `KFM_QUAD_AREA_M2 = 2`. |
| KFM (MBON quad/swath) | Macrocystis | 1 m² (1982-84), 2 m² (1985+ paired-diver), or 10 m² (1996+ paired-diver Macro-specific 5 m x 1 m protocol) | **Pipeline filters to 10 m² subset only** (the modern 1996+ Macro-targeted protocol). MBON's `area = 10` is the paired-diver convention: each row is one diver's 5 m² half-quad, and the CountA + CountB pair = 10 m². See Davis et al. 1997 KFM Handbook Vol 1 pp.13, 38. |
| KFM (MBON quad/swath) | Lobster (P. interruptus) | 40 m² (1983-1984), 60 m² (1985+) | Davis et al. 1997 p.39: 1983-84 used 10 transects/site at 2 m × 20 m. 1985-present uses 12 transects/site at 3 m × 20 m. Density divides by row-level `area`, so the temporal change is handled correctly. KFM does not record lobster sizes, so KFM contributes density only, no biomass. |
| KFM (MBON fish file) | Sheephead | Visual Fish Transect (since 1985, transect length 100 m -> 50 m in 1996) + Roving Diver Fish Count (RDFC, since 1996, not 2003) | KFM does not record fish sizes, so KFM contributes density only, no biomass. Per Davis et al. 1997 p.45-48. |
| LTER fish | Sheephead | **80 m² (40 m × 2 m, "BIG fish" protocol)** | LTER survey actually uses two transect sizes: 80 m² for big mobile fish, 20 m² for cryptic small fish. Sheephead is uniformly on 80 m². Cryptic fish are not in this analysis. Defensive check at runtime warns if SPUL transect area ever changes. |
| LTER lobster | Lobster | 300 m² (60 m × 5 m, 2012+) | Constant by protocol (verified). |
| LTER kelp | Macrocystis | 20 or 40 m² | Density normalization uses row-level `AREA`. LTER stores stipe counts in the `FRONDS` column. Per Li (LTER), fronds = stipes for this purpose. |
| Landsat (Bell/SBC LTER) | Macrocystis (canopy biomass only) | Per-MPA polygon | Annual values are **MAX** (peak) kelp canopy biomass within the year, in kg wet weight. Wide format: one row per MPA × status × replicate, one column per year (1984–2021). |

### Size-frequency pool composition (`data/ALL_sizefreq_2024.csv`)

The size-frequency pool that drives biomass bootstrapping was audited against the file's `method` column on 2026-05-04. **The urchin and lobster size pools are PISCO-only.** Method-code coverage:

- `SBTL_SIZEFREQ_PISCO`: PISCO UCSB measurements
- `SBTL_SIZEFREQ_VRG`: PISCO VRG measurements
- `SBTL_SWATH`: direct-on-swath measurements (lobster, UCSB protocol)
- `SBTL_SIZEFREQ_KFM`: present as a code but **0 observations** for urchins or lobster in the current 2024 release
- `SBTL_SIZEFREQ_LTER`: present as a code but **0 observations** for urchins

Urchin pool by method (MESFRA + STRPUR):

| Method | Observations |
|---|---|
| SBTL_SIZEFREQ_PISCO (UCSB) | 3,848 |
| SBTL_SIZEFREQ_VRG | 16,810 |
| SBTL_SIZEFREQ_KFM | 0 |
| SBTL_SIZEFREQ_LTER | 0 |

Lobster pool by method (PANINT):

| Method | Observations |
|---|---|
| SBTL_SIZEFREQ_PISCO (UCSB) | 122 |
| SBTL_SIZEFREQ_VRG | 1,940 |
| SBTL_SWATH (direct on transect, UCSB) | 2,109 |

**Implication:** KFM and LTER urchin biomass estimates rely on PISCO size data via cross-program bootstrap. KFM does collect its own urchin size-frequency data (~200 individuals per species per site per year since 1985, per Davis et al. 1997 p.28). Those records are not in `ALL_sizefreq_2024.csv` but are publicly available for the period 1985-2007 via NOAA ERDDAP at `https://coastwatch.pfeg.noaa.gov/erddap/tabledap/erdCinpKfmSFNH` (dataset "Channel Islands, Kelp Forest Monitoring, Size and Frequency, Natural Habitat, 1985-2007"). Post-2007 KFM size-frequency data require direct request from the NPS Channel Islands KFM program. Adding these data to `ALL_sizefreq_2024.csv` under method = `SBTL_SIZEFREQ_KFM` would strengthen the urchin biomass estimates from cross-program bootstrap to within-program estimation. This is a known future-work item, deferred for the current manuscript.

### Allometric biomass conversions

| Function | Used for | Source |
|---|---|---|
| `bio_lobster()` | Lobster CL (mm) → biomass (g) | SBC LTER allometric calibration (cf. Barsky 2001, CDFG) |
| `bio_redurch()` | Red urchin diameter (mm) → biomass (g) | SBC LTER allometric calibration (cf. Ebert 2010) |
| `bio_purpurch()` | Purple urchin diameter (mm) → biomass (g) | SBC LTER allometric calibration (cf. Ebert 2010) |
| `bio_sheephead()` | Sheephead TL (cm) → biomass (g) | SBC LTER allometric calibration (cf. Cowen 1990) |
| `bio_macro()` | Macrocystis stipe count → biomass (g) | Mean of six site-specific slopes from SBC LTER kelp allometry, **fit to May–October survey data only** (the LTER canopy-biomass calibration period). |

Sheephead biomass is only available from PISCO (size data) and LTER (size data). KFM lacks fish sizes and contributes density only.

### MPA exclusion lists

Maintained in `01_utils.R`:

- **`EXCLUDED_MPAS`**: MPAs excluded entirely (Carrington Point SMR, Anacapa Island SMR 1978 [pre-MLPA protection], etc.).
- **`EXCLUDED_LTER_SITES` / `EXCLUDED_LTER_SITE_CODES`**: Arroyo Hondo (AHND / a-l-02). Excluded because the site does not meet kelp-forest habitat criteria. It is dominated by surfgrass / sand with sparse, ephemeral kelp coverage and lacks a persistent giant-kelp forest community. Updated 2026-04 from earlier "data quality" framing.
- **`EXCLUDED_REFERENCE_SITES`**: 16 PISCO reference sites flagged by Emily during the 2026-04 review (e.g., ANACAPA_BLACK_SEA_BASS, SMI_PRINCE_ISLAND_*).
- **`EXCLUDED_KFM_SITES`**: 14 KFM transects within the Anacapa Island 1978 SMR boundary, excluded to prevent pre-MLPA protection from confounding the MLPA-era before/after analysis.
- **`SHEEPHEAD_ONLY_MPAS`**: 11 MPAs with sheephead-only data. Removed from cross-taxa pooled analyses to prevent over-representation, selectively re-included for sheephead-targeted analyses via `reinclude_sheephead_mpas()`.
- **`SHEEPHEAD_REINCLUSION_MPAS`** (updated 2026-04): MPAs whose regulations actually protect sheephead (no take of finfish), so the MPA functions like an SMR from sheephead's perspective. Final list is **Blue Cavern Onshore SMCA, Farnsworth Onshore SMCA, Swamis SMCA, Long Point SMR**. Dana Point SMCA and Cat Harbor SMCA were removed (Dana Point permits take of lobster, finfish, and urchins; Cat Harbor permits commercial take of lobster and urchin and recreational take of finfish).

### Proportion-based lnRR (non-standard)

This pipeline computes lnRR on **proportions of time-series maxima**, not on raw means. The proportion is the value at year *t* divided by the per-MPA × taxon × status maximum across the full time series, with a small zero-correction (`adaptive` mode = half the minimum non-zero proportion within the group, `fixed` mode = +0.01 for sensitivity). This non-standard approach standardizes across programs that use different transect areas, survey methods, and spatial extents, at the cost of giving up an absolute-scale lnRR. The choice is invariant to multiplicative scaling within a group, which is why the KFM urchin quad-area fix has zero impact on response ratios (response ratios are unchanged even though absolute densities are now correct for 1982–1984).

### Sheephead level filter (PISCO fish, audit 2026-05-04)

Distribution of all fish-level codes in the PISCO fish file: BOT (359,695) / MID (104,438) / CNMD (10,856) / CAN (42,120). For sheephead specifically: BOT (29,415) / MID (323) / CNMD (47) / CAN (7). The pipeline default is `SHEEPHEAD_BOT_ONLY = FALSE` (BOT/MID/CNMD summed; CAN excluded), which captures sheephead in the water column. Switch to `TRUE` for BOT-only (Hamilton replication / sensitivity).

### Notes

- **lnRR interpretation:** Positive lnDiff = higher abundance inside MPA. lnDiff = 0 means no MPA effect. `exp(lnDiff)` gives the response ratio (e.g., lnDiff = 0.69 means MPA site has ~2× the reference site value).
- **Time = 0:** Implementation year. Years before are negative, years after are positive.
- **Sources:** PISCO = Partnership for Interdisciplinary Studies of Coastal Oceans. KFM = NPS Kelp Forest Monitoring (delivered through MBON synthesis). LTER = Santa Barbara Coastal LTER. Landsat = Bell/SBC LTER kelp-watch satellite product.
