# Data Dictionary — Harmonized Output

These CSVs are produced by the data-processing pipeline and consumed by the analysis repo. They are also the analysis-ready dataset published on Dryad.

## harmonized_response_ratios.csv

Log response ratios comparing MPA vs reference sites for 5 focal species across 3 monitoring programs.

| Column | Type | Description |
|--------|------|-------------|
| CA_MPA_Name_Short | text | Standardized MPA name |
| year | int | Survey year |
| y | text | Full scientific name (5 focal species) |
| lnDiff | float | Log response ratio: ln(MPA/Reference) |
| mpa | float | Corrected proportion value at MPA site |
| reference | float | Corrected proportion value at reference site |
| Diff | float | Raw ratio (mpa / reference) |
| resp | text | Response type: "Bio" (biomass) or "Den" (density) |
| time | int | Years since MPA implementation (0 = before/implementation year) |
| type | text | MPA designation type (SMR, SMCA, Special Closure) |
| Location | text | Region code |
| Hectares | float | MPA area in hectares |
| source | text | Monitoring program: PISCO, KFM, or LTER |
| BA | text | Period: "Before" or "After" MPA establishment |

### Species (y column)
- Panulirus interruptus (California spiny lobster — predator)
- Semicossyphus pulcher (California sheephead — predator)
- Strongylocentrotus purpuratus (purple sea urchin — herbivore)
- Mesocentrotus franciscanus (red sea urchin — herbivore)
- Macrocystis pyrifera (giant kelp — producer)

## harmonized_raw_responses.csv

Raw density and biomass measurements inside and outside MPAs.

| Column | Type | Description |
|--------|------|-------------|
| CA_MPA_Name_Short | text | MPA name |
| year | int | Survey year |
| taxon_name | text | Full scientific name |
| source | text | Monitoring program: PISCO, KFM, or LTER |
| status | text | Site status: "Inside" or "Outside" MPA |
| value | float | Density (individuals/m2) or biomass (g/m2) measurement |
| resp | text | Response type: "Bio" or "Den" |
| BA | text | Period: "Before" or "After" |
| time | int | Years since MPA implementation |

## harmonized_landsat_rr.csv

Satellite-derived kelp canopy response ratios from Landsat imagery. This dataset bypasses the main monitoring program pipeline (scripts 04-06) and is consumed directly by the effect size calculation.

| Column | Type | Description |
|--------|------|-------------|
| CA_MPA_Name_Short | text | MPA name |
| year | int | Survey year |
| lnDiff | float | Log response ratio |
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

## Notes

- **lnRR interpretation:** Positive lnDiff = higher abundance inside MPA. lnDiff = 0 means no MPA effect. exp(lnDiff) gives the response ratio (e.g., lnDiff = 0.69 means MPA site has ~2x the reference site value).
- **Time = 0:** Includes the MPA implementation year and all years before it.
- **Sources:** PISCO = Partnership for Interdisciplinary Studies of Coastal Oceans; KFM = NPS Kelp Forest Monitoring / MBON; LTER = Santa Barbara Coastal LTER.
