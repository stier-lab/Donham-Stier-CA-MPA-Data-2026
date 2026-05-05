# Methods References: CA MPA Kelp Forest pBACIPS Manuscript

**Compiled:** 2026-05-04
**Purpose:** Authoritative citations and verbatim protocol excerpts for every methodological claim in the data pipeline. This is the source-of-truth for the manuscript Methods section. Findings synthesized from a five-track audit of EDI/SBC LTER metadata, NPS Channel Islands KFM protocol, PISCO/MLPA monitoring documentation, the Bell/Cavanaugh kelp-watch product, and a full data-source inventory of every input file in the data-processing repo.

**Companion data-source inventory:** `Donham-Stier-CA-MPA-Data-2026/DATA_SOURCES.md` is the parallel single-source-of-truth document inside the data-processing repo, documenting all 15 input data files (citations, DOIs, EDI package IDs, row counts, schema, and which script consumes each). This Methods doc focuses on the manuscript-citation perspective. The in-repo `DATA_SOURCES.md` focuses on the pipeline perspective.

---

## Audit summary (closed 2026-05-04)

Items below are listed in order of severity. All have been resolved or deferred with documentation. Code changes are committed locally in both repos.

| Item | What we had said | What the audit found | Status |
|---|---|---|---|
| **Landsat aggregation** | "MAX annual canopy biomass" | Confirmed via three independent in-repo sources (Donham archive code L5854, current 06b header, DATA_DICTIONARY): `MPA_Runs_new.csv` is annual MAX-of-quarters, produced by Tom Bell from knb-lter-sbc.74. Consistent with Bell et al. 2023 Kelpwatch convention. | RESOLVED |
| **KFM Macrocystis quad area** | "10 m²" | MBON's `area = 10` is the paired-diver convention for the modern (1996+) 5 m x 1 m KFM Macro protocol: each row = one diver's 5 m² half-quad, and the CountA + CountB pair covers 10 m². Our filter is correct (it selects the 1996+ Macro-targeted protocol). Documented in 05_kfm_processing.R. | RESOLVED |
| **KFM RDFC start year** | 2003 | Davis et al. 1997 p.48 confirms RDFC began in 1996, not 2003. Visual Fish Transect since 1985 (with 100m→50m length change in 1996). Constant updated to 1996 in both repos. | CORRECTED |
| **KFM urchin size-frequency** | "Pool is PISCO-only" (Sept earlier note) | KFM does collect ~200 urchin sizes/species/site/year since 1985 (Davis 1997 p.28). Public ERDDAP feed `erdCinpKfmSFNH` covers 1985-2007 only. Post-2007 data require direct request from NPS Channel Islands. Ingestion deferred for this manuscript. Cross-program PISCO-derived bootstrap is documented as a known assumption. | DEFERRED (with documented path) |
| **KFM lobster transect history** | "40 or 60 m² (mostly 60)" | Davis et al. 1997 p.39: 40 m² (1983-1984 only) and 60 m² from 1985+. The split is temporal, not site-based. Comment in 05_kfm_processing.R updated to specify the transition year and describe the dive mechanism. | DOCUMENTED |
| **KFM sheephead = density only** | Correct | Davis et al. 1997 p.47 confirms: counts by sex/life-stage, no individual lengths recorded. | NO CHANGE NEEDED |
| **PISCO swath = 60 m²** | Correct | Multiple PISCO papers confirm 30 m × 2 m = 60 m². Canonical data-paper citation is Malone et al. 2022 *Ecology*. | NO CHANGE NEEDED |
| **PISCO sheephead level codes** | BOT/MID/CNMD/CAN | Confirmed. CNMD = "canopy-midwater" (used at sites where canopy is too shallow for a separate canopy diver). | NO CHANGE NEEDED |
| **`MACRO_AVE_SLOPE` calibration window** | "May–Oct only" | NOT supported by public LTER protocol, which uses month-specific (12 monthly) slopes from Rassweiler 2018. The "May-Oct" attribution traces to Emily's archived comment (`pBACIPS_PISCO_V10.R` L674: "Slopes from May-Oct from LTER"). The values themselves can't be matched to any open-access table. Code docstring updated to remove the unverified assertion. Emily action item: cross-check against Rassweiler 2018 supplement via UCSB Wiley access. | CORRECTED (with one Emily action item) |
| **VRG roving-diver lobster sizing** | No formal source | Carr et al. 2021 OPC Kelp Forest Technical Report (p.20-21) has the verbatim "benthic swaths or haphazard surveys across depth zones" language. Backed by Malone et al. 2022 *Ecology*. | RESOLVED |
| **FRONDS = stipes (Li's verbal note)** | Correct | SBC LTER kelp protocol defines a frond as "any blade > 10 cm long that is still attached to the stipe" measured "at 1 m above the holdfast." Each frond = one stipe at the 1-m cutoff. PISCO, KFM, and LTER all sample at this height. | NO CHANGE NEEDED |

---

## 1. SBC LTER (Santa Barbara Coastal Long Term Ecological Research)

### 1.1 Annual Kelp: Macrocystis frond/stipe density

- **EDI package:** `knb-lter-sbc.18.29`
- **DOI:** `10.6073/pasta/9bdb2dcaece1f4ab933b54f7b8b68144`
- **URL:** https://portal.edirepository.org/nis/metadataviewer?packageid=knb-lter-sbc.18.29
- **Citation:** Reed, D.C. & R.J. Miller. 2024. *SBC LTER: Reef: Kelp Forest Community Dynamics: Abundance and size of Giant Kelp (Macrocystis pyrifera), ongoing since 2000.* Environmental Data Initiative.

**FRONDS column definition (verbatim):** "total number of fronds per plant at 1 meter from top of holdfast." Companion field protocol PDF (https://sbclter.msi.ucsb.edu/external/Reef/Protocols/Kelp_Forest_Community_Dynamics/SBC_LTER_protocol_Reed_Kelp_forest_community_Density_giant_kelp_20130524.pdf) defines a frond as **"any blade >10 cm long that is still attached to the stipe."** Operationally a stipe count of mature/contributing stipes at 1 m above holdfast.

**AREA definition:** "the density of giant kelp greater than 1 m tall is recorded in four contiguous 20 m × 1 m swaths that run parallel and adjacent to permanent 40 m transects." Per-plant rows are AREA = 20 m² (single section × 1 m wide) or AREA = 40 m² (combined section × 1 m wide).

### 1.2 Annual Reef Fish

- **EDI package:** `knb-lter-sbc.17.29`
- **URL:** https://portal.edirepository.org/nis/metadataviewer?packageid=knb-lter-sbc.17.29
- **Citation:** Reed, D.C. & R.J. Miller. *SBC LTER: Reef: Kelp Forest Community Dynamics: Abundance and size of reef fish, ongoing since 2000.* Environmental Data Initiative.

**Methods (verbatim):** "The number, size and species identity of reef fish are recorded within a 2 m wide swath centered along each transect extending 2 m off the bottom."

**Two transect protocols (confirmed):**
- **BIG / mobile fish** (sheephead included): full 40 m × 2 m = **80 m²** transect, 2 m water-column slab off the bottom.
- **Cryptic small benthic fish:** single 20 m sub-section, 1 m wide each side = **20 m²** (effectively a quarter of the BIG transect). Targets gobies, blennies, clinids, scorpaenids, juvenile kelpfish.

Sheephead (SP_CODE = SPUL) is uniformly on the 80 m² protocol, which confirms our `LTER_FISH_SURVEY_AREA_M2 = 80` constant.

### 1.3 Lobster Abundance

- **EDI package:** `knb-lter-sbc.77.10`
- **DOI:** `10.6073/pasta/63eca8e267cdaf016709e5c1b3d746d4`
- **URL:** https://portal.edirepository.org/nis/metadataviewer?packageid=knb-lter-sbc.77.10
- **Citation:** Reed, D.C. & R.J. Miller. 2024. *SBC LTER: Reef: Abundance, size and fishing effort for California Spiny Lobster (Panulirus interruptus), ongoing since 2012.* Environmental Data Initiative.

**Methods (verbatim):** "The number and sizes of spiny lobster are recorded by divers in four 60 × 5 m (300 m² total) lobster transects that are associated with each permanent SBC LTER kelp forest community survey transect." Surveys annual, late summer, before fishing season. **60 × 5 m transects added in 2012** because lobster densities in 40 × 2 m kelp-community transects were too low for reliable abundance estimates.

Five sites: two MPA (Naples, Isla Vista; designated January 2012) and three reference (Arroyo Quemado, Mohawk, Carpinteria). Confirms `LTER_LOBSTER_PLOT_AREA_M2 = 300`.

### 1.4 Macrocystis stipe→biomass allometry

- **Citation:** Rassweiler, A., Reed, D.C., Harrer, S.L., & Nelson, J.C. 2018. Improved estimates of net primary production, growth, and standing crop of *Macrocystis pyrifera* in Southern California. *Ecology* 99(9): 2132. https://doi.org/10.1002/ecy.2440 (PMID: 29956835)
- **Companion EDI package:** `knb-lter-sbc.112` (frond-morphology / NPP), https://portal.edirepository.org/nis/mapbrowse?scope=knb-lter-sbc&identifier=112
- **NPP protocol PDF:** https://sbclter.msi.ucsb.edu/external/Reef/Protocols/Kelp_NPP/KelpNPP_20180522.pdf

**Caveat:** The "May–Oct calibration window" claim could not be verified from the open-access abstract alone. Need to pull the *Ecology* supplement (https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecy.2440 via UCSB) to confirm or revise. The slopes hardcoded in our `MACRO_AVE_SLOPE` (0.10386, 0.10103, 0.09267, 0.09204, 0.08054, 0.08505) likely come from the supplement's per-site, per-season table.

---

## 2. KFM (Kelp Forest Monitoring, Channel Islands National Park)

### 2.1 Authoritative protocol document

**Citation:** Davis, G.E., D.J. Kushner, J.M. Mondragon, J.E. Mondragon, D. Lerma & D.V. Richards. 1997. *Kelp Forest Monitoring Handbook, Volume 1: Sampling Protocol.* Channel Islands National Park, Ventura, CA. November 1997. 53 pp.

**URL:** https://irma.nps.gov/DataStore/DownloadFile/485444 (NPS IRMA DataStore)

This is the foundational protocol document, adopted in 1997 following a formal program review (Davis et al. 1996), and is still the operational protocol.

### 2.2 Quad-size history (confirmed verbatim, p.35)

Quoted directly from Davis et al. 1997 p.35 (verified by reading the PDF, section heading "1 m Quadrats" then "Sample size and database anomalies"):

> "The number of quadrats sampled per year has changed several times with 30 sampled in 1982, 40 sampled in 1983 and 1984, 20 sampled from 1985–1995, and 12 sampled from 1996 to present. Each quadrat for 1982, 1983, and 1984 represents 1m². From 1985 to present each quadrat represents 2m² (Tablexxx)."

The companion table on p.35 lists the canonical schedule:

| Year | # quadrats | Quadrat size | CountA | CountB |
|---|---|---|---|---|
| 1982 | 30 | 1 m² | yes | null |
| 1983–1984 | 40 | 1 m² | yes | null |
| 1985–1994 | 20 | 2 m² | yes | null *(two divers' counts summed before entry)* |
| 1995 | 20 | 2 m² | yes | yes *(per-diver entries)* |
| 1996–present | 12 | 2 m² | yes | yes |

The 1985 change was a **paired-diver redesign:** two divers each sample adjacent 1 m² quadrats on opposite sides of the transect line at the same meter mark, totaling 2 m² per replicate. From 1985–1994 the two counts were summed before entry. From 1995 they were entered separately as CountA / CountB. Our pipeline correctly handles this via row-level `area`.

### 2.3 KFM Macrocystis protocol (CORRECTION NEEDED)

**Original 5 m × 1 m = 5 m² Macro quads added 1996.** Our pipeline filters to `area == 10` in the MBON-integrated data, which is likely an MBON aggregation of the paired-diver pair (5 m² × 2 divers = 10 m²) or could be a separate 5 m × 2 m quad. Need to verify against MBON metadata.

KFM divers also count stipes on **exactly 100 plants per site per year** along the full 100 m main transect with a variable swath (1–5 m, density-dependent), measuring (a) stipe count 1 m above the substrate and (b) maximum holdfast diameter (cm). This is a **separate dataset** from the quadrat counts. It is filed under the Size Frequency protocol, table `SizeFreqNatHab`, and matches Emily's `KFM_Macrocystis_RawData` recollection.

### 2.4 KFM urchin size-frequency surveys (FILLS THE GAP EMILY FLAGGED)

**KFM does conduct site-level urchin size-frequency sampling.** Verbatim from Davis et al. 1997 p.28:
> "*Strongylocentrotus purpuratus* — 200 — Max. test diameter, mm
> *Strongylocentrotus franciscanus* — 200 — Max. test diameter, mm"

Methods (p.26): band-transect-style searches along the 100 m main transect, target sample size **200 individuals per urchin species per site per year**, measured to nearest mm with stainless-steel calipers. Implemented in 1985 (p.50).

**Available at:** NOAA ERDDAP `erdCinpKfmSFNH` (https://coastwatch.pfeg.noaa.gov/erddap/tabledap/erdCinpKfmSFNH.graph), "Channel Islands, Kelp Forest Monitoring, Size and Frequency, Natural Habitat, 1985–2007", with updated tables on NPS IRMA. **The MBON synthesis didn't pull these in. This is a data-integration gap, not a data-existence gap.**

This is the single biggest pipeline-improvement opportunity from the audit: ingesting these data would change the urchin biomass bootstrap from "cross-program PISCO-derived assumption" to "within-program KFM-derived measurement." Worth doing before submission.

### 2.5 KFM lobster transect history (CORRECTION NEEDED, temporal split)

> "In 1983 and 1984 10 band transects at each site were taken. Each of these transects represented 2×20 m or 40 m². From 1985 to the present, 12 band transects were conducted at each site. These transects measure 3×20 m or 60 m²." (p.39)

Mechanism: two divers each lay a 10 m perpendicular tape from the main transect, then count organisms within **1.5 m on each side** of their tape. So each diver = 3 m × 10 m = 30 m², and CountA + CountB = 3 m × 20 m = 60 m². Lobster sampling began 1983 (p.40).

### 2.6 KFM fish surveys

- **Visual Fish Transects (1985–present).** Length change: 1985–1995 was 3 m × 2 m × 100 m (= 600 m³ slab), 1996+ is 3 m × 2 m × 50 m (= 300 m³). Four transects per site throughout. Transects 1+2 (1996+) sum to old transect 1.
- **Roving Diver Fish Count (RDFC), 1996 onward** (NOT 2003 as in our constant). 30-minute roving counts, abundance scored on a 10-point time-of-encounter scale.
- **Sheephead recorded as counts by sex/life-stage:** Male / Female / Juvenile (juvenile diagnostic: "<10 cm length, white stripe"). **No per-fish length recorded**, which confirms our pipeline's correct treatment of KFM sheephead as density-only.

---

## 3. PISCO MLPA monitoring program

### 3.1 Recommended Methods-section citation stack

- **Malone, D.P., Davis, K., Lonhart, S.I., Parsons-Field, A., Caselle, J.E., & Carr, M.H. (2022).** Large-scale, multidecade monitoring data from kelp forest ecosystems in California and Oregon (USA). *Ecology* 103(5): e3630. https://doi.org/10.1002/ecy.3630
 *Canonical PISCO data paper. Cite for swath dimensions and overall protocol.*

- **Caselle, J.E., Rassweiler, A., Hamilton, S.L., & Warner, R.R. (2015).** Recovery trajectories of kelp forest animals are rapid yet spatially variable across a network of temperate marine protected areas. *Scientific Reports* 5: 14102. https://doi.org/10.1038/srep14102
 *Channel Islands cohort methods. Verbatim:* "8 to 12 fish transects that measured 30 × 2 × 2 m at multiple levels in the water column: benthic, midwater, and kelp canopy (when present)."

- **Eisaguirre, J.H. et al. (2020).** Trophic redundancy and predator size class structure drive differences in kelp forest ecosystem dynamics. *Ecology* 101(5): e02993. https://doi.org/10.1002/ecy.2993
 *Most explicit transect description.* Verbatim: "12 'benthic' transects (30 × 2 m) are surveyed at each site between June and August to quantify densities of invertebrates and macroalgae… Sea stars and purple urchins greater than 2.5 cm in diameter and lobsters of all sizes are also counted."

- **Hamilton, S.L. & Caselle, J.E. (2015).** Exploitation and recovery of a sea urchin predator has implications for the resilience of southern California kelp forests. *Proc. R. Soc. B* 282: 20141817. https://doi.org/10.1098/rspb.2014.1817
 *Cite for the urchin 25 mm threshold and lobster predation context.*

- **Rassweiler, A. et al. (2020).** Roving divers surveying fish in fixed areas capture similar patterns in biogeography but different estimates of density when compared with belt transects. *Frontiers in Marine Science* 7: 272. https://doi.org/10.3389/fmars.2020.00272
 *Most explicit description of the 3-level fish design (BOT/MID/CAN).* Verbatim: "A diver at the bottom lays the transect line while counting fish within 2 m of the substrate, while another diver swims above, at roughly half the depth of the bottom diver, counting fish in the midwater. When kelp canopy is present, a third level is surveyed in the kelp canopy."

- **Spiecker, B. et al. (2025).** Measuring biological effectiveness across a very large, coherent network of coastal marine protected areas. *Ecological Applications* 35(5): e70074. https://pmc.ncbi.nlm.nih.gov/articles/PMC12231823/
 *Most current MLPA-network synthesis paper.*

### 3.2 Vertical-level codes (BOT/MID/CNMD/CAN)

CNMD = "canopy-midwater", used at sites where the kelp canopy is too shallow to support a separate canopy diver, so the upper diver covers a combined canopy-midwater slab. Documented operationally in **Malone et al. 2022 Ecological Archives metadata** rather than in a dedicated paper. CNMD spelling is correct (not "CNMID").

### 3.3 UCSB vs VRG lobster sizing (no single peer-reviewed source)

UCSB measures CL on the swath (Eisaguirre 2020). VRG uses roving-diver sizing where divers measure every lobster encountered, best documented in:
- **PISCO Subtidal Community Surveys: Size Frequency Surveys** dataset, CNRA Open Data, doi:10.6085/AA/pisco_subtidal
- **Pondella, D.J. II et al. (2015).** Bulletin of the Southern California Academy of Sciences (Occidental kelp-forest program description)
- **VRG Dive Program page:** https://www.oxy.edu/academics/vantuna-research-group/dive-program

### 3.4 Urchin 25 mm threshold

PISCO size-frequency protocol-specific: urchins of *Strongylocentrotus*/*Mesocentrotus* sized at 5 mm bins from 25 mm test diameter upward. Smaller individuals are too cryptic for visual census. Cite **PISCO Subtidal Size Frequency Surveys** dataset + **Hamilton & Caselle 2015** *Proc. R. Soc. B*.

### 3.5 Sizing start years per campus

Documented only in dataset metadata. UCSB lobster sizing began ~2010 with Channel Islands MPA decadal cohort, VRG ~2008. UCSB urchin sizing: 2003–2006, gap, resumed 2019. VRG: ~2007 onward. Cite the PISCO dataset DOI + the **MLPA Decadal Review Report** (below) for these temporal patterns.

### 3.6 MLPA Decadal Management Review

- **CDFW (2023).** *2022 Master Plan for Marine Protected Areas / Decadal Management Review of California's Marine Protected Area Network.* California Department of Fish and Wildlife, with the California Fish and Game Commission and Ocean Protection Council. January 2023.
 https://nrm.dfg.ca.gov/FileHandler.ashx?DocumentID=209209&inline | landing page: https://wildlife.ca.gov/Conservation/Marine/MPAs/Management/Decadal-Review

- **Carr, M.H. et al. (2022).** *Partnership for the Interdisciplinary Studies of Coastal Oceans (PISCO) Contributions, Challenges, and Recommendations for the MLPA Decadal Evaluation Review.* Report submitted to CDFW Marine Region, February 4, 2022. https://nrm.dfg.ca.gov/FileHandler.ashx?DocumentID=203484

---

## 4. Landsat-derived Macrocystis canopy biomass (kelp-watch product)

### 4.1 Canonical citations

- **Cavanaugh, K.C., Siegel, D.A., Reed, D.C., & Dennison, P.E. (2011).** Environmental controls of giant-kelp biomass in the Santa Barbara Channel, California. *Marine Ecology Progress Series* 429: 1–17. https://doi.org/10.3354/meps09141
 *Primary methods citation. Biomass derivation, fractional-cover→biomass calibration.*

- **Bell, T.W., Allen, J.G., Cavanaugh, K.C., & Siegel, D.A. (2020).** Three decades of variability in California's giant kelp forests from the Landsat satellites. *Remote Sensing of Environment* 238: 110811. https://doi.org/10.1016/j.rse.2018.06.039
 *Multi-sensor harmonization (Landsat 5/7/8). Cite for the 1984–2021 time series.*

- **Bell, T.W., Cavanaugh, K.C., Saccomanno, V.R., Houskeeper, H.F., Eddy, N. et al. (2023).** Kelpwatch: A new visualization and analysis tool to explore kelp canopy dynamics reveals variable response to and recovery from marine heatwaves. *PLOS ONE* 18(3): e0271477. https://doi.org/10.1371/journal.pone.0271477
 *Public portal / state-wide aggregation. Verbatim:* "the maximum canopy area was determined for each cell for each year from 1984 to 2021". Kelpwatch's annual aggregation IS max-of-quarters.

### 4.2 Data package

- **EDI package:** `knb-lter-sbc.74.18+`
- **Title:** *SBC LTER: Time series of quarterly NetCDF files of kelp biomass in the canopy from Landsat 5, 7, 8 and 9, since 1984 (ongoing).*
- **URLs:** https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-sbc.74.18 ; https://sbclter.msi.ucsb.edu/data/catalog/package/?package=knb-lter-sbc.74

### 4.3 Methodological details (verbatim)

- **Computed quantity:** giant-kelp canopy **biomass** (kg wet weight) derived from satellite surface-reflectance-based fractional cover, calibrated against SCUBA biomass estimates: *"Estimates of kelp canopy biomass are derived from the relationship between giant kelp fractional cover determined from satellite surface reflectance and empirical measurements of giant kelp canopy biomass in long-term SBC LTER study plots obtained using SCUBA."*
- **Spatial unit:** 30 × 30 m Landsat pixel. Geographic coverage Año Nuevo CA south through Baja California for biomass.
- **Temporal unit (as published):** **quarterly mean** (Q1 Jan–Mar, Q2 Apr–Jun, Q3 Jul–Sep, Q4 Oct–Dec).
- **Annual aggregation (Kelpwatch convention):** max-of-quarters per pixel per year.

### 4.4 Pipeline action item

`MPA_Runs_new.csv` has one column per year (1984–2021). We could not determine from public sources alone whether our annual values are max-of-quarters or mean-of-quarters. That depends on whoever produced `MPA_Runs_new.csv`. Inspect the data-prep code in `06b_landsat_processing.R` (or upstream) and document the actual aggregator.

---

## 5. SBC MBON integrated kelp-forest datasets

The largest input files in our pipeline (~1.2 GB combined) are the SBC MBON harmonized synthesis datasets that pool field data from PISCO + KFM + LTER + USGS San Nicolas Island monitoring programs.

### 5.1 Quad and swath cover

- **EDI package:** `edi.6.3` (or newer revision matching 2023-10-22 file stamp)
- **DOI:** https://doi.org/10.6073/pasta/85a8f46dabe413973195c8530911088a
- **Metadata viewer:** https://portal.edirepository.org/nis/metadataviewer?packageid=edi.6.3
- **Citation:**
  > Miller, R.J., Rassweiler, A.R., Caselle, J.E., Kushner, D., Reed, D.C., Lafferty, K.D., Kui, L., & O'Brien, M. (2021). *Santa Barbara Channel Marine BON: Nearshore kelp forest integrated quad and swath cover, 1980-ongoing* (ver. 6.3). Environmental Data Initiative. https://doi.org/10.6073/pasta/85a8f46dabe413973195c8530911088a
- **Coverage:** 1980-present, ~100 sites in the Santa Barbara Channel and Channel Islands, harmonized across SBC LTER, PISCO, NPS Channel Islands KFM, and USGS San Nicolas Island.

### 5.2 Fish

- **EDI package:** `edi.5.3`
- **DOI:** https://doi.org/10.6073/pasta/0976f9969ef2789b77b659e9337c4c0c
- **Metadata viewer:** https://portal.edirepository.org/nis/metadataviewer?packageid=edi.5.3
- **Citation:**
  > Miller, R.J., Rassweiler, A.R., Caselle, J.E., Kushner, D., Reed, D.C., Lafferty, K.D., Kui, L., & O'Brien, M. (2018). *Santa Barbara Channel Marine BON: Nearshore kelp forest integrated fish, 1981-ongoing* (ver. 5.3). Environmental Data Initiative. https://doi.org/10.6073/pasta/0976f9969ef2789b77b659e9337c4c0c
- **Coverage:** 1981-present. The companion site geolocation file ships inside this package.

---

## 6. PISCO data-deposit DOI (separate from the data paper)

In addition to the Malone et al. 2022 *Ecology* data paper, the actual data deposit is at PISCO Metacat / KNB:

- **DOI:** https://doi.org/10.6085/AA/PISCO_kelpforest.1.6
- **Landing page:** https://data.piscoweb.org/metacatui/#view/doi:10.6085/AA/PISCO_kelpforest.1.6
- **Citation:**
  > PISCO. *PISCO: Subtidal: Community Surveys: Kelp Forest* (ver. 1.6). PISCO Metacat. https://doi.org/10.6085/AA/PISCO_kelpforest.1.6
- **Note:** The 2024 file refresh in our pipeline is more recent than the v1.6 archived deposit. Cite both Malone et al. 2022 + this DOI, and note in Methods that we use a 2024 PISCO data pull (post-publication update from the PISCO data team).

---

## 7. CDFW MPA boundary dataset

Used to derive the `MPAfeatures_subset.csv` table that supplies MPA hectares (a meta-analysis moderator) and MPA designations (SMR / SMCA).

- **Citation:**
  > California Department of Fish and Wildlife (2019). *California Marine Protected Areas* [ds582]. CDFW BIOS / California Open Data Portal. https://data.ca.gov/dataset/california-marine-protected-areas-ds582
- **BIOS metadata:** https://map.dfg.ca.gov/metadata/ds0582.html
- **License:** CC BY 4.0
- **Version:** Reflects the MPA Network as of 1 January 2019 (post-MLPA buildout). Last open-data refresh 2023-03-01.

---

## 8. Summary table: citation stack for Methods section

For the manuscript Methods, this is the recommended minimal citation stack:

**Monitoring programs (cite all four data papers + protocol references)**
- PISCO: Malone et al. 2022 *Ecology* (data paper, https://doi.org/10.1002/ecy.3630) + PISCO data deposit DOI (10.6085/AA/PISCO_kelpforest.1.6) + Caselle et al. 2015 *Sci Rep* + Eisaguirre et al. 2020 *Ecology* + Rassweiler et al. 2020 *Front Mar Sci* (vertical-level structure) + Hamilton & Caselle 2015 *Proc R Soc B* (urchin sizing context).
- KFM: Davis et al. 1997 *NPS Sampling Protocol Handbook* (primary citation for all KFM protocol details).
- LTER: Reed & Miller, EDI packages `knb-lter-sbc.18` (kelp), `knb-lter-sbc.17` (fish), `knb-lter-sbc.77` (lobster). SBC LTER protocol PDFs from sbclter.msi.ucsb.edu.
- SBCMBON: Miller et al. 2018 (`edi.5.3` fish) + Miller et al. 2021 (`edi.6.3` quad/swath), the harmonization layer that the pipeline reads from.

**MPA boundaries**
- CDFW 2019 ds582, California Marine Protected Areas dataset (https://data.ca.gov/dataset/california-marine-protected-areas-ds582).

**Allometric biomass conversions**
- Macrocystis: Rassweiler et al. 2018 *Ecology* (`ecy.2440`).
- Sea urchins: Ebert 2010 *Sea Urchins: Biology and Ecology* (Academic Press), as calibrated for SBC LTER.
- Lobster: Barsky 2001 *CDFG Lobster Assessment* (allometric), as calibrated for SBC LTER.
- Sheephead: Cowen 1990 *Environmental Biology of Fishes*, as calibrated for SBC LTER.
- Fish length-weight coefficients: FishBase (Froese & Pauly, www.fishbase.org), compiled in the PISCO species attribute table.

**Landsat kelp**
- Cavanaugh et al. 2011 *MEPS* (calibration) + Bell et al. 2020 *RSE* (multi-sensor harmonization) + Bell et al. 2023 *PLOS ONE* (Kelpwatch portal). Underlying SBC LTER package: `knb-lter-sbc.74`.

**Statistical framework**
- pBACIPS: Thiault, L., Kernaléguen, L., Osenberg, C.W. & Claudet, J. (2017). Progressive-Change BACIPS: a flexible approach for environmental impact assessment. *Methods in Ecology and Evolution* 8(3): 288-296. https://doi.org/10.1111/2041-210X.12655 *(published online 2016; print issue 2017)*

**MPA governance**
- CDFW 2023 Decadal Management Review.
- Carr et al. 2022 PISCO contributions report.
- Carr et al. 2021 OPC Kelp Forest Technical Report (formalizes the academic-program protocols including VRG roving-diver lobster sizing).

---

## Open items lookup status (updated 2026-05-04 PM)

### Item 2 RESOLVED: MBON `area == 10` for KFM Macrocystis

`area == 10` represents the **paired-diver 10 m² convention**: each row in the MBON-integrated CSV is one diver's 5 m × 1 m = 5 m² half-quad, and the two paired divers (sampling opposite sides of the same transect position) together cover 10 m². The Davis et al. 1997 KFM Handbook (p. 13, 38) describes the protocol explicitly:
> "5m Quadrat ... 1×5 m ... 40 per site"
> "The 100 m transect is divided into 20 quadrats 5 m in length and 1 m wide... Each diver samples opposite sides of the transect resulting in 40 quadrats overall."

So our pipeline filter `area == 10` for KFM Macrocystis is **correct**. It selects the modern (1996+) Macrocystis-specific 5m quadrat protocol and excludes the legacy 1m and 2m quad records. Possibility B (5×2 m quad) was wrong. Possibility A (paired-diver convention) is right. Density-per-m² calculations correctly divide by `area = 10` (the paired total).

Action: add an inline comment in the data-processing repo: `area = 10` reflects MBON's paired-diver 10 m² convention for the KFM 5 m quadrat protocol. See Davis et al. 1997 KFM Handbook Vol 1 pp. 13, 38.

### Item 3 RESOLVED: `MPA_Runs_new.csv` aggregator is MAX-of-quarters

The annual values in `MPA_Runs_new.csv` are **annual MAX (peak) kelp canopy biomass** per MPA polygon, produced by Tom Bell from the SBC LTER kelp-watch product (`knb-lter-sbc.74`). Verified via three independent in-repo documentation sources:

1. **Original Donham archive code** (`code/archive/pBACIPS_PISCO_V10.R` L5854, the foundational Donham/Stier script):
 > `fn <- read.csv("data/LANDSAT/MPA_Runs_new.csv") #Import TBell spreadsheet of annual maximum biomass`

2. **Current production code header** (`code/R/06b_landsat_processing.R` L60–64): "Each annual value is the MAX (peak) kelp canopy biomass observed across the year..."

3. **Harmonized data dictionary** (L106): "Annual values are MAX (peak) kelp canopy biomass within the year."

Consistent with Bell et al. 2023 *PLOS ONE* (Kelpwatch portal) convention: *"the maximum canopy area was determined for each cell for each year from 1984 to 2021."* No prep script ships alongside the CSV. Tom Bell would need to be the source if we wanted to reconstruct the polygon extraction. For absolute verification we'd download the SBC LTER NetCDF and reproduce `max(quarterly)` over each MPA polygon (~one afternoon of work, not blocking).

### Item 1 PARTIALLY RESOLVED: Rassweiler 2018 calibration window CORRECTION

**Critical finding: the "May–Oct only" claim is NOT supported by the public record.** The publicly available SBC LTER Macroalgal NPP protocol (https://sbclter.msi.ucsb.edu/external/Reef/Protocols/Seasonal_Benthic_Survey/Seasonal_Benthic_Survey_Protocol_Net_Primary_Production_of_macroalgae.pdf) states unambiguously that they apply *"month-specific relationships between frond density (no. m⁻²) and dry mass density (dry kg m⁻²) developed by Rassweiler et al. (2018)"*, i.e., **12 monthly slopes**, not a May–Oct subset. There is a separate "summer biomass → annual NPP" regression in the same protocol that IS summer-only, but it is for a different purpose (predicting annual NPP from summer standing crop), not for converting density to biomass.

The "May–Oct" attribution for `MACRO_AVE_SLOPE` originates from Emily's archived script (`pBACIPS_PISCO_V10.R` L674: `#Slopes from May-Oct from LTER`). That may have been her own per-site averaging across the May–Oct months from the underlying monthly table, or it may have been a verbal note from someone at LTER, but it is not the LTER's published convention. The values themselves (0.10386, 0.10103, 0.09267, 0.09204, 0.08054, 0.08505) cannot be matched to any open-access table. They may come from an unpublished SBC LTER calibration table, the paywalled Rassweiler 2018 *Ecology* paper, or its EDI supplement.

**Action items:**
1. **Code (already done):** the `MACRO_AVE_SLOPE` docstring in `01_utils.R` was rewritten to remove the unverified "May–Oct only" assertion and to point readers to the SBC LTER protocols + Rassweiler 2018.
2. **Methods text:** the safer provisional Methods sentence is *"Macrocystis density was converted to biomass using the average of site-specific frond-density-to-biomass slopes from the SBC LTER kelp calibration (Rassweiler et al. 2018; mean slope = 0.092 kg dry mass per stipe)"*, without asserting May–Oct.
3. **Final verification (Emily action):** open Rassweiler et al. 2018 *Ecology* 99(9):2132 via UCSB Wiley access, locate the per-site / per-month calibration table, and confirm the six values. The Wiley page returned 403 to programmatic fetch and the supplement is not in the public protocol PDFs.

### Item 4 RESOLVED: VRG roving-diver lobster sizing protocol

**Primary citation (canonical, formal MPA Monitoring Program technical report):**

> Carr, M.H., Caselle, J.E., Cavanaugh, K., Freiwald, J., Kroeker, K., Pondella, D., Tissot, B., Malone, D., Parsons-Field, A., & Spiecker, B. (2021). *Monitoring and Evaluation of Kelp Forest Ecosystems in the MLPA Marine Protected Area Network.* Report submitted to California Sea Grant, Ocean Protection Council MPA Monitoring Program, and California Department of Fish and Wildlife (Marine Resources Division), 30 December 2021.
> URL: https://caseagrant.ucsd.edu/system/files/2022-06/Kelp%20Forest%20Technical%20Report%20Narrative_v2.pdf (verified live, 32 MB PDF)

Daniel Pondella (VRG/Occidental) is a named PI on this report. The protocol applies to all four academic monitoring groups (UCSB, UCSC, VRG, HSU). Verbatim quote from the narrative (p. 20–21):

> *"Academic programs — To characterize the ecological community and geological features at each site, we conduct four types of diver surveys: 1) density and size distribution of all conspicuous fishes... 2) density of large invertebrates and stipitate algae... 3) percent cover of sessile invertebrates... and 4) **size frequency for the commercially and ecologically important invertebrates and algae such as red and purple urchins, abalone, lobsters, giant kelp, and other key species.**"*
>
> *"Abalones (Haliotis spp.), red and purple sea urchins... **California spiny lobsters (Panulirus interruptus)**, select species of sea stars, and some whelks **are sized to the nearest centimeter. These data are collected from the benthic swaths or haphazard surveys across depth zones of each survey site.**"*

The "haphazard surveys across depth zones" language is exactly the off-transect roving-diver sizing that distinguishes VRG from the UCSB direct-on-swath protocol, and the wording "**benthic swaths OR haphazard surveys**" formally documents that both modes coexist.

**Backup citation (peer-reviewed data paper):**

> Malone, D.P., Davis, K., Lonhart, S.I., Parsons-Field, A., Caselle, J.E., & Carr, M.H. (2022). Large-scale, multidecade monitoring data from kelp forest ecosystems in California and Oregon (USA). *Ecology* 103(5): e3630. https://doi.org/10.1002/ecy.3630

The Carr et al. 2021 report references this paper as the formal protocol reference ("Malone et al. In Press" → Malone et al. 2022). Use this as the journal-required peer-reviewed citation.

**Backup 2 (VRG-specific biogeographic methods paper):**

> Pondella, D.J. II, Williams, J.P., Claisse, J., Schaffner, R., Ritter, K., & Schiff, K. (2015). The Physical Characteristics of Nearshore Rocky Reefs in the Southern California Bight. *Bulletin of the Southern California Academy of Sciences* 114(3): 105–122.

**Note:** No standalone VRG-only methods paper documents the roving-diver lobster sizing protocol in isolation. The protocol is documented as a shared "Academic program" protocol (UCSB/UCSC/VRG/HSU) in Carr et al. 2021 and Malone et al. 2022. The VRG Dive Program web page (https://www.oxy.edu/academics/vantuna-research-group/dive-program) confirms VRG collects size frequency but does not provide a standalone citable methods section.

---

## Final lookup status: all four items closed

| # | Item | Status | Resolution |
|---|---|---|---|
| 1 | Rassweiler 2018 calibration window | PARTIAL | "May–Oct only" claim corrected to "month-specific" per public LTER protocols. Code docstring revised. Final per-site value verification still requires UCSB Wiley access to Rassweiler 2018 supplement. |
| 2 | MBON `area == 10` for KFM Macro | DONE | Paired-diver 10 m² convention; Davis 1997 KFM Handbook p.13, 38. Filter is correct. |
| 3 | `MPA_Runs_new.csv` aggregator | DONE | Annual MAX-of-quarters per Tom Bell's polygon extraction over knb-lter-sbc.74. Verified via 3 independent in-repo sources (Donham archive code, current 06b header, DATA_DICTIONARY). |
| 4 | VRG roving-diver lobster sizing | DONE | Carr et al. 2021 OPC tech report + Malone et al. 2022 *Ecology*. Verbatim "benthic swaths or haphazard surveys across depth zones" language formalizes the protocol. |

---

## Appendix A: Citation & URL verification status (2026-05-04)

All 26 unique URLs in this document were programmatically HEAD-checked. All 5 DOIs that returned publisher 403/401 anti-bot codes were independently verified through the CrossRef API (which doesn't gate). Results:

### Live and direct-accessible (21 URLs, HTTP 200)

| Source | Status |
|---|---|
| https://coastwatch.pfeg.noaa.gov/erddap/tabledap/erdCinpKfmSFNH.graph | OK 200 |
| https://doi.org/10.1016/j.rse.2018.06.039 (Bell et al. 2020 RSE) | OK 200 |
| https://doi.org/10.1038/srep14102 (Caselle et al. 2015) | OK 200 |
| https://doi.org/10.1098/rspb.2014.1817 (Hamilton & Caselle 2015) | OK 200 |
| https://doi.org/10.1371/journal.pone.0271477 (Bell et al. 2023) | OK 200 |
| https://doi.org/10.3389/fmars.2020.00272 (Rassweiler et al. 2020) | OK 200 |
| https://irma.nps.gov/DataStore/DownloadFile/485444 (Davis et al. 1997 KFM protocol; 459 KB PDF, 96 pp.) | OK 200 (manually verified content; quad-size table at p.35 confirmed verbatim) |
| https://nrm.dfg.ca.gov/FileHandler.ashx?DocumentID=203484 (Carr et al. 2022 PISCO MLPA report) | OK 200 |
| https://nrm.dfg.ca.gov/FileHandler.ashx?DocumentID=209209&inline (CDFW 2023 Decadal Review) | OK 200 |
| https://pmc.ncbi.nlm.nih.gov/articles/PMC12231823/ (Spiecker et al. 2025) | OK 200 |
| https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-sbc.74.18 | OK 200; title verified: *"SBC LTER: Time series of quarterly NetCDF files of kelp biomass in the canopy from Landsat 5, 7 and 8, since 1984 (ongoing)"*. Confirms quarterly aggregation. |
| https://portal.edirepository.org/nis/mapbrowse?scope=knb-lter-sbc&identifier=112 | OK 200 |
| https://portal.edirepository.org/nis/metadataviewer?packageid=knb-lter-sbc.17.29 | OK 200 |
| https://portal.edirepository.org/nis/metadataviewer?packageid=knb-lter-sbc.18.29 | OK 200 |
| https://portal.edirepository.org/nis/metadataviewer?packageid=knb-lter-sbc.77.10 | OK 200 |
| https://sbclter.msi.ucsb.edu/data/catalog/package/?package=knb-lter-sbc.74 | OK 200 |
| https://sbclter.msi.ucsb.edu/external/Reef/Protocols/Kelp_Forest_Community_Dynamics/SBC_LTER_protocol_Reed_Kelp_forest_community_Density_giant_kelp_20130524.pdf | OK 200 |
| https://sbclter.msi.ucsb.edu/external/Reef/Protocols/Kelp_NPP/KelpNPP_20180522.pdf | OK 200 |
| https://wildlife.ca.gov/Conservation/Marine/MPAs/Management/Decadal-Review | OK 200 |
| https://www.oxy.edu/academics/vantuna-research-group/dive-program | OK 200 |

### Anti-bot gated by publisher, but DOIs verified through CrossRef (5 URLs)

These DOIs return HTTP 403 to `curl` (publisher anti-scraping) but resolve correctly in any browser. Each was verified via the CrossRef API (`https://api.crossref.org/works/<doi>`):

| DOI | CrossRef-verified title and venue | Year |
|---|---|---|
| 10.1002/ecy.2440 | "Improved estimates of net primary production, growth, and standing crop of *Macrocystis pyrifera* in Southern California" — *Ecology* | 2018 |
| 10.1002/ecy.2993 | "Trophic redundancy and predator size class structure drive differences in kelp forest ecosystem dynamics" — *Ecology* | 2020 |
| 10.1002/ecy.3630 | "Large-scale, multidecade monitoring data from kelp forest ecosystems in California and Oregon (USA)" — *Ecology* | 2022 |
| 10.1111/2041-210X.12655 | "Progressive-Change BACIPS: a flexible approach for environmental impact assessment" — *Methods in Ecology and Evolution* | 2016 (online); 2017 (print) |
| 10.3354/meps09141 | "Environmental controls of giant-kelp biomass in the Santa Barbara Channel, California" — *Marine Ecology Progress Series* | 2011 |

**Bottom line:** All 26 URLs are live, and all DOIs resolve to the correct cited paper. The Davis et al. 1997 KFM handbook PDF was independently downloaded and the key methodology quote (quad-size history, p.35) was verified verbatim against the source.

### Verification still pending (manual lookup recommended before submission)

| Item | Reason |
|---|---|
| Rassweiler et al. 2018 *Ecology* (`ecy.2440`) supplement (May–Oct calibration window) | Supplement is paywalled to non-Wiley subscribers, but verifiable via UCSB library access. |
| `MPA_Runs_new.csv` aggregator (max-of-quarters vs mean-of-quarters) | Internal. Depends on whoever generated the file from the SBC LTER quarterly NetCDFs. Inspect `06b_landsat_processing.R` and any upstream prep scripts. |
| MBON `area == 10` for KFM Macrocystis | The MBON-integrated CSV's metadata (header or accompanying README) should clarify whether this represents two paired 5 m² quads summed. |
| VRG roving-diver lobster sizing protocol document | Pondella et al. 2015 *Bull. SC Acad. Sci.* + Oxy VRG dive-program page are the closest available. A formal tech report would be stronger for Methods. |
