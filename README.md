# Monterey Bay Wildfires â€“ Phytoplankton Community Shifts

**Capone, D., Daniel, P., Kudela, R., Barton, A.D., Kahru, M., & DÃ©cima, M. (in review).**  
*Extreme wildfire conditions shift coastal phytoplankton community structure in California.*  
*Limnology & Oceanography.*

---

## Overview

<<<<<<< HEAD
- data/                     # Raw input datasets (e.g., IFCB, satellite, shore station, aerosols)

- processed_data/            # Processed datasets generated for analysis

- scripts/                   # Analysis and figure generation scripts (MATLAB & R)

- Monterey-Bay-Wildfires-Shift/ # Manuscript text, figures, and supporting information
=======
This repository contains the code, data access links, and processed outputs associated with the study examining how the 2020 California Lightning Complex Fires influenced coastal phytoplankton communities in Monterey Bay.  
The analysis integrates satellite observations, in situ shore station and atmospheric monitoring data, and high-resolution Imaging FlowCytobot (IFCB) time-series to identify wildfire-linked shifts in phytoplankton size and taxonomic composition.

---

## Repository Structure

```
data/              # Links and metadata for all input datasets (see below)
processed_data/    # Derived and cleaned data products used for figures and models
scripts/           # MATLAB and R scripts for analysis, statistics, and figure generation
manuscript/        # Text, figures, and Supporting Information files for publication
```
>>>>>>> 14b88b4 (Docs: replace README with README_v2 and add Key Scripts; track environmental_GAMs, figure_6_NMDS, preprocessing and satellite chl scripts)

---

## Data Sources

All datasets used in this analysis are publicly available.  
For transparency and reproducibility, each dataset is linked to its original hosting source or DOI.

| Dataset | Description | Source / Access |
|----------|--------------|----------------|
| **Imaging FlowCytobot (IFCB)** | Automated imaging of phytoplankton community composition, abundance, and biovolume at the Santa Cruz Municipal Wharf (20 min intervals, 2015â€“present). | https://ifcb.caloos.org/timeline?dataset=santa-cruz-municipal-wharf |
| **Satellite Chlorophyll-a & SST** | Merged multi-sensor ocean color and temperature data for Monterey Bay (300 m, daily/5-day, 2002â€“present). | https://spg-satdata.ucsd.edu/MBay/MBay.htm |
| **CalHABMAP** | Weekly chlorophyll-a, dissolved nutrients, and harmful algal bloom (HAB) taxa from Santa Cruz Municipal Wharf. | https://calhabmap.org/santa-cruz-wharf |
| **Shore Stations (SCMW & MLML)** | Continuous seawater and meteorological observations including temperature, salinity, nitrate, and PAR. | SCMW: https://www.cencoos.org/observations/sensor-platforms/shore-stations/  â€¢  MLML: http://pubdata.mlml.calstate.edu/seawater/index.php |
| **AirNow PM2.5** | Daily ground-level particulate matter concentrations (â‰¤2.5 Âµm) from Santa Cruz AMS and San Lorenzo Valley stations. | https://www.epa.gov/outdoor-air-quality-data/download-daily-data |
| **AERONET** | Aerosol Optical Depth (340â€“1640 nm) for the Monterey site, Level 1.5â€“2.0, NASA-PHOTONS network. | https://aeronet.gsfc.nasa.gov/new_web/data.html |
| **MERRA-2 Black Carbon Reanalysis** | NASA reanalysis (0.5Â° Ã— 0.65Â°) of global black carbon aerosols (1980â€“present). | https://doi.org/10.5067/KLICLTZ8EM9D |
| **USGS San Lorenzo River Outflow** | Continuous discharge data from USGS gauge 11161000. | https://waterdata.usgs.gov/monitoring-location/11161000 |
| **NOAA NDBC Buoy 46042** | Wind speed and direction observations off Monterey Bay (1987â€“present). | https://www.ndbc.noaa.gov/station_page.php?station=46042 |
| **BEUTI Upwelling Index** | Daily Biologically Effective Upwelling Transport Index at 37Â°N, estimating vertical nitrate flux. | https://oceanview.pfeg.noaa.gov/products/upwelling/beuti |
| **HYSPLIT Trajectories** | NOAA Hybrid Single-Particle Lagrangian Integrated Trajectory (HYSPLIT) model runs for CZU, SCU, and August Complex Fires. | https://ready.arl.noaa.gov/HYSPLIT.php |

---

## Software Requirements

**MATLAB:** R2023bâ€“R2024a  
Used for: time-series processing, climatology computation, and cross-correlation analysis.  

**R:** 4.4.1 or higher  
Required packages:  
- `mgcv` â€” Generalized Additive Models (Wood, 2017)  
- `vegan` â€” NMDS and community analysis  
- `tidyverse` â€” Data manipulation and visualization  
- `patchwork` â€” Figure assembly  
- `here`, `readr`, `lubridate` â€” Reproducible file and date handling  

---

## Reproducibility Workflow

1. **Clone the repository:**
   ```bash
   git clone https://github.com/DanteCapone/MB_Wildfire_Obs.git
   cd MB_Wildfire_Obs
   ```

2. **Access and organize data:**  
   - Run `scripts/data_download.R` to automatically fetch open datasets and build `data/` directory.  
   - For any datasets requiring manual download, follow URLs in the table above.

3. **Process IFCB and environmental data:**  
   - MATLAB script: `scripts/ifcb_processing.m`  
   - R script: `scripts/environment_merge.R`

4. **Run analyses:**  
   - Cross-correlation analysis: `scripts/cross_correlation_analysis.m`  
   - Generalized Additive Models (GAMs): `scripts/gam_analysis.R`  

5. **Reproduce all figures:**  
   - `scripts/figure_generation.R` outputs publication-ready plots to `/manuscript/figures/`.

---

## Citation

If you use or adapt materials from this repository, please cite:

> Capone, D., Daniel, P., Kudela, R., Barton, A.D., Kahru, M., & DÃ©cima, M. (in review).  
> *Extreme wildfire conditions shift coastal phytoplankton community structure in California.*  
> *Limnology & Oceanography.*  
> Available at: https://github.com/DanteCapone/MB_Wildfire_Obs

---

## Funding and Acknowledgements

This work was supported by:  
- **NSF Graduate Research Fellowship Program** (No. DGE-2038238)  
- **NOAA Marine Sensors Transitions Program** (NA14NOS0120148)  
- **CeNCOOS / IOOS** (NA16NOS0120021)  
- **NOAA Saltonstall-Kennedy Program** (NA16NMF4270263)  

Data provided by CeNCOOS, MLML, NASA, NOAA, USGS, EPA, and UC Santa Cruz IFCB Network.  
We acknowledge the use of ChatGPT (versions 3.5â€“5) for code generation and editorial refinement.

---

## Key Scripts

- scripts/main_analysis/MB_Wildfire_Obs_main_analysis.m — Orchestrates the main analysis pipeline and figure generation.
- scripts/figure_6_NMDS.m — NMDS analysis and figure panels comparing community composition during key wildfire-adjacent periods.
- scripts/environmental_GAMs.R — Fits environmental Generalized Additive Models linking drivers (e.g., aerosols, upwelling, river outflow) to phytoplankton metrics; outputs diagnostics and tables.
- scripts/pre-processing/ifcb_compile_all.m — Compiles IFCB time series across years and taxa into analysis-ready tables.
- scripts/pre-processing/ifcb_compile.m — Helper for compiling subsets of IFCB data.
- scripts/pre-processing/slr_outflow_data_retrieval.R — Retrieves and cleans USGS San Lorenzo River outflow data.
- scripts/satellite_chl_analysis/MB_chla_analysis.m — Computes chlorophyll-a anomalies/climatologies and key visualizations for Monterey Bay.
- scripts/satellite_chl_analysis/MB_chla_analysis_pt1_explore_and_vis.m — Exploratory visualizations for satellite chlorophyll-a.
- scripts/satellite_chl_analysis/anomaly_climatology_analysis.m — Anomaly/climatology computation utilities.
- scripts/satellite_chl_analysis/SST_compile.m — Compiles SST time series and fields.
- scripts/satellite_chl_analysis/iop_analysis.m — Inherent optical properties analysis and plots.

---

## Key Scripts

- scripts/main_analysis/MB_Wildfire_Obs_main_analysis.m — Orchestrates the main analysis pipeline and figure generation.
- scripts/figure_6_NMDS.m — NMDS analysis and figure panels comparing community composition during key wildfire-adjacent periods.
- scripts/environmental_GAMs.R — Fits environmental Generalized Additive Models linking drivers (e.g., aerosols, upwelling, river outflow) to phytoplankton metrics; outputs diagnostics and tables.
- scripts/pre-processing/ifcb_compile_all.m — Compiles IFCB time series across years and taxa into analysis-ready tables.
- scripts/pre-processing/ifcb_compile.m — Helper for compiling subsets of IFCB data.
- scripts/pre-processing/slr_outflow_data_retrieval.R — Retrieves and cleans USGS San Lorenzo River outflow data.
- scripts/satellite_chl_analysis/MB_chla_analysis.m — Computes chlorophyll-a anomalies/climatologies and key visualizations for Monterey Bay.
- scripts/satellite_chl_analysis/MB_chla_analysis_pt1_explore_and_vis.m — Exploratory visualizations for satellite chlorophyll-a.
- scripts/satellite_chl_analysis/anomaly_climatology_analysis.m — Anomaly/climatology computation utilities.
- scripts/satellite_chl_analysis/SST_compile.m — Compiles SST time series and fields.
- scripts/satellite_chl_analysis/iop_analysis.m — Inherent optical properties analysis and plots.
