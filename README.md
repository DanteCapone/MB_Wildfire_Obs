# Monterey Bay Wildfires Phytoplankton Shifts

This repository contains the code, data, and processed outputs associated with the manuscript:

Capone et al. (in review). Extreme wildfire conditions shift coastal phytoplankton community structure in California.

The project investigates how the 2020 California Lightning Complex Fires influenced phytoplankton communities in Monterey Bay using a combination of satellite, in situ, and Imaging FlowCytobot (IFCB) time-series data.

├── data/                     # Raw input datasets (e.g., IFCB, satellite, shore station, aerosols)

├── processed_data/            # Processed datasets generated for analysis

├── scripts/                   # Analysis and figure generation scripts (MATLAB & R)

├── Monterey-Bay-Wildfires-Shift/ # Manuscript text, figures, and supporting information

### Requirements
- **MATLAB R2023b–2024a**  
  For time-series analysis, climatology calculations, and cross-correlations.  
- **R 4.4.1** with the following packages:  
  - `mgcv` (Generalized Additive Models)  
  - `vegan` (NMDS and clustering)  
  - `tidyverse` (data manipulation/plotting)  
  - `patchwork` (figure assembly)  

### Suggested Workflow
1. Clone the repository:  
   ```bash
   git clone https://github.com/<your-username>/Monterey-Bay-Wildfires-Shift.git
   cd MB_Wildfire_Obs

### Citation

If you use this repository, please cite:

Capone, D., Daniel, P., Kudela, R., Barton, A.D., Kahru, M., & Décima, M. (in review). Extreme wildfire conditions shift coastal phytoplankton community structure in California. Limnology & Oceanography.
