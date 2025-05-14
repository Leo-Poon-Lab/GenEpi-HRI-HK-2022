# Genomic epidemiology of hospital related infections in Hong Kong 2022

This project focuses on analyzing hospital-related infections in Hong Kong, including case classifications, genomic surveillance, and mobility data. The repository contains scripts, data, and results for visualizing and interpreting the findings.

## Project Structure

### 1. **Data**
The `data/` directory contains raw and processed data files used in the analysis. Key subdirectories include:
- **`case_curve/`**: Contains JSON and Excel files for aggregated case data.
  - `HK_case_data.json`: JSON file with aggregated confirmed case classifications.
  - `hk_epicurve.xlsx`: Processed case data in Excel format.
- **`hospital_data/`**: Contains metadata for hospitals and districts.
  - `hong_kong_hospitals_with_districts.csv`: Metadata for hospitals, including names, clusters, and geographic coordinates.
  - `Hong_Kong_Districts_with_HA_Clusters.csv`: Metadata for districts and their corresponding hospital clusters.
  - `metadata_mt.xlsx`: Metadata for hospital-acquired infections.
- **`mobility_data/`**: Contains mobility and transport data.
  - `2022_HK_Region_Mobility_Report.csv`: Google Community Mobility Report for Hong Kong.
  - `statistics_on_daily_passenger_traffic.csv`: Cross-border passenger traffic data.
  - `public_transport.xlsx`: Local transport data.

### 2. **Scripts**
The `scripts/` directory contains R scripts for data processing and visualization. Key scripts include:
- **`GISAID_data_processing/QC.R`**:
  - Processes and filters GISAID genomic data.
  - Steps include running `pangolin` and `nextclade`, filtering low-quality samples, and exporting filtered metadata.
- **`1.metadata_visualization.R`**: 
  - Visualizes aggregated case classifications, genomic surveillance, and hospital-acquired infections.
  - Generates figures for case trends, genomic data, and hospital-related policies.
  - Figure 1 panels:
    - **Panel A**: Aggregated confirmed case classifications and sequenced samples.
    - **Panel B**: Hospital-acquired infections by date and ward.
    - **Panel C**: Mobility and hospital-related policy data.
    - **Panel D**: The map.
  - Figure 2 panel a:
  - **Panel A**: Officially reported hospital related infections in Hong Kong.
- **`2.tree_visualization.R`**:
  - Visualizes phylogenetic trees for genomic data.
  - Used for panels b and c for Figure 2.

### 3. **Results**
The `results/` directory contains output figures and processed data

## How to Run
1. **Install Required Libraries**:
   Ensure the following R libraries are installed:
   - `tidyverse`, `jsonlite`, `lubridate`, `ggplot2`, `ggforce`, `ggrepel`, `MetBrewer`, `ggbump`, `scales`, `sp`, `gridExtra`, `RColorBrewer`, `readxl`, `writexl`.

2. **Run Scripts**:
   - Execute `1.metadata_visualization.R` to generate Figures 1a, 1b, 1c, and 1d.
   - Execute `2.tree_visualization.R` to generate phylogenetic tree visualizations.
   - Execute `GISAID_data_processing/QC.R` to process and filter GISAID genomic data. (Details refer to [here](scripts/GISAID_data_processing/README.md))

3. **View Results**:
   - Figures are saved in the `results/` directory.

## Citation
Pending

---