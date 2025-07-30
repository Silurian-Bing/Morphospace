      
# Code and Data for "Evolutionary innovation underpins resilience to mass extinction"

This repository contains the data and R scripts required to reproduce all analyses and figures presented in the manuscript.

The analytical workflow is structured as a sequence of six R scripts, designed to be run in numerical order.

---

### 1. System Requirements

*   **Software**: R version 4.2 or newer.
*   **Operating Systems**: The code has been tested on Windows 10/11.
*   **R Packages**: The following R packages are required. They can be installed by running the command in the "Installation" section below.
    *   `readxl`
    *   `dplyr`
    *   `tidyr`
    *   `ggplot2`
    *   `vegan`
    *   `viridis`
    *   `mgcv`
    *   `gridExtra`
    *   `RColorBrewer`
    *   `ape`
    *   `phytools`
*   **Hardware**: No special hardware is needed. A standard desktop or laptop computer with at least 8 GB of RAM is sufficient.

---

### 2. Installation and Setup

1.  **Download Repository**: Clone or download this repository to your local machine.
2.  **Set Working Directory**: All scripts use relative paths and **must be run from the root directory of this repository**. Before running any script, open R or RStudio and set the working directory. For example:
    ```r
    # On Windows
    # setwd("C:/path/to/your/downloaded-repo-folder")
    # On macOS/Linux
    # setwd("/path/to/your/downloaded-repo-folder")
    ```
3.  **Install Packages**: If you do not have the required packages, run the following command in the R console:
    ```r
    install.packages(c("readxl", "dplyr", "tidyr", "ggplot2", "vegan", "viridis", "mgcv", "gridExtra", "RColorBrewer", "ape", "phytools"))
    ```

*   **Estimated Setup Time**: 5-10 minutes.

---

### 3. How to Reproduce the Analysis

To reproduce all results from the manuscript, run the scripts sequentially from `01` to `06`. Each script loads the necessary data and saves its outputs (figures, processed data) before the next script is run.

*   `01_PCA_and_Morphospace_Visualization.R`: Performs the primary Principal Component Analysis (PCA) and generates basic morphospace visualizations, including centroid trajectories.
*   `02_Diversity_Disparity_Comparison.R`: Calculates taxonomic richness (from `OtoD_atrypidesB516.xlsx`) and morphological disparity, and compares their trends over geological time.
*   `03_Disparity_Through_Time_Analysis.R`: Conducts the Disparity Through Time (DTT) analysis using geological time bins derived from the first and last appearance data.
*   `04_Environmental_Covariation_Analysis.R`: Integrates paleoenvironmental proxies (temperature and sea level) and performs cross-correlation analyses with the biological trends.
*   `05_Phylogenetic_Clustering_and_Signal.R`: Integrates the phylogenetic tree (`T29h tree.nex`), performs morphospace clustering, and quantifies phylogenetic signal (Blomberg's K).
*   `06_Phylomorphospace_Integration.R`: Generates the final integrated phylomorphospace visualizations, combining phylogeny, morphology, and temporal data.

*   **Expected Total Run Time**: Approximately 30-45 minutes on a standard desktop computer.
*   **Expected Output**: The scripts will generate all figures and statistical results reported in the manuscript. Figures will be saved in a `_outputs/figures/` sub-directory, and processed data tables will be saved in `_outputs/data/`.

---

### 4. Data Description

All necessary data files are located in the `_data/` directory:

*   `morphospace.xlsx`: The primary multi-sheet morphological character matrix for all genera.
*   `taxon_age_data.csv`: First and Last Appearance Datum (FAD/LAD) for each genus used for temporal binning.
*   `OtoD_atrypidesB516.xlsx`: The comprehensive occurrence dataset used for calculating taxonomic diversity.
*   `tem470_370.csv`: Paleo-temperature proxy data.
*   `sea470_370.csv`: Eustatic sea-level proxy data.
*   `T29h tree.nex`: The phylogenetic tree of the Atrypida in NEXUS format.

---

### 5. Citation





# Morphospace Evolution Analysis Pipeline

A comprehensive R-based analytical pipeline for quantifying morphological evolution patterns through geological time, integrating multiple approaches from basic morphospace visualization to phylogenetically-informed analyses and environmental correlations.

Pipeline Overview
This collection of six R scripts provides a complete workflow for morphospace evolution analysis, designed with increasing complexity and analytical sophistication. Each script builds upon previous analyses while maintaining independence for modular usage.
Script Execution Order
## 1. morphospace_evolution_analysis.R
Foundation: Basic Morphospace Visualization
The entry point for morphospace analysis, performing Principal Component Analysis (PCA) of morphological datasets across geological time periods.
Core Functions:
•	Strict Excel data validation and processing (A:AV column enforcement)
•	Multi-period morphological data integration with metadata preservation
•	Standardized PCA implementation with comprehensive error handling
•	Advanced visualization including point clouds, confidence ellipses, and temporal trajectories
Key Outputs:
•	Basic morphospace plots with time-period color coding
•	Confidence ellipse representations of morphological variance
•	Centroid trajectory analysis showing evolutionary pathways
•	Character loadings analysis for morphological drivers identification
________________________________________
## 2. morphospace_diversity_comparison.R
Expansion: Diversity-Disparity Integration
Comparative analysis between morphological disparity (morphospace occupancy) and taxonomic diversity through geological time, incorporating multiple diversity estimation methods.
Core Functions:
•	Multi-method diversity estimation (Chao1, SQS, Combined estimators)
•	Bootstrap confidence interval generation for all metrics
•	PCA-based morphological disparity quantification
•	GAM modeling for trend detection and relationship quantification
Key Outputs:
•	Diversity method comparison plots (Chao1 vs SQS vs Combined)
•	Dual-axis temporal plots showing diversity-disparity relationships
•	Standardized multi-variate time series analysis
•	Comprehensive correlation matrices between variables
________________________________________
## 3. dtt_timebin_analysis.R
Temporal Patterns: Disparity Through Time Analysis
Time-binned Disparity Through Time (DTT) analysis using geological time windows rather than phylogenetic branch lengths, providing robust macroevolutionary pattern detection.
Core Functions:
•	Automated geological time bin assignment based on FAD/LAD data
•	Multi-metric disparity calculation (variance, range, distance-based measures)
•	Temporal trend detection with linear regression analysis
•	Bootstrap confidence intervals for disparity estimates
Key Outputs:
•	Total morphological variance evolution through time
•	PC-axis specific disparity trajectories
•	Taxonomic diversity curves with genus-level richness
•	Statistical trend analysis with significance testing
________________________________________
## 4. dtt_environment_covariation.R
Environmental Integration: Climate-Evolution Relationships
Quantification of covariation patterns between morphological evolution (DTT) and paleoenvironmental changes, integrating high-resolution climate proxy data.
Core Functions:
•	GAM smoothing framework for noise reduction across time series
•	Environmental data integration (temperature, sea level) with temporal alignment
•	Cross-correlation analysis between environmental and evolutionary variables
•	High-resolution temporal synchronization (0.5 Ma intervals)
Key Outputs:
•	Comprehensive time series showing normalized environmental and evolutionary trajectories
•	Full correlation matrix heatmaps with statistical significance
•	Environmental-evolutionary covariation significance tests
•	Publication-ready multi-panel visualizations
________________________________________
## 5. phylo_morphospace_clustering.R
Phylogenetic Context: Tree-Independent Clustering
Phylogenetically-informed morphospace analysis using intelligent clustering methods to identify evolutionary patterns without requiring time-calibrated trees.
Core Functions:
•	Phylogenetic integration with morphometric datasets
•	Automated optimal cluster detection using silhouette analysis
•	Phylogenetic signal quantification (Blomberg's K) for morphological traits
•	Temporal pattern detection in morphospace occupation
Key Outputs:
•	Morphospace clustering with confidence ellipses and phylogenetic distribution
•	Temporal trajectory analysis with age-coded morphospace patterns
•	PC1 temporal evolution with group-specific trends
•	Phylogenetic tree visualization with morphological group assignments
________________________________________
## 6. treepca_phylomorphospace.R
Advanced Integration: Full Phylogenetic-Morphospace Analysis
Comprehensive analysis integrating phylogenetic relationships, morphological data, and temporal information for understanding macroevolutionary patterns with full phylogenetic context.
Core Functions:
•	Three-way data integration (phylogeny, morphology, temporal data)
•	Intelligent clustering with phylogenetic constraints
•	Phylogenetic signal analysis across morphospace axes
•	Comprehensive data validation and intersection analysis
Key Outputs:
•	Integrated phylogenetic-morphospace visualization
•	Morphological group characterization with phylogenetic context
•	Temporal disparity analysis with phylogenetic constraints
•	Comprehensive statistical reports with evolutionary interpretations

## Data Requirements
Essential Files:
•	morphospace.xlsx - Multi-sheet Excel file with morphological character matrices
•	taxon_age_data.csv - First/Last Appearance Datum data for temporal analysis
•	T29h tree.nex - Phylogenetic tree in NEXUS format (scripts 5-6)
•	tem470_370.csv - Temperature proxy data (script 4)
•	sea470_370.csv - Sea level data (script 4)
•	OtoD_atrypidesB516.xlsx - Diversity occurrence data (script 2)
Working Directory:
All scripts expect data files in: F:/HB/DataBase/Morphospace/
Installation and Usage
Prerequisites:

## Core packages required across all scripts
install.packages(c("readxl", "dplyr", "tidyr", "ggplot2", "vegan", "viridis", 
                   "mgcv", "gridExtra", "RColorBrewer", "ape", "phytools"))
Execution:
1.	Set working directory to data location
2.	Ensure all required data files are present
3.	Execute scripts in numerical order for complete pipeline
4.	Each script can be run independently for specific analyses
Output Structure
Each script generates:
•	High-resolution figures (PNG/PDF format, 300 DPI)
•	Processed datasets (CSV format for downstream analysis)
•	Statistical summaries (Text reports with comprehensive results)
•	R objects (RDS format for reproducibility)
Analytical Philosophy
This pipeline emphasizes:
•	Data integrity through comprehensive validation
•	Methodological rigor with bootstrap confidence intervals
•	Reproducibility through detailed documentation
•	Scalability with modular design for different datasets
•	Publication readiness with professional visualization standards
Citation and Usage
When using this pipeline, please acknowledge the analytical framework and cite relevant methodological papers for PCA, DTT analysis, phylogenetic comparative methods, and environmental correlation techniques as appropriate for your specific research application.

