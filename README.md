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

