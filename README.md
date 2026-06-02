# Salmonella AMR Prediction & AEI Analysis
## Project Overview
This repository is intended to support the study on antimicrobial resistance (AMR) prediction and Antibiotic Resistance Epidemic Index (AEI) evaluation of Salmonella enterica. All resources are provided to ensure complete reproducibility of our analyses, model training, statistical tests and visualization.

## Overview of Subdirectories
For detailed information, please refer to command.sh.
- **ARG_combine**: Reference databases and annotation scripts for antimicrobial resistance gene annotation.
- **RefSeq_Anno**: Raw data and computational codes for *Salmonella* RefSeq genome annotation, alongside scripts used to generate finalized gene annotation outputs.
- **ARG_MIC**: Scripts for dataset construction and XGBoost model training based solely on antimicrobial resistance genes.
- **ARGSNP**: Annotation pipelines, dataset generation codes and modeling scripts for SNP/INDEL-driven AMR prediction.
- **ARG_with_SNP**: Codes to construct integrated datasets and train prediction models combining resistance genes and SNP/INDEL variants.
> Within each above folder, the subfolder **RS_811** contains matched phenotypic information and feature matrices for corresponding modeling workflows.
- **MGE**: Source codes for mobile genetic element (MGE) annotation.
- **Feature_analysis**: Codes and supporting datasets for the calculation of the **Antibiotic Resistance Epidemic Index (AEI)**.
- **Plot**: All scripting codes used for figure generation and statistical visualization.

## Repository Contents
### 1. Genome Processing Pipeline
**Full pipeline** for genome processing, annotation and feature construction. 

### 2. Machine Learning Scripts (Python & R)
All scripts for:
- Model training
- Feature selection
- Hyperparameter optimization
- Cross-validation
- Feature importance analysis

### 3. Modeling Datasets
Processed **feature matrices and metadata tables** used for modeling across 13 antibiotics.

### 4. Trained XGBoost Models
**Complete trained model files**, along with detailed model architecture and parameter configurations for full reproduction of XGBoost models.
**Finalized XGBoost models** are deposited at: Feature_analysis/use_topfeature/XGBoost_XXX_topfeature_finalmodel_traindata.txt.
**Corresponding feature** importance results of each model: Feature_analysis/use_topfeature/XGBoost_XXX_topfeature_finalmodel_feature_importances.txt.
**Model application instructions and optimized hyperparameters** for all 13 antibiotics are documented in: Feature_analysis/13_XGBoost_topfeature_final_module_XXX_predict.py.

### 5. Analytical & Visualization Code
Scripts for:
- Calculation of the **AEI**
- Statistical analyses
- Figure generation and plotting

## Environment & Dependencies
Please install the required Python and R packages before running the scripts.
