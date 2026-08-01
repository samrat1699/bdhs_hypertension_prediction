# BDHS 2022 Hypertension Analysis & Survey-Aware Machine Learning

## 📌 Project Overview
This repository contains a comprehensive statistical and machine learning pipeline for analyzing hypertension prevalence, its socio-demographic determinants, and developing predictive models using the **Bangladesh Demographic and Health Survey (BDHS) 2022** dataset. 

The project bridges traditional epidemiological survey analysis (using R) with advanced, survey-aware machine learning techniques (using Python) to provide robust, publication-ready insights, geospatial visualizations, and predictive tools.

## 📂 Repository Structure
- **`main_1.R`, `2nd_main.R`, `3d.R`**: R scripts for complex survey design, weighted statistical modeling (Logistic Regression), and visualizations (Forest Plots, Geospatial Maps).
- **`Hypertension_universityProject.ipynb`, `Final_hp.ipynb`**: Python notebooks for data preprocessing, survey-aware machine learning modeling, probability calibration, and SHAP-based interpretability.
- **`BDHS_2022_HYPERTENSION_SELECTED_VARIABLES.csv`**: The cleaned and processed dataset used for analysis.

## 📊 Dataset
- **Source**: Bangladesh Demographic and Health Survey (BDHS) 2022.
- **Target Variable**: `hypertension` (Binary: 1 = Yes, 0 = No).
- **Features**: 
  - *Socio-demographics*: Age, Sex, Education, Marital Status, Wealth Index, Division.
  - *Health Metrics*: BMI, Diabetes Status.
  - *Environmental/Household*: Crowding Index, Household Size, Electricity, Water Source, Sanitation.

## 🛠️ Methodology

### 1. Statistical & Epidemiological Analysis (R)
- **Complex Survey Design**: Utilized the `survey` package to account for BDHS sampling weights (`sample_weight`), stratification (`strata_id`), and clustering (`cluster_id`) to ensure nationally representative estimates.
- **Geospatial Visualization**: Mapped hypertension prevalence across 8 administrative divisions using `geom_sf` and `viridis` color scales.
- **Multivariate Analysis**: Fitted weighted logistic regression models to calculate Adjusted Odds Ratios (AOR) with 95% Confidence Intervals, visualized via publication-ready Forest Plots.

### 2. Machine Learning Pipeline (Python)
- **Survey-Aware Data Splitting**: Used `GroupShuffleSplit` based on `cluster_id` to prevent data leakage between train, validation, and test sets, respecting the complex survey design.
- **Models Implemented**: XGBoost, Random Forest, Gradient Boosting, Extra Trees, Logistic Regression, SVC, and Decision Trees.
- **Class Imbalance Handling**: Applied `scale_pos_weight` (XGBoost) and `class_weight='balanced'` to address the skewed distribution of hypertension cases.
- **Probability Calibration**: Applied Platt Scaling (Logistic Regression on logits) to ensure ML predicted probabilities reflect true epidemiological risks.
- **Evaluation Metrics**: 
  - *Performance*: ROC-AUC, PR-AUC, F1-Score, Sensitivity, Specificity, Cohen's Kappa.
  - *Robustness*: Bootstrap resampling (1000 iterations) for 95% Confidence Intervals.
  - *Calibration*: Brier Score and Calibration Slope/Intercept.
- **Interpretability**: SHAP (SHapley Additive exPlanations) values for global and local feature importance.

## ⚙️ Requirements

### R Dependencies
```R
install.packages(c("survey", "car", "dplyr", "broom", "ggplot2", "forcats", "stringr", "readr", "scales", "sf", "patchwork"))