"""
Assessing the Prevalence and Determinants of Hypertension in Bangladesh
Core package for survey analysis and machine learning pipeline.
"""

__version__ = "1.0.0"
__author__ = "Your Name"
__email__ = "your.email@juniv.edu"

# Import core configurations
from src.config import (
    RAW_DATA_PATH, PROCESSED_DATA_PATH, OUTPUT_DIR,
    TARGET, DESIGN_COLS, CATEGORICAL_VARS, RANDOM_SEED,
    FIGURE_DPI, FIGURE_WIDTH, FIGURE_HEIGHT
)

# Import core functions for easy access
from src.data_processing import load_and_preprocess_data
from src.survey_stats import calculate_weighted_prevalence, run_survey_logistic_regression
from src.ml_pipeline import cluster_aware_split, get_models, calibrate_and_tune
from src.evaluation import evaluate_metrics, bootstrap_ci
from src.visualization import plot_roc_pr_curves, plot_model_performance