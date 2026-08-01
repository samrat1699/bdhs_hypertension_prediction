import pandas as pd
from src.data_preprocessing import load_and_preprocess_data
from src.survey_analysis import calculate_weighted_prevalence, run_survey_logistic_regression
from src.ml_pipeline import cluster_aware_split, get_models, calibrate_and_tune
from src.evaluation import evaluate_metrics, bootstrap_ci
from src.visualization import plot_roc_pr_curves, plot_model_performance
from src.config import PROCESSED_DATA_PATH

import sys
import os

# Add the root directory to the system path to ensure imports work
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# Import directly from the src package (thanks to __init__.py)
from src import (
    load_and_preprocess_data,
    calculate_weighted_prevalence,
    cluster_aware_split,
    get_models,
    calibrate_and_tune,
    evaluate_metrics,
    bootstrap_ci,
    plot_roc_pr_curves,
    plot_model_performance,
    TARGET,
    FIGURE_DPI
)

def main():
    print("="*60)
    print("HYPERTENSION IN BANGLADESH: ML & SURVEY ANALYSIS PIPELINE")
    print("="*60)
    
    # 1. Data Preprocessing
    df = load_and_preprocess_data()
    
    # 2. Survey Analysis
    calculate_weighted_prevalence(df)
    
    # 3. ML Pipeline Setup
    X_train, y_train, w_train, X_val, y_val, w_val, X_test, y_test, w_test = cluster_aware_split(df)
    class_ratio = (y_train == 0).sum() / (y_train == 1).sum()
    models = get_models(class_ratio)
    
    master_results = []
    plot_data = {}
    
    print("\n🤖 Training and Evaluating Machine Learning Models...")
    for name, model in models.items():
        print(f"  -> Training {name}...")
        train_prob, val_prob, best_t = calibrate_and_tune(model, X_train, y_train, w_train, X_val, y_val, w_val)
        
        # Get Test Probabilities
        test_prob = model.predict_proba(X_test)[:, 1]
        
        # Evaluate & Bootstrap
        boot_results = bootstrap_ci(y_test, test_prob, w_test, best_t)
        master_results.append([f"{name} (Test)"] + boot_results)
        plot_data[name] = test_prob
        
    # 4. Visualizations
    plot_roc_pr_curves(y_test, w_test, plot_data, cv_summary=None)
    
    print(f"\n✅ Pipeline completed successfully. Check the 'outputs/' directory.")

if __name__ == "__main__":
    main()