import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.stats.weightstats import DescrStatsW
from src.config import PROCESSED_DATA_PATH, TARGET

def calculate_weighted_prevalence(df: pd.DataFrame) -> None:
    """Calculates overall and age-stratified weighted prevalence."""
    print("\n📊 Calculating Weighted Prevalence...")
    y = df[TARGET].values
    w = df['sample_weight'].values
    
    p_hat = np.sum(w * y) / np.sum(w)
    n_eff = (np.sum(w) ** 2) / np.sum(w ** 2)
    weighted_var = np.sum(w * (y - p_hat) ** 2) / ((np.sum(w) - 1) * np.sum(w) / np.sum(w))
    se = np.sqrt(weighted_var / n_eff)
    
    print(f"Weighted Prevalence: {p_hat*100:.2f}%")
    print(f"95% CI: {(p_hat - 1.96*se)*100:.2f}% – {(p_hat + 1.96*se)*100:.2f}%")
    print(f"Effective Sample Size: {n_eff:,.0f}")

def run_survey_logistic_regression(df: pd.DataFrame) -> pd.DataFrame:
    """Runs survey-weighted logistic regression to get Adjusted Odds Ratios (AOR)."""
    print("\n Running Survey-Weighted Logistic Regression...")
    
    # Prepare data for statsmodels (using sample weights as frequency weights for approximation)
    X = pd.get_dummies(df.drop(columns=[TARGET, 'cluster_id', 'strata_id']), drop_first=True).astype(float)
    y = df[TARGET]
    w = df['sample_weight']
    
    # Using statsmodels GLM with Poisson/Binomial family and weights
    # Note: For complex survey designs, Python's statsmodels is an approximation of R's 'survey' package.
    model = sm.GLM(y, sm.add_constant(X), family=sm.families.Binomial(), freq_weights=w)
    results = model.fit()
    
    summary_df = pd.DataFrame({
        'Variable': X.columns,
        'AOR': np.exp(results.params[1:]),
        'Lower_CI': np.exp(results.conf_int()[1:, 0]),
        'Upper_CI': np.exp(results.conf_int()[1:, 1]),
        'P_Value': results.pvalues[1:]
    }).sort_values(by='P_Value')
    
    summary_df.to_csv("outputs/tables/logistic_regression_results.csv", index=False)
    print("✅ Logistic Regression results saved.")
    return summary_df