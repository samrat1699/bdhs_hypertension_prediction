import pandas as pd
import numpy as np
from src.config import RAW_DATA_PATH, PROCESSED_DATA_PATH, TARGET, CATEGORICAL_VARS

def load_and_preprocess_data() -> pd.DataFrame:
    """Loads raw data, creates bins, and formats categorical variables."""
    print("🔄 Loading and preprocessing data...")
    df = pd.read_csv(RAW_DATA_PATH)
    
    # 1. Create Household Size Categories
    df['household_size_cat'] = pd.cut(
        df['household_size'],
        bins=[0, 3, 6, np.inf],
        labels=['Small (1-3)', 'Medium (4-6)', 'Large (7+)']
    )
    
    # 2. Create Crowding Categories
    df['crowding_category'] = pd.cut(
        df['crowding_index'],
        bins=[0, 2, 3, np.inf],
        labels=['Low (<2)', 'Moderate (2-3)', 'High (>3)'],
        include_lowest=True
    )
    
    # 3. Ensure Target is Integer
    df[TARGET] = df[TARGET].astype(int)
    
    # 4. Format Categorical Variables as Strings for One-Hot Encoding
    for col in CATEGORICAL_VARS:
        if col in df.columns:
            df[col] = df[col].astype(str)
            
    df.to_csv(PROCESSED_DATA_PATH, index=False)
    print(f"✅ Data saved to {PROCESSED_DATA_PATH} | Shape: {df.shape}")
    return df

if __name__ == "__main__":
    load_and_preprocess_data()