import numpy as np
from sklearn.model_selection import GroupShuffleSplit
from sklearn.ensemble import (RandomForestClassifier, GradientBoostingClassifier, ExtraTreesClassifier)
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from xgboost import XGBClassifier
from src.config import RANDOM_SEED, TEST_SIZE, VAL_SIZE

def cluster_aware_split(X, y, w, clusters):
    """Splits data respecting cluster_id to prevent leakage."""
    # 80% Train+Val / 20% Test
    gss1 = GroupShuffleSplit(n_splits=1, test_size=TEST_SIZE, random_state=RANDOM_SEED)
    train_val_idx, test_idx = next(gss1.split(X, y, groups=clusters))
    
    # 75% Train / 25% Val
    gss2 = GroupShuffleSplit(n_splits=1, test_size=VAL_SIZE, random_state=RANDOM_SEED)
    train_idx, val_idx = next(gss2.split(
        X.iloc[train_val_idx], y.iloc[train_val_idx], groups=clusters.iloc[train_val_idx]
    ))
    
    return (
        X.iloc[train_idx], y.iloc[train_idx], w.iloc[train_idx],
        X.iloc[val_idx], y.iloc[val_idx], w.iloc[val_idx],
        X.iloc[test_idx], y.iloc[test_idx], w.iloc[test_idx]
    )

def get_models(class_ratio):
    """Returns dictionary of models with specific hyperparameters from your notebook."""
    return {
        "XGBoost": XGBClassifier(
            n_estimators=1000, max_depth=4, learning_rate=0.03, 
            subsample=0.8, colsample_bytree=0.8, scale_pos_weight=class_ratio,
            eval_metric="auc", tree_method="hist", random_state=RANDOM_SEED
        ),
        "Random Forest": RandomForestClassifier(
            n_estimators=800, max_depth=None, min_samples_leaf=10,
            class_weight="balanced", random_state=RANDOM_SEED
        ),
        # ... (Add Gradient Boosting, Extra Trees, LR, SVC, DT as per notebook)
    }