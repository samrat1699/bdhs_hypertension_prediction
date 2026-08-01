import numpy as np
from sklearn.metrics import (roc_auc_score, accuracy_score, cohen_kappa_score, 
                             precision_score, recall_score, f1_score, brier_score_loss, confusion_matrix)
from sklearn.linear_model import LogisticRegression

def calibration_slope(y_true, p, w):
    """Calculates calibration slope using weighted logistic regression."""
    p = np.clip(p, 1e-6, 1-1e-6)
    logit = np.log(p / (1 - p)).reshape(-1, 1)
    lr = LogisticRegression(fit_intercept=False).fit(logit, y_true, sample_weight=w)
    return lr.coef_[0][0]

def evaluate_metrics(y_true, p, w, t):
    """Calculates all metrics including Sensitivity/Specificity."""
    pred = (p >= t).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, pred, sample_weight=w).ravel()
    
    return {
        "Acc": accuracy_score(y_true, pred, sample_weight=w),
        "AUC": roc_auc_score(y_true, p, sample_weight=w),
        "Kappa": cohen_kappa_score(y_true, pred, sample_weight=w),
        "Sens": tp / (tp + fn + 1e-9),
        "Spec": tn / (tn + fp + 1e-9),
        "Brier": brier_score_loss(y_true, p, sample_weight=w),
        "Slope": calibration_slope(y_true, p, w)
    }

def bootstrap_ci(y, p, w, t, B=1000):
    """Calculates 95% Confidence Intervals via Bootstrap."""
    rng = np.random.RandomState(42)
    n = len(y)
    boot_data = []
    
    for _ in range(B):
        idx = rng.choice(np.arange(n), size=n, replace=True)
        y_b, p_b, w_b = y.iloc[idx], p[idx], w.iloc[idx]
        if len(np.unique(y_b)) < 2: continue
        boot_data.append(list(evaluate_metrics(y_b, p_b, w_b, t).values()))
        
    boot_data = np.array(boot_data)
    means = np.mean(boot_data, axis=0)
    lower = np.percentile(boot_data, 2.5, axis=0)
    upper = np.percentile(boot_data, 97.5, axis=0)
    
    return [f"{m:.3f} ({l:.3f}-{u:.3f})" for m, l, u in zip(means, lower, upper)]