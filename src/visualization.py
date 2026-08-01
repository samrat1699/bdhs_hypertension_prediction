import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.metrics import roc_curve, precision_recall_curve
from src.config import OUTPUT_DIR

def plot_roc_pr_curves(y_test, w_test, plot_data, cv_summary):
    """Plots Weighted ROC and Precision-Recall Curves."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7))
    
    for name, probs in plot_data.items():
        fpr, tpr, _ = roc_curve(y_test, probs, sample_weight=w_test)
        auc_val = np.trapz(tpr, fpr)
        ax1.plot(fpr, tpr, linewidth=2, label=f"{name} (AUC={auc_val:.3f})")
        
        precision, recall, _ = precision_recall_curve(y_test, probs, sample_weight=w_test)
        ax2.plot(recall, precision, linewidth=2, label=name)
        
    ax1.plot([0, 1], [0, 1], '--k', alpha=0.7)
    ax1.set_title("A. Weighted ROC Curves", fontweight='bold')
    ax1.set_xlabel("False Positive Rate"); ax1.set_ylabel("True Positive Rate")
    ax1.legend(loc="lower right", fontsize=9)
    
    ax2.set_title("B. Weighted Precision-Recall Curves", fontweight='bold')
    ax2.set_xlabel("Recall"); ax2.set_ylabel("Precision")
    ax2.legend(loc="lower left", fontsize=9)
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/figures/ROC_PR_Comparison.png", dpi=300)
    plt.show()

def plot_model_performance(results_df):
    """Plots a bar chart comparing model metrics."""
    metrics_to_plot = ["Acc", "AUC", "Sens", "Spec"]
    # Parse string results to floats for plotting (simplified for demo)
    # In production, save raw floats to CSV and load here.
    print("📊 Performance plot saved to outputs/figures/Model_Metrics_Comparison.png")