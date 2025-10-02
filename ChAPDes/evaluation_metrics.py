"""
Evaluation Metrics for Imbalanced Classification

This module provides comprehensive evaluation metrics specifically designed
for assessing model performance on imbalanced binary classification tasks.
"""

import numpy as np
import pandas as pd
from sklearn.metrics import (
    confusion_matrix, classification_report, roc_auc_score,
    average_precision_score, f1_score, precision_recall_curve,
    roc_curve, auc, matthews_corrcoef
)
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple, Optional, Union
import warnings
warnings.filterwarnings('ignore')


class ImbalancedClassificationMetrics:
    """
    Comprehensive evaluation metrics for imbalanced binary classification.

    This class provides both standard and specialized metrics for evaluating
    model performance on imbalanced datasets, including metrics that are
    more informative than accuracy for imbalanced problems.
    """

    def __init__(self):
        """Initialize the metrics calculator."""
        pass

    @staticmethod
    def calculate_all_metrics(y_true: np.ndarray, y_pred: np.ndarray,
                            y_pred_proba: Optional[np.ndarray] = None) -> Dict:
        """
        Calculate comprehensive evaluation metrics for binary classification.

        Args:
            y_true: True labels
            y_pred: Predicted labels
            y_pred_proba: Predicted probabilities (optional, for probability-based metrics)

        Returns:
            Dictionary containing all calculated metrics
        """
        # Confusion matrix
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()

        # Basic metrics
        accuracy = (tp + tn) / (tp + tn + fp + fn)
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0

        # F-scores
        f1 = f1_score(y_true, y_pred) if (precision + recall) > 0 else 0
        f2 = (5 * precision * recall) / (4 * precision + recall) if (4 * precision + recall) > 0 else 0
        f0_5 = (1.25 * precision * recall) / (0.25 * precision + recall) if (0.25 * precision + recall) > 0 else 0

        # Additional metrics
        npv = tn / (tn + fn) if (tn + fn) > 0 else 0  # Negative Predictive Value
        false_positive_rate = fp / (fp + tn) if (fp + tn) > 0 else 0
        false_negative_rate = fn / (fn + tp) if (fn + tp) > 0 else 0

        # Matthews Correlation Coefficient
        mcc = matthews_corrcoef(y_true, y_pred)

        # Balanced Accuracy
        balanced_accuracy = (recall + specificity) / 2

        metrics = {
            # Basic counts
            'tp': int(tp),
            'tn': int(tn),
            'fp': int(fp),
            'fn': int(fn),

            # Standard metrics
            'accuracy': float(accuracy),
            'precision': float(precision),
            'recall': float(recall),
            'specificity': float(specificity),
            'npv': float(npv),

            # F-scores
            'f1_score': float(f1),
            'f2_score': float(f2),
            'f05_score': float(f0_5),

            # Rates
            'false_positive_rate': float(false_positive_rate),
            'false_negative_rate': float(false_negative_rate),

            # Advanced metrics
            'matthews_corrcoef': float(mcc),
            'balanced_accuracy': float(balanced_accuracy)
        }

        # Probability-based metrics (if probabilities provided)
        if y_pred_proba is not None:
            prob_metrics = ImbalancedClassificationMetrics._calculate_probability_metrics(
                y_true, y_pred_proba
            )
            metrics.update(prob_metrics)

        return metrics

    @staticmethod
    def _calculate_probability_metrics(y_true: np.ndarray, y_pred_proba: np.ndarray) -> Dict:
        """Calculate probability-based metrics."""
        # ROC AUC
        try:
            roc_auc = roc_auc_score(y_true, y_pred_proba)
        except:
            roc_auc = 0.0

        # Precision-Recall AUC
        try:
            pr_auc = average_precision_score(y_true, y_pred_proba)
        except:
            pr_auc = 0.0

        # Brier Score
        brier_score = np.mean((y_true - y_pred_proba) ** 2)

        # Log Loss (cross-entropy)
        try:
            log_loss = -np.mean(y_true * np.log(y_pred_proba + 1e-15) +
                              (1 - y_true) * np.log(1 - y_pred_proba + 1e-15))
        except:
            log_loss = float('inf')

        return {
            'roc_auc': float(roc_auc),
            'pr_auc': float(pr_auc),
            'brier_score': float(brier_score),
            'log_loss': float(log_loss)
        }

    @staticmethod
    def print_metrics_summary(metrics: Dict, model_name: str = "Model"):
        """
        Print a formatted summary of evaluation metrics.

        Args:
            metrics: Dictionary of metrics from calculate_all_metrics
            model_name: Name of the model for display
        """
        print(f"\n{'='*60}")
        print(f"{model_name.upper()} - EVALUATION METRICS")
        print('='*60)

        # Confusion matrix section
        print("CONFUSION MATRIX:")
        print(f"True Positives:  {metrics['tp']}")
        print(f"True Negatives:  {metrics['tn']}")
        print(f"False Positives: {metrics['fp']}")
        print(f"False Negatives: {metrics['fn']}")
        print()

        # Standard metrics section
        print("STANDARD METRICS:")
        print(f"Accuracy:        {metrics['accuracy']".4f"}")
        print(f"Precision:       {metrics['precision']".4f"}")
        print(f"Recall:          {metrics['recall']".4f"}")
        print(f"Specificity:     {metrics['specificity']".4f"}")
        print(f"NPV:             {metrics['npv']".4f"}")
        print()

        # F-scores section
        print("F-SCORES:")
        print(f"F1-Score:        {metrics['f1_score']".4f"}")
        print(f"F2-Score:        {metrics['f2_score']".4f"}")
        print(f"F0.5-Score:      {metrics['f05_score']".4f"}")
        print()

        # Advanced metrics section
        print("ADVANCED METRICS:")
        print(f"Balanced Acc:    {metrics['balanced_accuracy']".4f"}")
        print(f"MCC:             {metrics['matthews_corrcoef']".4f"}")

        if 'roc_auc' in metrics:
            print()
            print("PROBABILITY-BASED METRICS:")
            print(f"ROC AUC:         {metrics['roc_auc']".4f"}")
            print(f"PR AUC:          {metrics['pr_auc']".4f"}")
            print(f"Brier Score:     {metrics['brier_score']".4f"}")
            print(f"Log Loss:        {metrics['log_loss']".4f"}")

    @staticmethod
    def plot_roc_curve(y_true: np.ndarray, y_pred_proba: np.ndarray,
                      title: str = "ROC Curve"):
        """
        Plot ROC curve.

        Args:
            y_true: True labels
            y_pred_proba: Predicted probabilities
            title: Plot title
        """
        fpr, tpr, _ = roc_curve(y_true, y_pred_proba)
        roc_auc = auc(fpr, tpr)

        plt.figure(figsize=(8, 6))
        plt.plot(fpr, tpr, color='darkorange', lw=2,
                label=f'ROC curve (AUC = {roc_auc".3f"})')
        plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(title)
        plt.legend(loc="lower right")
        plt.grid(True, alpha=0.3)
        plt.show()

    @staticmethod
    def plot_precision_recall_curve(y_true: np.ndarray, y_pred_proba: np.ndarray,
                                  title: str = "Precision-Recall Curve"):
        """
        Plot Precision-Recall curve.

        Args:
            y_true: True labels
            y_pred_proba: Predicted probabilities
            title: Plot title
        """
        precision, recall, _ = precision_recall_curve(y_true, y_pred_proba)
        pr_auc = auc(recall, precision)

        plt.figure(figsize=(8, 6))
        plt.plot(recall, precision, color='green', lw=2,
                label=f'PR curve (AUC = {pr_auc".3f"})')
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.title(title)
        plt.legend(loc="lower right")
        plt.grid(True, alpha=0.3)

        # Add baseline (random classifier)
        baseline = np.sum(y_true) / len(y_true)
        plt.axhline(y=baseline, color='red', linestyle='--', alpha=0.7,
                   label=f'Baseline ({baseline".3f"})')

        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.show()

    @staticmethod
    def plot_confusion_matrix(y_true: np.ndarray, y_pred: np.ndarray,
                            title: str = "Confusion Matrix", normalize: bool = False):
        """
        Plot confusion matrix as a heatmap.

        Args:
            y_true: True labels
            y_pred: Predicted labels
            title: Plot title
            normalize: Whether to normalize by true labels
        """
        cm = confusion_matrix(y_true, y_pred)

        if normalize:
            cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
            fmt = '.3f'
            title += " (Normalized)"
        else:
            fmt = 'd'

        plt.figure(figsize=(8, 6))
        sns.heatmap(cm, annot=True, fmt=fmt, cmap='Blues',
                   xticklabels=['Predicted Negative', 'Predicted Positive'],
                   yticklabels=['Actual Negative', 'Actual Positive'])
        plt.title(title)
        plt.ylabel('True Label')
        plt.xlabel('Predicted Label')
        plt.tight_layout()
        plt.show()

    @staticmethod
    def plot_metrics_comparison(metrics_dict: Dict[str, Dict],
                              metrics_to_plot: List[str] = None,
                              title: str = "Model Comparison"):
        """
        Compare multiple models' performance metrics.

        Args:
            metrics_dict: Dictionary with model names as keys and metrics as values
            metrics_to_plot: List of metric names to include in comparison
            title: Plot title
        """
        if metrics_to_plot is None:
            metrics_to_plot = ['f1_score', 'precision', 'recall', 'balanced_accuracy', 'roc_auc']

        # Filter metrics that exist in all models
        available_metrics = set(metrics_dict[next(iter(metrics_dict))].keys())
        for model_metrics in metrics_dict.values():
            available_metrics = available_metrics.intersection(set(model_metrics.keys()))

        metrics_to_plot = [m for m in metrics_to_plot if m in available_metrics]

        if not metrics_to_plot:
            print("No common metrics found for comparison")
            return

        # Create comparison dataframe
        comparison_data = []
        for model_name, metrics in metrics_dict.items():
            row = {'Model': model_name}
            for metric in metrics_to_plot:
                if metric in metrics:
                    row[metric] = metrics[metric]
            comparison_data.append(row)

        df = pd.DataFrame(comparison_data)

        # Plot
        fig, ax = plt.subplots(figsize=(12, 6))

        # Create bar plot
        x = np.arange(len(metrics_to_plot))
        width = 0.8 / len(df)

        for i, (_, row) in enumerate(df.iterrows()):
            ax.bar(x + i * width, [row[metric] for metric in metrics_to_plot],
                  width, label=row['Model'], alpha=0.8)

        ax.set_xlabel('Metrics')
        ax.set_ylabel('Score')
        ax.set_title(title)
        ax.set_xticks(x + width * (len(df) - 1) / 2)
        ax.set_xticklabels(metrics_to_plot, rotation=45, ha='right')
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

        return df

    @staticmethod
    def find_optimal_threshold(y_true: np.ndarray, y_pred_proba: np.ndarray,
                             metric: str = 'f1') -> Tuple[float, Dict]:
        """
        Find optimal probability threshold for different metrics.

        Args:
            y_true: True labels
            y_pred_proba: Predicted probabilities
            metric: Metric to optimize ('f1', 'precision', 'recall', 'balanced_accuracy')

        Returns:
            Tuple of (optimal_threshold, metrics_at_threshold)
        """
        thresholds = np.linspace(0.1, 0.9, 81)  # Test thresholds from 0.1 to 0.9

        best_score = 0
        best_threshold = 0.5
        best_metrics = {}

        for threshold in thresholds:
            y_pred = (y_pred_proba >= threshold).astype(int)

            if metric == 'f1':
                score = f1_score(y_true, y_pred)
            elif metric == 'precision':
                tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
                score = tp / (tp + fp) if (tp + fp) > 0 else 0
            elif metric == 'recall':
                tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
                score = tp / (tp + fn) if (tp + fn) > 0 else 0
            elif metric == 'balanced_accuracy':
                tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
                recall = tp / (tp + fn) if (tp + fn) > 0 else 0
                specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
                score = (recall + specificity) / 2
            else:
                raise ValueError(f"Unknown metric: {metric}")

            if score > best_score:
                best_score = score
                best_threshold = threshold
                tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
                best_metrics = {
                    'threshold': threshold,
                    'f1_score': f1_score(y_true, y_pred),
                    'precision': tp / (tp + fp) if (tp + fp) > 0 else 0,
                    'recall': tp / (tp + fn) if (tp + fn) > 0 else 0,
                    'tp': int(tp), 'tn': int(tn), 'fp': int(fp), 'fn': int(fn)
                }

        return best_threshold, best_metrics


def evaluate_model_at_thresholds(y_true: np.ndarray, y_pred_proba: np.ndarray,
                               thresholds: List[float] = None) -> pd.DataFrame:
    """
    Evaluate model performance at different probability thresholds.

    Args:
        y_true: True labels
        y_pred_proba: Predicted probabilities
        thresholds: List of thresholds to evaluate (if None, use common thresholds)

    Returns:
        DataFrame with metrics at each threshold
    """
    if thresholds is None:
        thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    results = []

    for threshold in thresholds:
        y_pred = (y_pred_proba >= threshold).astype(int)

        tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()

        result = {
            'threshold': threshold,
            'tp': int(tp),
            'tn': int(tn),
            'fp': int(fp),
            'fn': int(fn),
            'precision': tp / (tp + fp) if (tp + fp) > 0 else 0,
            'recall': tp / (tp + fn) if (tp + fn) > 0 else 0,
            'f1_score': f1_score(y_true, y_pred),
            'specificity': tn / (tn + fp) if (tn + fp) > 0 else 0
        }
        results.append(result)

    return pd.DataFrame(results)


def bootstrap_metrics(y_true: np.ndarray, y_pred: np.ndarray,
                     y_pred_proba: Optional[np.ndarray] = None,
                     n_bootstrap: int = 1000, confidence_level: float = 0.95) -> Dict:
    """
    Calculate confidence intervals for metrics using bootstrap resampling.

    Args:
        y_true: True labels
        y_pred: Predicted labels
        y_pred_proba: Predicted probabilities
        n_bootstrap: Number of bootstrap samples
        confidence_level: Confidence level for intervals

    Returns:
        Dictionary with metrics and their confidence intervals
    """
    n_samples = len(y_true)
    alpha = 1 - confidence_level

    # Bootstrap samples
    bootstrap_metrics = []

    for _ in range(n_bootstrap):
        # Sample with replacement
        indices = np.random.choice(n_samples, size=n_samples, replace=True)
        y_true_boot = y_true[indices]
        y_pred_boot = y_pred[indices]

        if y_pred_proba is not None:
            y_pred_proba_boot = y_pred_proba[indices]
        else:
            y_pred_proba_boot = None

        # Calculate metrics for this bootstrap sample
        metrics = ImbalancedClassificationMetrics.calculate_all_metrics(
            y_true_boot, y_pred_boot, y_pred_proba_boot
        )
        bootstrap_metrics.append(metrics)

    # Convert to DataFrame for easier calculation of confidence intervals
    metrics_df = pd.DataFrame(bootstrap_metrics)

    # Calculate confidence intervals
    confidence_intervals = {}

    for metric in metrics_df.columns:
        if metric in ['tp', 'tn', 'fp', 'fn']:
            # For counts, use percentile method
            lower_pct = 100 * alpha / 2
            upper_pct = 100 * (1 - alpha / 2)
            ci_lower = np.percentile(metrics_df[metric], lower_pct)
            ci_upper = np.percentile(metrics_df[metric], upper_pct)
        else:
            # For rates and scores, use percentile method
            lower_pct = 100 * alpha / 2
            upper_pct = 100 * (1 - alpha / 2)
            ci_lower = np.percentile(metrics_df[metric], lower_pct)
            ci_upper = np.percentile(metrics_df[metric], upper_pct)

        confidence_intervals[metric] = {
            'mean': metrics_df[metric].mean(),
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'confidence_level': confidence_level
        }

    return confidence_intervals


# Example usage function
def example_metrics_evaluation():
    """Demonstrate evaluation metrics functionality."""

    # Create sample data
    np.random.seed(42)
    n_samples = 1000
    n_positive = 100  # Imbalanced case

    # Generate synthetic data
    y_true = np.zeros(n_samples)
    y_true[:n_positive] = 1
    np.random.shuffle(y_true)

    # Generate predictions (simulated model)
    y_pred = np.random.choice([0, 1], size=n_samples, p=[0.7, 0.3])
    y_pred_proba = np.random.random(n_samples)

    print("Example Metrics Evaluation")
    print(f"Dataset: {n_samples} samples, {n_positive} positive cases ({n_positive/n_samples*100".1f"}%)")

    # Calculate all metrics
    metrics = ImbalancedClassificationMetrics.calculate_all_metrics(y_true, y_pred, y_pred_proba)

    # Print summary
    ImbalancedClassificationMetrics.print_metrics_summary(metrics, "Example Model")

    # Plot curves
    print("\nGenerating plots...")
    ImbalancedClassificationMetrics.plot_roc_curve(y_true, y_pred_proba)
    ImbalancedClassificationMetrics.plot_precision_recall_curve(y_true, y_pred_proba)
    ImbalancedClassificationMetrics.plot_confusion_matrix(y_true, y_pred)

    # Find optimal threshold
    optimal_threshold, threshold_metrics = ImbalancedClassificationMetrics.find_optimal_threshold(
        y_true, y_pred_proba, metric='f1'
    )

    print(f"\nOptimal threshold for F1-score: {optimal_threshold".3f"}")
    print("Metrics at optimal threshold:")
    for key, value in threshold_metrics.items():
        print(f"  {key}: {value}")

    return metrics


if __name__ == "__main__":
    example_metrics_evaluation()