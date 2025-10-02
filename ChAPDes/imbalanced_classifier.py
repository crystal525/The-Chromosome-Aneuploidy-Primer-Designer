"""
Imbalanced Binary Classifier for Genomic Data

This module provides tools for handling imbalanced binary classification problems
commonly encountered in genomic and biomedical data analysis.
"""

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, StratifiedKFold, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    classification_report, confusion_matrix, roc_auc_score,
    precision_recall_curve, average_precision_score, f1_score,
    roc_curve, auc
)
from sklearn.utils.class_weight import compute_class_weight
from imblearn.over_sampling import SMOTE, ADASYN
from imblearn.under_sampling import RandomUnderSampler
from imblearn.ensemble import BalancedRandomForestClassifier, EasyEnsembleClassifier
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple, Optional, Union
import warnings
warnings.filterwarnings('ignore')


class ImbalancedBinaryClassifier:
    """
    A comprehensive class for handling imbalanced binary classification problems.

    This class provides multiple strategies for dealing with class imbalance:
    - Class weighting
    - Oversampling (SMOTE, ADASYN)
    - Undersampling
    - Ensemble methods (Balanced Random Forest, Easy Ensemble)
    """

    def __init__(self, method: str = 'class_weight', random_state: int = 42):
        """
        Initialize the imbalanced classifier.

        Args:
            method: Strategy for handling imbalance
                - 'class_weight': Use class weights
                - 'smote': Use SMOTE oversampling
                - 'adasyn': Use ADASYN oversampling
                - 'undersample': Use random undersampling
                - 'balanced_rf': Use Balanced Random Forest
                - 'easy_ensemble': Use Easy Ensemble
            random_state: Random seed for reproducibility
        """
        self.method = method
        self.random_state = random_state
        self.model = None
        self.scaler = StandardScaler()
        self.feature_names = None
        self.is_fitted = False

        # Initialize model based on method
        self._initialize_model()

    def _initialize_model(self):
        """Initialize the appropriate model based on the chosen method."""
        if self.method == 'class_weight':
            self.model = LogisticRegression(random_state=self.random_state, max_iter=1000)
        elif self.method == 'smote':
            self.model = LogisticRegression(random_state=self.random_state, max_iter=1000)
        elif self.method == 'adasyn':
            self.model = LogisticRegression(random_state=self.random_state, max_iter=1000)
        elif self.method == 'undersample':
            self.model = LogisticRegression(random_state=self.random_state, max_iter=1000)
        elif self.method == 'balanced_rf':
            self.model = BalancedRandomForestClassifier(
                random_state=self.random_state, n_estimators=100
            )
        elif self.method == 'easy_ensemble':
            self.model = EasyEnsembleClassifier(
                random_state=self.random_state, n_estimators=10
            )
        else:
            raise ValueError(f"Unknown method: {self.method}")

    def fit(self, X: Union[pd.DataFrame, np.ndarray],
            y: Union[pd.Series, np.ndarray],
            validation_split: float = 0.2) -> 'ImbalancedBinaryClassifier':
        """
        Fit the classifier to the training data.

        Args:
            X: Feature matrix
            y: Target vector (binary: 0 or 1)
            validation_split: Fraction of data to use for validation

        Returns:
            Self for method chaining
        """
        # Convert to numpy arrays if needed
        if isinstance(X, pd.DataFrame):
            self.feature_names = X.columns.tolist()
            X = X.values
        if isinstance(y, pd.Series):
            y = y.values

        # Split the data
        X_train, X_val, y_train, y_val = train_test_split(
            X, y, test_size=validation_split, random_state=self.random_state,
            stratify=y
        )

        # Scale features
        X_train_scaled = self.scaler.fit_transform(X_train)
        X_val_scaled = self.scaler.transform(X_val)

        # Apply balancing technique if needed
        if self.method in ['smote', 'adasyn']:
            X_train_balanced, y_train_balanced = self._apply_oversampling(
                X_train_scaled, y_train
            )
        elif self.method == 'undersample':
            X_train_balanced, y_train_balanced = self._apply_undersampling(
                X_train_scaled, y_train
            )
        else:
            X_train_balanced, y_train_balanced = X_train_scaled, y_train

        # Fit the model
        if self.method == 'class_weight':
            # Calculate class weights
            class_weights = compute_class_weight('balanced', classes=np.unique(y_train), y=y_train)
            class_weight_dict = dict(zip(np.unique(y_train), class_weights))
            self.model.set_params(class_weight=class_weight_dict)

        self.model.fit(X_train_balanced, y_train_balanced)
        self.is_fitted = True

        # Store validation data for evaluation
        self.X_val = X_val_scaled
        self.y_val = y_val

        return self

    def _apply_oversampling(self, X: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Apply oversampling technique."""
        if self.method == 'smote':
            sampler = SMOTE(random_state=self.random_state)
        elif self.method == 'adasyn':
            sampler = ADASYN(random_state=self.random_state)
        else:
            return X, y

        return sampler.fit_resample(X, y)

    def _apply_undersampling(self, X: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Apply undersampling technique."""
        sampler = RandomUnderSampler(random_state=self.random_state)
        return sampler.fit_resample(X, y)

    def predict(self, X: Union[pd.DataFrame, np.ndarray]) -> np.ndarray:
        """Make predictions on new data."""
        if not self.is_fitted:
            raise ValueError("Model must be fitted before making predictions")

        if isinstance(X, pd.DataFrame):
            X = X.values

        X_scaled = self.scaler.transform(X)
        return self.model.predict(X_scaled)

    def predict_proba(self, X: Union[pd.DataFrame, np.ndarray]) -> np.ndarray:
        """Predict class probabilities."""
        if not self.is_fitted:
            raise ValueError("Model must be fitted before making predictions")

        if isinstance(X, pd.DataFrame):
            X = X.values

        X_scaled = self.scaler.transform(X)
        return self.model.predict_proba(X_scaled)

    def evaluate(self, X: Union[pd.DataFrame, np.ndarray] = None,
                 y: Union[pd.Series, np.ndarray] = None,
                 plot: bool = True) -> Dict:
        """
        Evaluate the model performance.

        Args:
            X: Test features (if None, uses validation data)
            y: Test targets (if None, uses validation data)
            plot: Whether to create evaluation plots

        Returns:
            Dictionary containing evaluation metrics
        """
        if not self.is_fitted:
            raise ValueError("Model must be fitted before evaluation")

        if X is None or y is None:
            X, y = self.X_val, self.y_val

        if isinstance(X, pd.DataFrame):
            X = X.values
        if isinstance(y, pd.Series):
            y = y.values

        X_scaled = self.scaler.transform(X)
        y_pred = self.predict(X_scaled)
        y_pred_proba = self.predict_proba(X_scaled)[:, 1]

        # Calculate metrics
        metrics = self._calculate_metrics(y, y_pred, y_pred_proba)

        if plot:
            self._plot_evaluation(y, y_pred, y_pred_proba, metrics)

        return metrics

    def _calculate_metrics(self, y_true: np.ndarray, y_pred: np.ndarray,
                          y_pred_proba: np.ndarray) -> Dict:
        """Calculate comprehensive evaluation metrics."""
        # Basic metrics
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()

        metrics = {
            'accuracy': (tp + tn) / (tp + tn + fp + fn),
            'precision': tp / (tp + fp) if (tp + fp) > 0 else 0,
            'recall': tp / (tp + fn) if (tp + fn) > 0 else 0,
            'f1_score': f1_score(y_true, y_pred),
            'specificity': tn / (tn + fp) if (tn + fp) > 0 else 0,
            'npv': tn / (tn + fn) if (tn + fn) > 0 else 0,  # Negative Predictive Value
            'roc_auc': roc_auc_score(y_true, y_pred_proba),
            'pr_auc': average_precision_score(y_true, y_pred_proba)
        }

        # Confusion matrix counts
        metrics.update({
            'tp': int(tp),
            'tn': int(tn),
            'fp': int(fp),
            'fn': int(fn)
        })

        return metrics

    def _plot_evaluation(self, y_true: np.ndarray, y_pred: np.ndarray,
                        y_pred_proba: np.ndarray, metrics: Dict):
        """Create evaluation plots."""
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # Confusion Matrix
        cm = confusion_matrix(y_true, y_pred)
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=axes[0, 0])
        axes[0, 0].set_title('Confusion Matrix')
        axes[0, 0].set_xlabel('Predicted')
        axes[0, 0].set_ylabel('Actual')

        # ROC Curve
        fpr, tpr, _ = roc_curve(y_true, y_pred_proba)
        roc_auc = auc(fpr, tpr)
        axes[0, 1].plot(fpr, tpr, color='darkorange', lw=2,
                       label=f'ROC curve (AUC = {roc_auc:.3f}')
        axes[0, 1].plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
        axes[0, 1].set_xlim([0.0, 1.0])
        axes[0, 1].set_ylim([0.0, 1.05])
        axes[0, 1].set_xlabel('False Positive Rate')
        axes[0, 1].set_ylabel('True Positive Rate')
        axes[0, 1].set_title('Receiver Operating Characteristic')
        axes[0, 1].legend(loc="lower right")

        # Precision-Recall Curve
        precision, recall, _ = precision_recall_curve(y_true, y_pred_proba)
        pr_auc = auc(recall, precision)
        axes[1, 0].plot(recall, precision, color='green', lw=2,
                       label=f'PR curve (AUC = {pr_auc:.3f}')
        axes[1, 0].set_xlim([0.0, 1.0])
        axes[1, 0].set_ylim([0.0, 1.05])
        axes[1, 0].set_xlabel('Recall')
        axes[1, 0].set_ylabel('Precision')
        axes[1, 0].set_title('Precision-Recall Curve')
        axes[1, 0].legend(loc="lower right")

        # Feature Importance (if available)
        if hasattr(self.model, 'feature_importances_'):
            if self.feature_names is None:
                feature_names = [f'Feature_{i}' for i in range(len(self.model.feature_importances_))]
            else:
                feature_names = self.feature_names

            importance = self.model.feature_importances_
            indices = np.argsort(importance)[::-1]

            axes[1, 1].bar(range(min(10, len(importance))),
                          importance[indices][:10])
            axes[1, 1].set_xticks(range(min(10, len(importance))))
            axes[1, 1].set_xticklabels([feature_names[i] for i in indices][:10],
                                      rotation=45, ha='right')
            axes[1, 1].set_title('Top 10 Feature Importances')
            axes[1, 1].set_ylabel('Importance')

        plt.tight_layout()
        plt.show()

        # Print metrics
        print("\n" + "="*60)
        print("CLASSIFICATION REPORT")
        print("="*60)
        print(classification_report(y_true, y_pred))

        print("\n" + "="*60)
        print("ADDITIONAL METRICS")
        print("="*60)
        for metric, value in metrics.items():
            if metric not in ['tp', 'tn', 'fp', 'fn']:
                print(f"{metric.upper():<15}: {value:.4f}")

        print(f"\nTP: {metrics['tp']}, TN: {metrics['tn']}, FP: {metrics['fp']}, FN: {metrics['fn']}")


def create_synthetic_imbalanced_data(n_samples: int = 1000,
                                   n_features: int = 20,
                                   imbalance_ratio: float = 0.1,
                                   random_state: int = 42) -> Tuple[pd.DataFrame, pd.Series]:
    """
    Create synthetic imbalanced dataset for testing.

    Args:
        n_samples: Total number of samples
        n_features: Number of features
        imbalance_ratio: Ratio of minority class to majority class
        random_state: Random seed

    Returns:
        Tuple of (features, labels)
    """
    np.random.seed(random_state)

    # Calculate samples per class
    n_minority = int(n_samples * imbalance_ratio)
    n_majority = n_samples - n_minority

    # Generate features
    X_minority = np.random.randn(n_minority, n_features) + 1  # Shifted for separation
    X_majority = np.random.randn(n_majority, n_features)

    X = np.vstack([X_minority, X_majority])

    # Generate informative features
    # Add some correlation to make classification easier
    for i in range(min(5, n_features)):
        X[:, i] += 0.5 * (X[:, i] > 0).astype(float)

    # Create labels (1 for minority, 0 for majority)
    y = np.hstack([np.ones(n_minority), np.zeros(n_majority)])

    # Shuffle the data
    indices = np.random.permutation(n_samples)
    X = X[indices]
    y = y[indices]

    # Convert to DataFrame for better handling
    feature_names = [f'feature_{i}' for i in range(n_features)]
    X_df = pd.DataFrame(X, columns=feature_names)
    y_series = pd.Series(y, name='target')

    return X_df, y_series


# Example usage function
def example_usage():
    """Demonstrate how to use the ImbalancedBinaryClassifier."""

    print("Creating synthetic imbalanced dataset...")
    X, y = create_synthetic_imbalanced_data(
        n_samples=1000,
        n_features=10,
        imbalance_ratio=0.1  # 10% minority class
    )

    print(f"Dataset shape: {X.shape}")
    print(f"Class distribution:\n{y.value_counts()}")

    # Test different methods
    methods = ['class_weight', 'smote', 'balanced_rf']

    results = {}

    for method in methods:
        print(f"\n{'='*60}")
        print(f"Testing method: {method}")
        print('='*60)

        # Create and train classifier
        clf = ImbalancedBinaryClassifier(method=method)
        clf.fit(X, y, validation_split=0.2)

        # Evaluate
        metrics = clf.evaluate()
        results[method] = metrics

    # Compare results
    print(f"\n{'='*80}")
    print("COMPARISON OF METHODS")
    print('='*80)

    comparison_df = pd.DataFrame({
        method: {k: v for k, v in results[method].items()
                if k in ['f1_score', 'precision', 'recall', 'roc_auc', 'pr_auc']}
        for method in results.keys()
    }).T

    print(comparison_df.round(4))

    return results


if __name__ == "__main__":
    example_usage()