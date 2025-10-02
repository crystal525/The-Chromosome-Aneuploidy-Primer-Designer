"""
Data Preprocessing Utilities for Imbalanced Classification

This module provides comprehensive preprocessing tools specifically designed
for handling imbalanced datasets in genomic and biomedical applications.
"""

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
from sklearn.feature_selection import SelectKBest, f_classif, mutual_info_classif
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import warnings
warnings.filterwarnings('ignore')
from typing import Dict, List, Tuple, Optional, Union
import matplotlib.pyplot as plt
import seaborn as sns


class GenomicDataPreprocessor:
    """
    Preprocessing pipeline specifically designed for genomic data with class imbalance.

    This class handles common preprocessing tasks such as:
    - Feature scaling and normalization
    - Feature selection
    - Dimensionality reduction
    - Data cleaning and validation
    """

    def __init__(self, scaler_type: str = 'standard', feature_selection_method: str = None,
                 n_features: int = None, pca_components: int = None, random_state: int = 42):
        """
        Initialize the preprocessor.

        Args:
            scaler_type: Type of scaler ('standard', 'minmax', 'robust')
            feature_selection_method: Method for feature selection ('f_classif', 'mutual_info')
            n_features: Number of top features to select
            pca_components: Number of PCA components (if None, no PCA)
            random_state: Random seed for reproducibility
        """
        self.scaler_type = scaler_type
        self.feature_selection_method = feature_selection_method
        self.n_features = n_features
        self.pca_components = pca_components
        self.random_state = random_state

        self.scaler = self._get_scaler()
        self.feature_selector = None
        self.pca = None if pca_components is None else PCA(n_components=pca_components, random_state=random_state)
        self.selected_features = None
        self.is_fitted = False

    def _get_scaler(self):
        """Get the appropriate scaler object."""
        scalers = {
            'standard': StandardScaler(),
            'minmax': MinMaxScaler(),
            'robust': RobustScaler()
        }
        return scalers.get(self.scaler_type, StandardScaler())

    def fit(self, X: Union[pd.DataFrame, np.ndarray],
            y: Union[pd.Series, np.ndarray] = None) -> 'GenomicDataPreprocessor':
        """
        Fit the preprocessor to the data.

        Args:
            X: Feature matrix
            y: Target vector (required for feature selection)

        Returns:
            Self for method chaining
        """
        # Convert to DataFrame if needed
        if isinstance(X, np.ndarray):
            # Assume feature names if not provided
            feature_names = [f'feature_{i}' for i in range(X.shape[1])]
            X = pd.DataFrame(X, columns=feature_names)

        self.feature_names_in = X.columns.tolist()
        X_array = X.values

        # Feature selection (if specified)
        if self.feature_selection_method and y is not None:
            if isinstance(y, pd.Series):
                y = y.values

            X_selected, self.selected_features = self._select_features(X_array, y)
        else:
            X_selected = X_array

        # Fit scaler
        self.scaler.fit(X_selected)

        # Fit PCA (if specified)
        if self.pca is not None:
            self.pca.fit(X_selected)

        self.is_fitted = True
        return self

    def transform(self, X: Union[pd.DataFrame, np.ndarray]) -> np.ndarray:
        """
        Transform the data using fitted preprocessor.

        Args:
            X: Feature matrix to transform

        Returns:
            Transformed feature matrix
        """
        if not self.is_fitted:
            raise ValueError("Preprocessor must be fitted before transforming data")

        # Convert to DataFrame if needed
        if isinstance(X, np.ndarray):
            feature_names = [f'feature_{i}' for i in range(X.shape[1])]
            X = pd.DataFrame(X, columns=feature_names)

        X_array = X.values

        # Feature selection (if fitted)
        if self.selected_features is not None:
            X_array = X_array[:, self.selected_features]

        # Scale features
        X_scaled = self.scaler.transform(X_array)

        # Apply PCA (if specified)
        if self.pca is not None:
            X_scaled = self.pca.transform(X_scaled)

        return X_scaled

    def fit_transform(self, X: Union[pd.DataFrame, np.ndarray],
                     y: Union[pd.Series, np.ndarray] = None) -> np.ndarray:
        """
        Fit and transform in one step.

        Args:
            X: Feature matrix
            y: Target vector (required for feature selection)

        Returns:
            Transformed feature matrix
        """
        return self.fit(X, y).transform(X)

    def _select_features(self, X: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, List[int]]:
        """Select top features based on the specified method."""
        if self.feature_selection_method == 'f_classif':
            selector = SelectKBest(f_classif, k=self.n_features)
        elif self.feature_selection_method == 'mutual_info':
            selector = SelectKBest(mutual_info_classif, k=self.n_features)
        else:
            raise ValueError(f"Unknown feature selection method: {self.feature_selection_method}")

        X_selected = selector.fit_transform(X, y)
        selected_indices = selector.get_support(indices=True)

        return X_selected, selected_indices.tolist()

    def get_feature_names_out(self) -> List[str]:
        """Get names of selected features."""
        if not self.is_fitted:
            raise ValueError("Preprocessor must be fitted to get feature names")

        if self.selected_features is not None:
            return [self.feature_names_in[i] for i in self.selected_features]
        else:
            return self.feature_names_in.copy()


def analyze_class_distribution(y: Union[pd.Series, np.ndarray]) -> Dict:
    """
    Analyze the class distribution in the target variable.

    Args:
        y: Target vector

    Returns:
        Dictionary containing class distribution analysis
    """
    if isinstance(y, pd.Series):
        class_counts = y.value_counts()
        class_props = y.value_counts(normalize=True)
    else:
        unique, counts = np.unique(y, return_counts=True)
        class_counts = pd.Series(counts, index=unique)
        class_props = pd.Series(counts / len(y), index=unique)

    minority_class = class_counts.idxmin()
    majority_class = class_counts.idxmax()
    imbalance_ratio = class_counts[minority_class] / class_counts[majority_class]

    analysis = {
        'class_counts': class_counts.to_dict(),
        'class_proportions': class_props.to_dict(),
        'minority_class': int(minority_class),
        'majority_class': int(majority_class),
        'imbalance_ratio': float(imbalance_ratio),
        'total_samples': len(y),
        'is_imbalanced': imbalance_ratio < 0.5  # Consider imbalanced if ratio < 0.5
    }

    return analysis


def detect_outliers_iqr(X: Union[pd.DataFrame, np.ndarray],
                      features: Optional[List[str]] = None,
                      multiplier: float = 1.5) -> Dict:
    """
    Detect outliers using IQR method.

    Args:
        X: Feature matrix
        features: List of feature names to analyze (if None, analyze all)
        multiplier: IQR multiplier for outlier detection

    Returns:
        Dictionary containing outlier information
    """
    if isinstance(X, np.ndarray):
        if features is None:
            features = [f'feature_{i}' for i in range(X.shape[1])]
        X = pd.DataFrame(X, columns=features)

    if features is None:
        features = X.columns.tolist()

    outlier_info = {}

    for feature in features:
        if feature in X.columns:
            Q1 = X[feature].quantile(0.25)
            Q3 = X[feature].quantile(0.75)
            IQR = Q3 - Q1

            lower_bound = Q1 - multiplier * IQR
            upper_bound = Q3 + multiplier * IQR

            outliers = X[(X[feature] < lower_bound) | (X[feature] > upper_bound)]
            outlier_info[feature] = {
                'n_outliers': len(outliers),
                'outlier_percentage': len(outliers) / len(X) * 100,
                'bounds': {'lower': lower_bound, 'upper': upper_bound}
            }

    return outlier_info


def remove_outliers_iqr(X: Union[pd.DataFrame, np.ndarray],
                       y: Union[pd.Series, np.ndarray] = None,
                       features: Optional[List[str]] = None,
                       multiplier: float = 1.5) -> Union[pd.DataFrame, Tuple[pd.DataFrame, pd.Series]]:
    """
    Remove outliers using IQR method.

    Args:
        X: Feature matrix
        y: Target vector (optional)
        features: List of feature names to process
        multiplier: IQR multiplier for outlier detection

    Returns:
        DataFrame with outliers removed, optionally with target vector
    """
    if isinstance(X, np.ndarray):
        if features is None:
            features = [f'feature_{i}' for i in range(X.shape[1])]
        X = pd.DataFrame(X, columns=features)

    if features is None:
        features = X.columns.tolist()

    # Create mask for non-outliers
    mask = np.ones(len(X), dtype=bool)

    for feature in features:
        if feature in X.columns:
            Q1 = X[feature].quantile(0.25)
            Q3 = X[feature].quantile(0.75)
            IQR = Q3 - Q1

            lower_bound = Q1 - multiplier * IQR
            upper_bound = Q3 + multiplier * IQR

            feature_mask = (X[feature] >= lower_bound) & (X[feature] <= upper_bound)
            mask &= feature_mask

    # Apply mask
    X_clean = X[mask]

    if y is not None:
        if isinstance(y, pd.Series):
            y_clean = y[mask]
        else:
            y_clean = y[mask]
        return X_clean, y_clean

    return X_clean


def visualize_class_distribution(y: Union[pd.Series, np.ndarray],
                               title: str = "Class Distribution"):
    """
    Visualize class distribution.

    Args:
        y: Target vector
        title: Plot title
    """
    if isinstance(y, np.ndarray):
        y = pd.Series(y)

    plt.figure(figsize=(10, 6))

    # Count plot
    plt.subplot(1, 2, 1)
    ax = y.value_counts().plot(kind='bar')
    plt.title(f'{title} - Counts')
    plt.xlabel('Class')
    plt.ylabel('Count')

    # Add value labels on bars
    for p in ax.patches:
        ax.annotate(str(p.get_height()), (p.get_x() + p.get_width() / 2., p.get_height()),
                   ha='center', va='center', xytext=(0, 10), textcoords='offset points')

    # Pie chart
    plt.subplot(1, 2, 2)
    y.value_counts().plot(kind='pie', autopct='%1.1f%%')
    plt.title(f'{title} - Proportions')
    plt.ylabel('')

    plt.tight_layout()
    plt.show()


def visualize_feature_distributions(X: Union[pd.DataFrame, np.ndarray],
                                  y: Union[pd.Series, np.ndarray],
                                  features: Optional[List[str]] = None,
                                  n_cols: int = 3):
    """
    Visualize feature distributions by class.

    Args:
        X: Feature matrix
        y: Target vector
        features: Features to visualize (if None, sample from available)
        n_cols: Number of columns in subplot grid
    """
    if isinstance(X, np.ndarray):
        if features is None:
            features = [f'feature_{i}' for i in range(min(9, X.shape[1]))]
        X = pd.DataFrame(X, columns=features)

    if features is None:
        features = X.columns.tolist()[:9]  # Limit to first 9 features

    n_features = len(features)
    n_rows = (n_features + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)

    if isinstance(y, np.ndarray):
        y = pd.Series(y)

    for i, feature in enumerate(features):
        row, col = i // n_cols, i % n_cols

        if feature in X.columns:
            # Box plot by class
            data_to_plot = [X[feature][y == cls] for cls in y.unique()]

            bp = axes[row, col].boxplot(data_to_plot, labels=[f'Class_{cls}' for cls in y.unique()])
            axes[row, col].set_title(f'{feature} Distribution by Class')
            axes[row, col].set_xlabel('Class')
            axes[row, col].set_ylabel(feature)

    # Hide empty subplots
    for i in range(n_features, n_rows * n_cols):
        row, col = i // n_cols, i % n_cols
        axes[row, col].set_visible(False)

    plt.tight_layout()
    plt.show()


def create_cross_validation_folds(X: Union[pd.DataFrame, np.ndarray],
                                y: Union[pd.Series, np.ndarray],
                                n_splits: int = 5,
                                stratify: bool = True) -> List[Tuple]:
    """
    Create stratified cross-validation folds for imbalanced data.

    Args:
        X: Feature matrix
        y: Target vector
        n_splits: Number of CV folds
        stratify: Whether to stratify folds

    Returns:
        List of (train_indices, val_indices) tuples
    """
    if isinstance(X, pd.DataFrame):
        X = X.values
    if isinstance(y, pd.Series):
        y = y.values

    if stratify:
        skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
        folds = list(skf.split(X, y))
    else:
        # Simple random split
        indices = np.random.permutation(len(X))
        fold_size = len(X) // n_splits
        folds = []

        for i in range(n_splits):
            start_idx = i * fold_size
            end_idx = (i + 1) * fold_size if i < n_splits - 1 else len(X)

            val_indices = indices[start_idx:end_idx]
            train_indices = np.concatenate([indices[:start_idx], indices[end_idx:]])

            folds.append((train_indices, val_indices))

    return folds


# Example usage function
def example_preprocessing():
    """Demonstrate preprocessing utilities."""

    # Create sample data
    from imbalanced_classifier import create_synthetic_imbalanced_data

    print("Creating synthetic imbalanced dataset...")
    X, y = create_synthetic_imbalanced_data(n_samples=1000, n_features=10, imbalance_ratio=0.1)

    print("Original class distribution:")
    print(y.value_counts())

    # Analyze class distribution
    analysis = analyze_class_distribution(y)
    print(f"\nClass distribution analysis: {analysis}")

    # Visualize class distribution
    visualize_class_distribution(y)

    # Detect outliers
    outlier_info = detect_outliers_iqr(X, multiplier=1.5)
    print("
Outlier analysis:")
    for feature, info in outlier_info.items():
        print(f"{feature}: {info['n_outliers']} outliers ({info['outlier_percentage']".2f"}%)")

    # Remove outliers
    X_clean, y_clean = remove_outliers_iqr(X, y, multiplier=1.5)
    print(f"\nData shape after outlier removal: {X_clean.shape}")

    # Feature selection and scaling
    preprocessor = GenomicDataPreprocessor(
        scaler_type='standard',
        feature_selection_method='f_classif',
        n_features=8
    )

    X_preprocessed = preprocessor.fit_transform(X_clean, y_clean)
    print(f"\nPreprocessed data shape: {X_preprocessed.shape}")

    selected_features = preprocessor.get_feature_names_out()
    print(f"Selected features: {selected_features}")

    # Visualize feature distributions
    visualize_feature_distributions(X_clean, y_clean, features=selected_features[:6])

    return X_preprocessed, y_clean


if __name__ == "__main__":
    example_preprocessing()