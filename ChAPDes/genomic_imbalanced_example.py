#!/usr/bin/env python3
"""
Comprehensive Example: Imbalanced Binary Classification for Genomic Data

This script demonstrates how to use the imbalanced classification toolkit
with genomic data, including data preprocessing, model training, and evaluation.
"""

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.datasets import make_classification

# Import our custom modules
from imbalanced_classifier import ImbalancedBinaryClassifier, create_synthetic_imbalanced_data
from data_preprocessing import (
    GenomicDataPreprocessor, analyze_class_distribution,
    visualize_class_distribution, detect_outliers_iqr, remove_outliers_iqr
)
from evaluation_metrics import ImbalancedClassificationMetrics

import warnings
warnings.filterwarnings('ignore')


def create_genomic_like_dataset(n_samples: int = 2000,
                              n_features: int = 50,
                              n_informative: int = 20,
                              n_redundant: int = 10,
                              imbalance_ratio: float = 0.15,
                              random_state: int = 42) -> tuple:
    """
    Create a synthetic dataset that mimics genomic data characteristics.

    Args:
        n_samples: Total number of samples
        n_features: Total number of features (e.g., genes, SNPs)
        n_informative: Number of informative features
        n_redundant: Number of redundant features
        imbalance_ratio: Ratio of minority class
        random_state: Random seed

    Returns:
        Tuple of (features, labels, feature_names)
    """
    # Calculate class samples
    n_minority = int(n_samples * imbalance_ratio)
    n_majority = n_samples - n_minority

    # Generate imbalanced dataset
    X, y = make_classification(
        n_samples=n_samples,
        n_features=n_features,
        n_informative=n_informative,
        n_redundant=n_redundant,
        n_clusters_per_class=1,
        weights=[1 - imbalance_ratio, imbalance_ratio],  # Majority, minority
        flip_y=0.05,  # Add some noise
        class_sep=1.0,  # Good separation for demonstration
        random_state=random_state
    )

    # Create meaningful feature names (like genomic coordinates or gene names)
    feature_names = []
    for i in range(n_features):
        if i < n_informative:
            feature_names.append(f"INFORMATIVE_GENE_{i+1"03d"}")
        elif i < n_informative + n_redundant:
            feature_names.append(f"REDUNDANT_GENE_{(i-n_informative)+1"03d"}")
        else:
            feature_names.append(f"NOISE_FEATURE_{(i-n_informative-n_redundant)+1"03d"}")

    # Convert to DataFrame for better handling
    X_df = pd.DataFrame(X, columns=feature_names)
    y_series = pd.Series(y, name='disease_status')

    return X_df, y_series, feature_names


def run_comprehensive_analysis():
    """Run a comprehensive analysis of imbalanced classification methods."""

    print("=" * 80)
    print("GENOMIC IMBALANCED CLASSIFICATION EXAMPLE")
    print("=" * 80)

    # Create genomic-like dataset
    print("Creating synthetic genomic dataset...")
    X, y, feature_names = create_genomic_like_dataset(
        n_samples=2000,
        n_features=50,
        n_informative=20,
        imbalance_ratio=0.15  # 15% minority class (e.g., disease cases)
    )

    print(f"Dataset shape: {X.shape}")
    print(f"Minority class ratio: {y.sum() / len(y)".1%"}")

    # Analyze class distribution
    print("\n" + "-" * 60)
    print("CLASS DISTRIBUTION ANALYSIS")
    print("-" * 60)

    class_analysis = analyze_class_distribution(y)
    print(f"Class counts: {class_analysis['class_counts']}")
    print(f"Imbalance ratio: {class_analysis['imbalance_ratio']".3f"}")
    print(f"Is imbalanced: {class_analysis['is_imbalanced']}")

    # Visualize class distribution
    visualize_class_distribution(y, "Disease Status Distribution")

    # Split data
    print("\n" + "-" * 60)
    print("DATA SPLITTING")
    print("-" * 60)

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )

    print(f"Training set: {X_train.shape}, {y_train.sum()} positive cases")
    print(f"Test set: {X_test.shape}, {y_test.sum()} positive cases")

    # Data preprocessing
    print("\n" + "-" * 60)
    print("DATA PREPROCESSING")
    print("-" * 60)

    # Detect and remove outliers
    outlier_info = detect_outliers_iqr(X_train, multiplier=1.5)
    outlier_features = [feat for feat, info in outlier_info.items()
                       if info['outlier_percentage'] > 5]  # Remove features with >5% outliers

    print(f"Features with significant outliers: {len(outlier_features)}")

    # Remove outliers if any significant ones found
    if outlier_features:
        X_train_clean, y_train_clean = remove_outliers_iqr(
            X_train, y_train, features=outlier_features, multiplier=1.5
        )
        print(f"Data shape after outlier removal: {X_train_clean.shape}")
    else:
        X_train_clean, y_train_clean = X_train.copy(), y_train.copy()

    # Set up preprocessor
    preprocessor = GenomicDataPreprocessor(
        scaler_type='standard',
        feature_selection_method='f_classif',
        n_features=30  # Select top 30 features
    )

    # Fit preprocessor on training data
    X_train_processed = preprocessor.fit_transform(X_train_clean, y_train_clean)

    # Transform test data
    X_test_processed = preprocessor.transform(X_test)

    print(f"Processed training data shape: {X_train_processed.shape}")
    print(f"Selected features: {len(preprocessor.get_feature_names_out())}")

    # Compare different classification methods
    print("\n" + "-" * 60)
    print("MODEL COMPARISON")
    print("-" * 60)

    methods = [
        'class_weight',     # Class weighting
        'smote',           # SMOTE oversampling
        'balanced_rf',     # Balanced Random Forest
        'easy_ensemble'    # Easy Ensemble
    ]

    results = {}

    for method in methods:
        print(f"\n--- Training {method.upper()} ---")

        # Create and train classifier
        clf = ImbalancedBinaryClassifier(method=method, random_state=42)
        clf.fit(X_train_processed, y_train_clean, validation_split=0.2)

        # Evaluate on test set
        test_metrics = clf.evaluate(X_test_processed, y_test, plot=False)

        results[method] = test_metrics

        # Print key metrics
        print(f"F1-Score: {test_metrics['f1_score']".4f"}")
        print(f"Precision: {test_metrics['precision']".4f"}")
        print(f"Recall: {test_metrics['recall']".4f"}")
        print(f"ROC AUC: {test_metrics.get('roc_auc', 'N/A')".4f"}")

    # Detailed evaluation of best method
    print("\n" + "-" * 60)
    print("DETAILED EVALUATION OF BEST METHOD")
    print("-" * 60)

    # Find best method based on F1-score
    best_method = max(results.keys(), key=lambda x: results[x]['f1_score'])
    print(f"Best method: {best_method}")

    # Retrain best method on full training data
    best_clf = ImbalancedBinaryClassifier(method=best_method, random_state=42)
    best_clf.fit(X_train_processed, y_train_clean, validation_split=0.2)

    # Comprehensive evaluation
    final_metrics = best_clf.evaluate(X_test_processed, y_test, plot=True)

    # Print detailed metrics
    ImbalancedClassificationMetrics.print_metrics_summary(final_metrics, best_method)

    # Additional analysis
    print("\n" + "-" * 60)
    print("ADDITIONAL ANALYSIS")
    print("-" * 60)

    # Threshold optimization
    y_test_proba = best_clf.predict_proba(X_test_processed)[:, 1]

    optimal_threshold, threshold_metrics = ImbalancedClassificationMetrics.find_optimal_threshold(
        y_test, y_test_proba, metric='f1'
    )

    print(f"Optimal threshold for F1-score: {optimal_threshold".3f"}")
    print(f"Metrics at optimal threshold: {threshold_metrics}")

    # Feature importance (if available)
    if hasattr(best_clf.model, 'feature_importances_'):
        feature_names_out = preprocessor.get_feature_names_out()
        importance = best_clf.model.feature_importances_

        # Get top 10 features
        top_indices = np.argsort(importance)[::-1][:10]
        print("\nTop 10 Most Important Features:")
        for i, idx in enumerate(top_indices):
            print(f"{i+1"2d"}. {feature_names_out[idx]}: {importance[idx]".4f"}")

    return results, final_metrics


def demonstrate_cross_validation():
    """Demonstrate cross-validation for imbalanced data."""

    print("\n" + "=" * 80)
    print("CROSS-VALIDATION EXAMPLE")
    print("=" * 80)

    # Create smaller dataset for faster CV
    X, y, _ = create_genomic_like_dataset(
        n_samples=1000,
        n_features=30,
        imbalance_ratio=0.1
    )

    # Simple preprocessing
    preprocessor = GenomicDataPreprocessor(scaler_type='standard')
    X_processed = preprocessor.fit_transform(X, y)

    # Create cross-validation folds
    from data_preprocessing import create_cross_validation_folds

    folds = create_cross_validation_folds(X, y, n_splits=5, stratify=True)

    print(f"Created {len(folds)} CV folds")

    # Evaluate different methods using CV
    methods = ['class_weight', 'smote', 'balanced_rf']
    cv_results = {}

    for method in methods:
        print(f"\n--- Cross-validating {method.upper()} ---")

        fold_metrics = []

        for fold_idx, (train_idx, val_idx) in enumerate(folds):
            # Split data
            X_train_fold = X_processed[train_idx]
            y_train_fold = y.iloc[train_idx]
            X_val_fold = X_processed[val_idx]
            y_val_fold = y.iloc[val_idx]

            # Train and evaluate
            clf = ImbalancedBinaryClassifier(method=method, random_state=42)
            clf.fit(X_train_fold, y_train_fold, validation_split=0.0)  # No internal validation

            # Get predictions on validation fold
            y_val_pred = clf.predict(X_val_fold)
            y_val_proba = clf.predict_proba(X_val_fold)[:, 1]

            # Calculate metrics
            fold_metric = ImbalancedClassificationMetrics.calculate_all_metrics(
                y_val_fold, y_val_pred, y_val_proba
            )
            fold_metrics.append(fold_metric['f1_score'])

        # Average across folds
        cv_results[method] = {
            'mean_f1': np.mean(fold_metrics),
            'std_f1': np.std(fold_metrics),
            'fold_scores': fold_metrics
        }

        print(f"CV F1-Score: {np.mean(fold_metrics)".4f"} ± {np.std(fold_metrics):".4f"")

    # Compare CV results
    print("\n" + "-" * 40)
    print("CROSS-VALIDATION COMPARISON")
    print("-" * 40)

    for method, result in cv_results.items():
        print(f"{method.upper()"<15"}: {result['mean_f1']".4f"} ± {result['std_f1']:".4f"")

    return cv_results


def main():
    """Main function to run the complete example."""

    try:
        # Run comprehensive analysis
        results, final_metrics = run_comprehensive_analysis()

        # Run cross-validation example
        cv_results = demonstrate_cross_validation()

        print("\n" + "=" * 80)
        print("EXAMPLE COMPLETED SUCCESSFULLY")
        print("=" * 80)
        print("The imbalanced classification toolkit is ready for genomic data analysis!")

        return results, final_metrics, cv_results

    except Exception as e:
        print(f"Error in example: {e}")
        import traceback
        traceback.print_exc()
        return None, None, None


if __name__ == "__main__":
    main()