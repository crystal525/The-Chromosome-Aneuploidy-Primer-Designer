# Machine Learning Toolkit for Imbalanced Genomic Data

This toolkit provides comprehensive tools for handling imbalanced binary classification problems commonly encountered in genomic and biomedical data analysis.

## Overview

The toolkit consists of three main modules:

1. **`imbalanced_classifier.py`** - Core classification algorithms for imbalanced data
2. **`data_preprocessing.py`** - Data preprocessing utilities for genomic data
3. **`evaluation_metrics.py`** - Comprehensive evaluation metrics for imbalanced classification

## Features

### Imbalanced Classification Methods

- **Class Weighting**: Automatic class weight balancing
- **SMOTE (Synthetic Minority Oversampling Technique)**: Oversampling minority class
- **ADASYN (Adaptive Synthetic Sampling)**: Adaptive oversampling
- **Random Undersampling**: Undersampling majority class
- **Balanced Random Forest**: Ensemble method with balanced bootstrapping
- **Easy Ensemble**: Ensemble of balanced classifiers

### Data Preprocessing

- **Feature Scaling**: Standard, MinMax, and Robust scaling
- **Feature Selection**: Statistical feature selection (F-test, mutual information)
- **Outlier Detection and Removal**: IQR-based outlier handling
- **Dimensionality Reduction**: PCA integration
- **Class Distribution Analysis**: Comprehensive imbalance assessment

### Evaluation Metrics

- **Standard Metrics**: Precision, Recall, F1-Score, Accuracy
- **Advanced Metrics**: Matthews Correlation Coefficient, Balanced Accuracy
- **Probability Metrics**: ROC AUC, Precision-Recall AUC, Brier Score
- **Threshold Optimization**: Find optimal decision thresholds
- **Bootstrap Confidence Intervals**: Statistical reliability assessment

## Installation

1. Install the required dependencies:
```bash
pip install -r requirements_ml.txt
```

2. Ensure you have the necessary system dependencies for the main ChAPDes pipeline (see main README.md)

## Quick Start

### Basic Usage

```python
from imbalanced_classifier import ImbalancedBinaryClassifier
from data_preprocessing import GenomicDataPreprocessor
from evaluation_metrics import ImbalancedClassificationMetrics

# Load your genomic data
# X should be a DataFrame with features (e.g., gene expressions, SNP data)
# y should be a Series with binary labels (0=normal, 1=disease)

# Preprocess data
preprocessor = GenomicDataPreprocessor(
    scaler_type='standard',
    feature_selection_method='f_classif',
    n_features=50
)
X_processed = preprocessor.fit_transform(X, y)

# Train classifier
clf = ImbalancedBinaryClassifier(method='smote')
clf.fit(X_processed, y)

# Make predictions
predictions = clf.predict(X_processed)
probabilities = clf.predict_proba(X_processed)

# Evaluate performance
metrics = clf.evaluate(X_processed, y)
```

### Advanced Usage

```python
# Comprehensive evaluation
from genomic_imbalanced_example import create_genomic_like_dataset

# Create synthetic genomic data for testing
X, y, feature_names = create_genomic_like_dataset(
    n_samples=2000,
    n_features=50,
    imbalance_ratio=0.15  # 15% minority class
)

# Compare multiple methods
methods = ['class_weight', 'smote', 'balanced_rf', 'easy_ensemble']
results = {}

for method in methods:
    clf = ImbalancedBinaryClassifier(method=method)
    clf.fit(X, y, validation_split=0.2)
    metrics = clf.evaluate()
    results[method] = metrics

# Find best performing method
best_method = max(results.keys(), key=lambda x: results[x]['f1_score'])
print(f"Best method: {best_method} with F1-score: {results[best_method]['f1_score']:.4f}")
```

## Integration with ChAPDes Pipeline

The ML toolkit can be integrated into the ChAPDes pipeline for enhanced primer design decisions:

```python
from ChAPDes import pipeline
from imbalanced_classifier import ImbalancedBinaryClassifier

# Run standard ChAPDes pipeline
pipeline(chrom='chr21', popul='EAS')

# Use ML to classify primer candidates based on genomic features
# (Example: classify primers as high/low specificity)

# Extract features from primer design results
primer_features = extract_primer_features(primer_results)

# Train classifier on historical primer performance data
clf = ImbalancedBinaryClassifier(method='balanced_rf')
clf.fit(primer_features_train, performance_labels_train)

# Predict performance of new primer candidates
predicted_performance = clf.predict(primer_features_new)
```

## Evaluation Metrics Guide

For imbalanced datasets, consider these key metrics:

1. **F1-Score**: Harmonic mean of precision and recall - good for imbalanced data
2. **Balanced Accuracy**: Average of sensitivity and specificity
3. **Precision-Recall AUC**: Better than ROC AUC for highly imbalanced data
4. **Matthews Correlation Coefficient**: Accounts for all confusion matrix elements

### Metric Selection by Imbalance Level

- **Mild Imbalance (15-30% minority)**: F1-score, Balanced Accuracy
- **Moderate Imbalance (5-15% minority)**: PR AUC, F1-score
- **Severe Imbalance (<5% minority)**: PR AUC, Matthews Correlation Coefficient

## Data Format Requirements

### Input Data Format

```python
# Features (X)
- pandas.DataFrame or numpy.ndarray
- Shape: (n_samples, n_features)
- Feature names should be meaningful (e.g., 'GENE_001', 'SNP_rs123')

# Labels (y)
- pandas.Series or numpy.ndarray
- Binary: 0 (majority/negative class), 1 (minority/positive class)
- Imbalanced: minority class should be < 50% of total samples
```

### Output Format

```python
# Predictions
- Binary classifications: numpy.ndarray of shape (n_samples,)
- Predicted probabilities: numpy.ndarray of shape (n_samples, 2)

# Evaluation Metrics
- Dictionary containing all calculated metrics
- Keys include: 'f1_score', 'precision', 'recall', 'roc_auc', etc.
```

## Performance Tips

1. **Feature Selection**: Use domain knowledge to select relevant genomic features
2. **Scaling**: Standard scaling works well for most genomic data
3. **Method Selection**:
   - Use SMOTE for moderate imbalance
   - Use class weighting for simple cases
   - Use ensemble methods for complex datasets
4. **Validation**: Always use stratified cross-validation
5. **Threshold Tuning**: Optimize probability thresholds for your specific use case

## Example Applications in Genomics

1. **Disease Classification**: Classify samples as disease/normal based on genomic markers
2. **Primer Specificity**: Predict primer specificity using sequence features
3. **SNP Classification**: Identify functional vs. neutral genetic variants
4. **Expression Analysis**: Classify tissue types or conditions from gene expression
5. **Drug Response**: Predict drug efficacy from genomic profiles

## Troubleshooting

### Common Issues

1. **Memory Issues**: For large genomic datasets, use undersampling or feature selection
2. **Overfitting**: Use cross-validation and regularization
3. **Poor Performance**: Check class distribution and consider different balancing methods

### Debugging Tips

```python
# Check class distribution
from data_preprocessing import analyze_class_distribution
analysis = analyze_class_distribution(y)
print(f"Imbalance ratio: {analysis['imbalance_ratio']:.3f}")

# Visualize feature distributions
from data_preprocessing import visualize_feature_distributions
visualize_feature_distributions(X, y, features=['important_gene_1', 'important_gene_2'])

# Check for outliers
from data_preprocessing import detect_outliers_iqr
outliers = detect_outliers_iqr(X)
print(f"Outlier summary: {outliers}")
```

## Citation

If you use this toolkit in your research, please cite:

```
ChAPDes: Chromosome Aneuploidy Primer Designer
[Your Name], [Year]
GitHub: https://github.com/yourusername/ChAPDes
```

## License

This toolkit follows the same license as the main ChAPDes project.

## Contributing

Contributions are welcome! Areas for improvement:
- Additional balancing techniques
- Genomic-specific preprocessing methods
- Integration with more ML frameworks
- Performance optimizations

---

*For more information about the main ChAPDes pipeline, see the main README.md file.*