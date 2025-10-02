#!/bin/bash

# Setup script for ChAPDes Machine Learning Environment
# This script helps set up the ML dependencies for imbalanced classification

echo "Setting up ChAPDes Machine Learning Environment..."
echo "=" * 60

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo "Python 3 is required but not installed."
    echo "Please install Python 3.6 or higher and try again."
    exit 1
fi

# Check if pip is available
if ! command -v pip3 &> /dev/null; then
    echo "pip3 is required but not installed."
    echo "Please install pip3 and try again."
    exit 1
fi

echo "Installing Python dependencies..."

# Install core ML packages
pip3 install numpy pandas scikit-learn imbalanced-learn matplotlib seaborn

# Install development tools (optional)
echo "Installing development tools..."
pip3 install pytest jupyter

echo "Installation completed!"
echo ""
echo "To test the installation, run:"
echo "  python3 -c 'from imbalanced_classifier import ImbalancedBinaryClassifier; print(\"ML modules imported successfully!\")'"
echo ""
echo "To run the example:"
echo "  python3 genomic_imbalanced_example.py"
echo ""
echo "To run the ML integration example with ChAPDes:"
echo "  python3 ChAPDes.py --ml-example"
echo ""
echo "For more information, see README_ML.md"