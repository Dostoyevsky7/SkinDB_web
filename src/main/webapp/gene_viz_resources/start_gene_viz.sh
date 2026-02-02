#!/bin/bash

# Start Gene Expression Visualization Server for scSAID
# This script starts the Dash server on port 8053

echo "======================================================"
echo "  scSAID Gene Expression Visualization Server"
echo "======================================================"
echo ""

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo "Error: Python 3 is not installed or not in PATH"
    exit 1
fi

# Check if required packages are installed
echo "Checking dependencies..."
python3 -c "import dash, plotly, pandas, numpy, scanpy, scipy, flask" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "Error: Required Python packages are missing"
    echo "Please install: pip install dash plotly pandas numpy scanpy scipy flask"
    exit 1
fi

echo "Dependencies OK"
echo ""

# Change to script directory
cd "$SCRIPT_DIR"

# Check if integrated dataset exists
if [ ! -f "integrated.h5ad" ]; then
    echo "Warning: integrated.h5ad not found"
    echo "The server will start but visualization will not work until you create the integrated dataset"
    echo "See README.md for instructions on creating integrated.h5ad"
    echo ""
fi

# Start the server
echo "Starting server on http://localhost:8053/gene-viz/"
echo "Press Ctrl+C to stop"
echo "======================================================"
echo ""

python3 gene_expression_viz.py
