#!/bin/bash

# Start Enrichment Analysis Server for scSAID
# This script starts the Dash server on port 8051

echo "======================================================"
echo "  scSAID Enrichment Analysis Server"
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
python3 -c "import dash, plotly, pandas, numpy, flask" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "Error: Required Python packages are missing"
    echo "Please install: pip install dash plotly pandas numpy flask"
    exit 1
fi

echo "Dependencies OK"
echo ""

# Change to script directory
cd "$SCRIPT_DIR"

# Check if data files exist
if [ ! -f "enrichment_results.json" ]; then
    echo "Error: enrichment_results.json not found"
    exit 1
fi

if [ ! -f "marker_genes_fine_map.json" ]; then
    echo "Error: marker_genes_fine_map.json not found"
    exit 1
fi

echo "Data files found"
echo ""

# Start the server
echo "Starting server on http://localhost:8051/enrichment/"
echo "Press Ctrl+C to stop"
echo "======================================================"
echo ""

python3 enrichment_web.py
