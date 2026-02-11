# Enrichment Analysis Integration

## Overview

The enrichment analysis feature provides gene set enrichment analysis for cell type marker genes using MSigDB v2026.1.Mm databases.

## Features

- **Interactive Cell Type Selection**: Choose from all available cell types
- **Multiple Databases**: GO Biological Process, GO Cellular Component, GO Molecular Function, Hallmark, etc.
- **Marker Gene Display**: Shows marker genes with their statistics
- **Enrichment Results Table**: Detailed table with p-values, FDR, scores, and overlapping genes
- **Visual Charts**: Bar charts showing top enriched terms
- **Warm Scientific Aesthetic**: Matches the scSAID design system

## Architecture

- **Backend**: Python Dash + Flask server (port 8051)
- **Frontend**: Embedded iframe in details.jsp
- **Data**: Pre-computed enrichment results and marker genes in JSON format

## Setup Instructions

### 1. Install Dependencies

```bash
pip install dash plotly pandas numpy flask
```

### 2. Start the Enrichment Server

Navigate to the enrichment resources directory:

```bash
cd src/main/webapp/enrichment_resources/gmt/web
```

Run the startup script:

```bash
./start_enrichment.sh
```

Or start Python directly:

```bash
python3 enrichment_web.py
```

The server will start on `http://localhost:8051/enrichment/`

### 3. Access in Browser

1. Start your Tomcat server
2. Start the enrichment server (as above)
3. Navigate to any dataset details page
4. Scroll to the "Enrichment Analysis" section

## Files Structure

```
enrichment_resources/
├── README.md (this file)
└── gmt/
    ├── m5.go.bp.v2026.1.Mm.symbols.gmt
    ├── m5.go.cc.v2026.1.Mm.symbols.gmt
    ├── m5.go.mf.v2026.1.Mm.symbols.gmt
    ├── m8.all.v2026.1.Mm.symbols.gmt
    ├── mh.all.v2026.1.Mm.symbols.gmt
    └── web/
        ├── enrichment_web.py          # Main Dash application
        ├── enrichment.py               # Original dark theme version
        ├── enrichment_results.json     # Pre-computed enrichment data
        ├── marker_genes_fine_map.json  # Marker genes data
        └── start_enrichment.sh         # Startup script
```

## Design Integration

The enrichment interface uses the scSAID warm scientific aesthetic:

### Colors
- Primary: `#e8927c` (Coral)
- Secondary: `#d4a574` (Bronze)
- Background: `#faf8f5` (Warm beige)
- Text: `#1a2332` (Dark navy)
- Border: `#e5e0d8` (Light taupe)

### Typography
- Headers: Cormorant Garamond (serif)
- Body: Source Sans 3 (sans-serif)
- Mono: JetBrains Mono (code/genes)

### Components
- Rounded corners (12px)
- Subtle shadows
- Clean borders
- Responsive layout

## Usage

1. **Select Cell Type**: Choose from the dropdown of available cell types
2. **Select Database**: Pick an enrichment database (GO_BP, GO_CC, etc.)
3. **View Results**:
   - Marker genes with statistics
   - Enrichment table with terms, p-values, FDR, and genes
   - Bar chart of top enriched terms
4. **Interpret Colors**:
   - Coral (●): p < 0.001 (highly significant)
   - Bronze (●): p < 0.01 (significant)
   - Taupe (●): p < 0.05 (marginally significant)
   - Gray (●): p ≥ 0.05 (not significant)

## Troubleshooting

### Server won't start
- Check Python 3 is installed: `python3 --version`
- Install dependencies: `pip install dash plotly pandas numpy flask`
- Check port 8051 is not in use: `lsof -i :8051`

### iframe not loading
- Verify enrichment server is running
- Check browser console for errors
- Ensure both servers are accessible from the same network

### No data showing
- Verify JSON files exist and are valid
- Check server console for error messages
- Ensure file paths are correct

## Development

To modify the enrichment interface:

1. Edit `enrichment_web.py` for backend changes
2. Restart the server to see changes
3. Edit `details.css` for iframe container styling
4. Edit `details.jsp` for iframe integration

## Production Deployment

For production deployment:

1. Use a process manager (systemd, supervisor, pm2)
2. Configure proper logging
3. Set up reverse proxy (nginx) if needed
4. Consider containerization (Docker)
5. Implement proper error handling and monitoring

## Credits

- MSigDB v2026.1.Mm from Broad Institute
- Dash framework from Plotly
- Design system matches scSAID warm scientific aesthetic
