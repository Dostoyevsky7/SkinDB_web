# Gene Expression UMAP Visualization

## Overview

The Gene Expression UMAP Visualization feature allows users to visualize gene expression patterns on an integrated UMAP plot across all datasets. Users can select one or multiple genes from search results and see their expression overlaid on the unified cell landscape.

## Features

### Core Functionality
- **Integrated UMAP Display**: Unified visualization of all datasets
- **Multi-Gene Selection**: Select multiple genes to visualize simultaneously
- **Expression Overlay**: Color-coded expression levels on UMAP
- **Interactive Exploration**:
  - Hover for cell-level details (cell type, dataset, expression)
  - Pan and zoom capabilities
  - Export as high-resolution PNG
- **Real-time Updates**: Instant visualization upon gene selection

### Visualization Details
- **Color Scale**: Gradient from light (low expression) to deep red (high expression)
- **Expression Aggregation**: When multiple genes selected, shows mean expression
- **Cell Metadata**: Displays cell type and dataset origin on hover
- **Publication-Ready**: High-resolution export with adjustable dimensions

## Data Requirements

### Integrated Dataset

The visualization requires an **integrated h5ad file** containing:
- All datasets merged into a single AnnData object
- UMAP coordinates (`X_umap` in `obsm`)
- Gene expression matrix (`X`)
- Cell type annotations (`cell_type` in `obs`)
- Dataset identifiers (`dataset` in `obs`)

**File Location**: `gene_viz_resources/integrated.h5ad`

### Creating the Integrated Dataset

If you don't have an integrated dataset yet, create one using this Python script:

```python
import scanpy as sc
import pandas as pd
import os

# Load dataset information
info_df = pd.read_excel("WEB-INF/IntegrateTable.xlsx")

# Path to individual h5ad files
h5ad_dir = "path/to/your/h5ad/files"

# Load all datasets
adatas = []
for idx, row in info_df.iterrows():
    said = row['SAID']
    h5ad_path = os.path.join(h5ad_dir, f"{said}.h5ad")

    if os.path.exists(h5ad_path):
        adata = sc.read_h5ad(h5ad_path)
        adata.obs['dataset'] = said
        adata.obs['gse'] = row['GSE']
        adatas.append(adata)
        print(f"Loaded {said}: {adata.n_obs} cells")

print(f"\nTotal datasets loaded: {len(adatas)}")

# Concatenate all datasets
integrated = sc.concat(adatas, join='outer', label='batch', keys=[a.obs['dataset'][0] for a in adatas])
print(f"Integrated dataset: {integrated.n_obs} cells, {integrated.n_vars} genes")

# Preprocessing for integration
sc.pp.normalize_total(integrated, target_sum=1e4)
sc.pp.log1p(integrated)
sc.pp.highly_variable_genes(integrated, n_top_genes=2000, batch_key='batch')

# Integration using Harmony (or other methods)
import scanpy.external as sce
sce.pp.harmony_integrate(integrated, 'batch')

# Compute UMAP on integrated space
sc.pp.neighbors(integrated, use_rep='X_pca_harmony')
sc.tl.umap(integrated)

# Save integrated dataset
integrated.write('gene_viz_resources/integrated.h5ad')
print(f"\nSaved integrated dataset to integrated.h5ad")
print(f"UMAP coordinates shape: {integrated.obsm['X_umap'].shape}")
```

**Alternative Integration Methods**:
- **Harmony**: Fast, works well for batch correction
- **Scanorama**: Good for diverse datasets
- **BBKNN**: K-nearest neighbor graph integration
- **scVI**: Deep learning-based integration

**Installation for Harmony**:
```bash
pip install scanpy-scripts harmonypy
```

### Required Metadata

Ensure your integrated dataset has these annotations:

- `adata.obs['cell_type']`: Cell type labels
- `adata.obs['dataset']`: Dataset identifier (SAID)
- `adata.obsm['X_umap']`: UMAP coordinates (n_cells × 2)

## Installation

### Dependencies

```bash
pip install dash plotly pandas numpy scanpy scipy flask
```

### File Structure

```
gene_viz_resources/
├── README.md (this file)
├── gene_expression_viz.py     # Main application
├── start_gene_viz.sh           # Startup script
└── integrated.h5ad             # Integrated dataset (you create this)
```

## Usage

### 1. Create Integrated Dataset

Follow the instructions above to create `integrated.h5ad`

### 2. Start the Visualization Server

```bash
cd src/main/webapp/gene_viz_resources
./start_gene_viz.sh
```

Or start Python directly:
```bash
python3 gene_expression_viz.py
```

Server starts on: `http://localhost:8053/gene-viz/`

### 3. Access via Gene Search Page

1. Start Tomcat server
2. Start visualization server (as above)
3. Navigate to gene search page: `http://localhost:8080/gene-search.jsp`
4. Search for genes (e.g., "KRT14", "COL1A1")
5. Select one or more genes from results
6. Click "Visualize Expression on UMAP"
7. View integrated UMAP with expression overlay

### 4. Using the Visualization

**Select Genes**:
- Check boxes next to genes in search results
- Select 1-10 genes (more genes = averaged expression)

**Interpret Plot**:
- **Color**: Expression level (light = low, red = high)
- **Hover**: See cell details (type, dataset, expression value)
- **Zoom**: Click and drag to zoom into regions
- **Export**: Use camera icon to save high-res image

**Multiple Genes**:
- Shows mean expression across selected genes
- Useful for gene signatures or pathway genes

## Design Integration

### Warm Scientific Aesthetic

**Colors**:
- Background: `#faf8f5` (Warm beige)
- Text: `#1a2332` (Dark navy)
- Borders: `#e5e0d8` (Light taupe)
- Expression gradient:
  - Low: `#f5f3f0` (Very light)
  - Mid: `#e8927c` (Coral - brand color)
  - High: `#d63031` (Deep red)

**Typography**:
- Headers: Cormorant Garamond (serif)
- Body: Urbanist (sans-serif)

**Components**:
- Rounded corners (8-12px)
- Subtle shadows
- Smooth hover effects

## API Endpoints

### `/gene-viz/?genes=GENE1,GENE2,...`
Main visualization interface

**Parameters**:
- `genes`: Comma-separated gene names (e.g., `KRT14,COL1A1`)

**Returns**: Interactive Dash application with UMAP

**Example**:
```
http://localhost:8053/gene-viz/?genes=KRT14,COL1A1,KRT5
```

### `/api/genes` (internal)
Get list of all available genes

**Returns**: JSON with gene names and count

### `/api/gene-info/<gene_name>` (internal)
Get statistics for a specific gene

**Parameters**:
- `gene_name`: Gene symbol (case-insensitive)

**Returns**: JSON with expression stats
```json
{
  "gene": "KRT14",
  "mean_expression": 2.34,
  "max_expression": 8.92,
  "pct_expressed": 45.2,
  "n_cells": 125000
}
```

## Troubleshooting

### Server won't start
- Check Python 3 installed: `python3 --version`
- Install dependencies: `pip install dash plotly pandas numpy scanpy scipy flask`
- Check port 8053 available: `lsof -i :8053`

### "Integrated dataset not available"
- Create `integrated.h5ad` following instructions above
- Ensure file is in `gene_viz_resources/` directory
- Check file permissions: `ls -l integrated.h5ad`

### "Genes not found"
- Verify gene names are uppercase (e.g., "KRT14" not "krt14")
- Check genes exist in integrated dataset
- Gene symbols should match mouse nomenclature

### Slow visualization
- Large datasets (>500k cells) may take 5-10 seconds to load
- Consider downsampling integrated dataset for faster performance
- Use representative subset for visualization if needed

### UMAP looks wrong
- Ensure integration was performed correctly
- Check that `X_umap` coordinates exist: `adata.obsm['X_umap']`
- Verify UMAP was computed on integrated space (not raw PCA)

## Performance Optimization

**For Large Datasets** (>500,000 cells):
- Consider downsampling to 100k-200k cells for visualization
- Use representative sampling across cell types
- Pre-compute UMAP on full dataset, sample for display

**Memory Usage**:
- Server loads integrated dataset into memory (~5-10 GB for large datasets)
- Restart server to clear cache if needed
- Monitor with: `ps aux | grep python`

**Speed Tips**:
- Use Scattergl (already implemented) for >100k cells
- Pre-filter genes to reduce dataset size
- Cache UMAP coordinates separately

## Integration Algorithm Details

### Harmony Integration (Recommended)

**Advantages**:
- Fast (minutes for 100k+ cells)
- Preserves biological variation
- Good batch effect removal
- No hyperparameter tuning needed

**How it works**:
1. Compute PCA on log-normalized data
2. Iteratively correct PC embeddings
3. Soft clustering and linear correction
4. UMAP on corrected space

### Alternative: Scanorama

```python
import scanorama

# After loading adatas
integrated = scanorama.integrate(adatas, return_dimred=True)
```

**Advantages**:
- Handles very different datasets
- Panorama matching approach
- Good for cross-study integration

### Alternative: scVI

```python
import scvi

adata.layers['counts'] = adata.X.copy()
scvi.model.SCVI.setup_anndata(adata, layer='counts', batch_key='batch')
model = scvi.model.SCVI(adata)
model.train()
adata.obsm['X_scvi'] = model.get_latent_representation()

sc.pp.neighbors(adata, use_rep='X_scvi')
sc.tl.umap(adata)
```

**Advantages**:
- Deep learning-based
- Best for complex batch effects
- Probabilistic framework

## Development

### Adding New Features

**Custom Color Scales**:
Edit the `colorscale` parameter in `create_umap_plot()`:
```python
colorscale=[
    [0, '#ffffff'],    # Low
    [0.5, '#e8927c'],  # Mid
    [1, '#d63031']     # High
]
```

**Additional Metadata Display**:
Modify hover text in `create_umap_plot()`:
```python
hover_text = [
    f"Cell: {i}<br>Type: {ct}<br>Dataset: {ds}<br>Expression: {exp:.3f}"
    for i, (ct, ds, exp) in enumerate(zip(...))
]
```

**Gene Set Visualization**:
Extend to show pathway scores instead of individual genes.

### Testing

**Test with sample genes**:
```bash
# Single gene
curl "http://localhost:8053/gene-viz/?genes=KRT14"

# Multiple genes
curl "http://localhost:8053/gene-viz/?genes=KRT14,COL1A1,KRT5"

# Get available genes
curl "http://localhost:8053/api/genes"

# Get gene info
curl "http://localhost:8053/api/gene-info/KRT14"
```

## Citation

When using this feature, please cite:

- **Scanpy**: Wolf et al. (2018) Genome Biology
- **Harmony**: Korsunsky et al. (2019) Nature Methods
- **UMAP**: McInnes et al. (2018) arXiv

## Future Enhancements

Planned features:
- Side-by-side comparison of multiple genes
- Cell density overlays
- Custom gene signatures
- Export expression values
- Statistical testing across cell types
- Integration quality metrics
- 3D UMAP visualization
- Animation between genes

## Support

For issues or questions:
1. Check this documentation
2. Verify integrated.h5ad exists and is valid
3. Review console logs for errors
4. Check GitHub issues

---

**Version**: 1.0.0
**Last Updated**: 2026-02-02
**Compatibility**: scSAID v1.0+
