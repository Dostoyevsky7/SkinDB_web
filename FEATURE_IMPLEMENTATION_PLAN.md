# scSAID Advanced Features - Comprehensive Implementation Plan

## Table of Contents
1. [Gene Information & External Links](#1-gene-information--external-links)
2. [Transcription Factor Activity Inference](#2-transcription-factor-activity-inference)
3. [Pathway Enrichment Analysis](#3-pathway-enrichment-analysis)
4. [Cell-Cell Communication Enhancement](#4-cell-cell-communication-enhancement)
5. [Trajectory & Pseudotime Analysis](#5-trajectory--pseudotime-analysis)
6. [Gene Regulatory Networks](#6-gene-regulatory-networks)
7. [Interactive Visualization Suite](#7-interactive-visualization-suite)
8. [Comparative Analysis Tools](#8-comparative-analysis-tools)
9. [Advanced Search & Filtering](#9-advanced-search--filtering)
10. [Cell Type Annotation & Prediction](#10-cell-type-annotation--prediction)
11. [REST API & Data Access](#11-rest-api--data-access)
12. [User Workspace & Saved Analyses](#12-user-workspace--saved-analyses)
13. [Custom Dataset Upload & Analysis](#13-custom-dataset-upload--analysis)
14. [Publication-Ready Figure Generation](#14-publication-ready-figure-generation)

---

## 1. Gene Information & External Links

### **Feature Description**
Provide comprehensive external database links for every gene, enabling researchers to quickly access detailed information.

### **Implementation Plan**

#### **1.1 Frontend Changes**

**File: `src/main/webapp/JSP/gene-details.jsp`** (NEW)
```jsp
<%@ page contentType="text/html;charset=UTF-8" %>
<!DOCTYPE html>
<html>
<head>
    <title>Gene: ${geneName} - scSAID</title>
    <link rel="stylesheet" href="CSS/gene-details.css">
</head>
<body>
    <div class="gene-header">
        <h1>${geneName}</h1>
        <div class="gene-links">
            <a href="https://www.ncbi.nlm.nih.gov/gene/?term=${geneName}" target="_blank">
                <img src="images/ncbi-logo.png" alt="NCBI Gene"/>
                NCBI Gene
            </a>
            <a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=${geneName}" target="_blank">
                <img src="images/genecards-logo.png" alt="GeneCards"/>
                GeneCards
            </a>
            <a href="https://www.ensembl.org/Multi/Search/Results?q=${geneName}" target="_blank">
                <img src="images/ensembl-logo.png" alt="Ensembl"/>
                Ensembl
            </a>
            <a href="https://www.proteinatlas.org/search/${geneName}" target="_blank">
                <img src="images/hpa-logo.png" alt="Human Protein Atlas"/>
                Protein Atlas
            </a>
            <a href="https://gtexportal.org/home/gene/${geneName}" target="_blank">
                <img src="images/gtex-logo.png" alt="GTEx"/>
                GTEx Portal
            </a>
            <a href="https://www.uniprot.org/uniprot/?query=${geneName}" target="_blank">
                <img src="images/uniprot-logo.png" alt="UniProt"/>
                UniProt
            </a>
            <a href="https://string-db.org/network/${geneName}" target="_blank">
                <img src="images/string-logo.png" alt="STRING"/>
                STRING
            </a>
        </div>
    </div>

    <div class="gene-expression-section">
        <h2>Expression Across Datasets</h2>
        <div id="expressionVisualization"></div>
    </div>

    <div class="gene-metadata">
        <h2>Gene Information</h2>
        <table>
            <tr><td>Gene Symbol:</td><td>${geneInfo.symbol}</td></tr>
            <tr><td>Gene ID:</td><td>${geneInfo.ensemblId}</td></tr>
            <tr><td>Description:</td><td>${geneInfo.description}</td></tr>
            <tr><td>Chromosome:</td><td>${geneInfo.chromosome}</td></tr>
            <tr><td>Location:</td><td>${geneInfo.location}</td></tr>
        </table>
    </div>
</body>
</html>
```

**Modify: `src/main/webapp/gene-search.jsp`**
- Add clickable gene names that link to gene details page
- Add quick-access external link icons next to each gene

#### **1.2 Backend Changes**

**File: `src/main/java/Servlet/GeneDetailsServlet.java`** (NEW)
```java
@WebServlet("/gene-details")
public class GeneDetailsServlet extends HttpServlet {
    protected void doGet(HttpServletRequest request, HttpServletResponse response) {
        String geneName = request.getParameter("gene");

        // Fetch gene information from Ensembl REST API
        GeneInfo geneInfo = fetchGeneInfo(geneName);

        // Get expression data across all datasets
        List<ExpressionData> expressionData = getGeneExpression(geneName);

        request.setAttribute("geneName", geneName);
        request.setAttribute("geneInfo", geneInfo);
        request.setAttribute("expressionData", expressionData);

        request.getRequestDispatcher("gene-details.jsp").forward(request, response);
    }

    private GeneInfo fetchGeneInfo(String geneName) {
        // Call Ensembl REST API
        String url = "https://rest.ensembl.org/lookup/symbol/homo_sapiens/" + geneName;
        // Parse JSON response and return GeneInfo object
    }
}
```

**File: `src/main/java/Entity/GeneInfo.java`** (NEW)
```java
public class GeneInfo {
    private String symbol;
    private String ensemblId;
    private String description;
    private String chromosome;
    private String location;
    private String biotype;
    // Getters and setters
}
```

**File: `src/main/java/Utils/ExternalAPIClient.java`** (NEW)
```java
public class ExternalAPIClient {
    public static GeneInfo fetchEnsemblData(String geneName) {}
    public static String generateNCBILink(String geneName) {}
    public static String generateGeneCardsLink(String geneName) {}
    // ... other link generators
}
```

#### **1.3 Dependencies to Add**
```xml
<!-- HTTP Client for API calls -->
<dependency>
    <groupId>org.apache.httpcomponents</groupId>
    <artifactId>httpclient</artifactId>
    <version>4.5.14</version>
</dependency>
```

#### **1.4 Database/Cache**
- Cache gene information in JSON files to reduce API calls
- Create `gene_cache/` directory in WEB-INF
- Implement TTL (time-to-live) of 7 days for cached data

---

## 2. Transcription Factor Activity Inference

### **Feature Description**
Calculate and visualize transcription factor (TF) activity scores across cell types using regulon-based methods.

### **Implementation Plan**

#### **2.1 Python Backend (pySCENIC Integration)**

**File: `src/main/webapp/WEB-INF/scenic_analysis.py`** (NEW)
```python
import scanpy as sc
import pandas as pd
import numpy as np
from pyscenic.rnkdb import FeatherRankingDatabase
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

def run_scenic_analysis(h5ad_path, output_path):
    """
    Run SCENIC TF activity inference
    """
    # Load data
    adata = sc.read_h5ad(h5ad_path)

    # Load ranking databases (download from SCENIC resources)
    db_fnames = ['hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather',
                 'hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather']

    dbs = [FeatherRankingDatabase(fname) for fname in db_fnames]

    # Step 1: GRN inference (co-expression modules)
    adjacencies = grnboost2(adata, tf_names=TF_NAMES)

    # Step 2: Regulon prediction (prune modules for TF-binding motifs)
    modules = modules_from_adjacencies(adjacencies)
    regulons = prune2df(dbs, modules, motif_annotations)

    # Step 3: AUCell scoring (calculate regulon activity)
    auc_mtx = aucell(adata, regulons)

    # Save results
    auc_mtx.to_csv(f'{output_path}/tf_activity_scores.csv')

    # Add to adata
    adata.obsm['X_scenic'] = auc_mtx
    adata.write(output_path + '/adata_with_scenic.h5ad')

    return auc_mtx

def get_top_regulons_per_celltype(adata, cell_type_key='cell_type', top_n=10):
    """
    Get top active regulons for each cell type
    """
    results = {}
    for ct in adata.obs[cell_type_key].unique():
        subset = adata[adata.obs[cell_type_key] == ct]
        mean_activity = subset.obsm['X_scenic'].mean(axis=0)
        top_tfs = mean_activity.nlargest(top_n)
        results[ct] = top_tfs.to_dict()

    return results
```

**File: `src/main/webapp/WEB-INF/tf_activity_dash.py`** (NEW)
```python
import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import plotly.graph_objects as go
import scanpy as sc
import pandas as pd

app = dash.Dash(__name__, url_base_pathname='/tf-activity/')

@app.callback(
    Output('tf-heatmap', 'figure'),
    Input('dataset-dropdown', 'value')
)
def update_heatmap(dataset_id):
    # Load SCENIC results
    adata = sc.read_h5ad(f'/data/{dataset_id}/adata_with_scenic.h5ad')

    # Compute mean TF activity per cell type
    tf_activity = pd.DataFrame(adata.obsm['X_scenic'],
                                columns=adata.uns['regulons'])
    tf_activity['cell_type'] = adata.obs['cell_type'].values

    heatmap_data = tf_activity.groupby('cell_type').mean()

    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=heatmap_data.values,
        x=heatmap_data.columns,
        y=heatmap_data.index,
        colorscale='RdBu_r'
    ))

    fig.update_layout(
        title='TF Activity Across Cell Types',
        xaxis_title='Transcription Factor',
        yaxis_title='Cell Type'
    )

    return fig

if __name__ == '__main__':
    app.run_server(port=8051)
```

#### **2.2 Java Servlet**

**File: `src/main/java/Servlet/TFActivityServlet.java`** (NEW)
```java
@WebServlet("/tf-activity")
public class TFActivityServlet extends HttpServlet {
    protected void doGet(HttpServletRequest request, HttpServletResponse response) {
        String datasetId = request.getParameter("dataset");

        // Check if SCENIC analysis exists
        String scenicPath = "/data/" + datasetId + "/scenic_results/";
        File scenicFile = new File(scenicPath + "tf_activity_scores.csv");

        if (!scenicFile.exists()) {
            // Trigger SCENIC analysis
            runScenicAnalysis(datasetId);
        }

        // Load TF activity data
        Map<String, Object> tfData = loadTFActivity(scenicPath);

        request.setAttribute("tfData", tfData);
        request.setAttribute("datasetId", datasetId);
        request.getRequestDispatcher("tf-activity.jsp").forward(request, response);
    }

    private void runScenicAnalysis(String datasetId) {
        try {
            ProcessBuilder pb = new ProcessBuilder(
                "python3",
                "WEB-INF/scenic_analysis.py",
                "--dataset", datasetId
            );
            Process p = pb.start();
            p.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
```

#### **2.3 Frontend**

**File: `src/main/webapp/tf-activity.jsp`** (NEW)
```jsp
<div class="tf-activity-container">
    <h1>Transcription Factor Activity - ${datasetId}</h1>

    <div class="tf-visualization">
        <iframe src="http://localhost:8051/tf-activity/?dataset=${datasetId}"
                width="100%" height="800px"></iframe>
    </div>

    <div class="tf-table">
        <h2>Top TFs per Cell Type</h2>
        <table id="tf-table">
            <thead>
                <tr>
                    <th>Cell Type</th>
                    <th>Top TFs</th>
                    <th>Activity Score</th>
                </tr>
            </thead>
            <tbody>
                <!-- Populated via AJAX -->
            </tbody>
        </table>
    </div>

    <div class="tf-network">
        <h2>TF Regulatory Network</h2>
        <div id="network-visualization"></div>
    </div>
</div>
```

#### **2.4 Dependencies**
```xml
<!-- Python packages -->
pyscenic==0.12.1
arboreto==0.1.6
loompy==3.0.7
```

#### **2.5 Data Requirements**
- Download SCENIC ranking databases (~10GB):
  - hg38 motif rankings from https://resources.aertslab.org/cistarget/
  - mm10 motif rankings for mouse data
- Store in `/data/scenic_resources/`

---

## 3. Pathway Enrichment Analysis

### **Feature Description**
Perform GO, KEGG, Reactome enrichment analysis on gene lists (e.g., DEGs, cluster markers).

### **Implementation Plan**

#### **3.1 Backend API Integration**

**File: `src/main/java/Utils/EnrichmentAPI.java`** (NEW)
```java
public class EnrichmentAPI {

    // Use Enrichr API
    public static JSONObject runEnrichrAnalysis(List<String> genes, String library) {
        String enrichrURL = "https://maayanlab.cloud/Enrichr/addList";

        // Submit gene list
        String listId = submitGeneList(genes);

        // Get enrichment results
        String resultsURL = "https://maayanlab.cloud/Enrichr/enrich?listId=" +
                           listId + "&backgroundType=" + library;

        // Parse and return results
        return fetchResults(resultsURL);
    }

    // Use g:Profiler API
    public static JSONObject runGProfilerAnalysis(List<String> genes, String organism) {
        String url = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/";

        JSONObject payload = new JSONObject();
        payload.put("organism", organism);
        payload.put("query", genes);
        payload.put("sources", Arrays.asList("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"));

        // POST request and return results
        return postRequest(url, payload);
    }
}
```

**File: `src/main/java/Servlet/EnrichmentServlet.java`** (NEW)
```java
@WebServlet("/enrichment")
public class EnrichmentServlet extends HttpServlet {
    protected void doPost(HttpServletRequest request, HttpServletResponse response) {
        // Get gene list from request
        String genesParam = request.getParameter("genes");
        String[] geneArray = genesParam.split(",");
        List<String> genes = Arrays.asList(geneArray);

        String analysisType = request.getParameter("type"); // "go", "kegg", "reactome"
        String organism = request.getParameter("organism"); // "hsapiens", "mmusculus"

        // Run enrichment
        JSONObject results = EnrichmentAPI.runGProfilerAnalysis(genes, organism);

        // Return JSON
        response.setContentType("application/json");
        response.getWriter().write(results.toString());
    }
}
```

#### **3.2 Python Alternative (for local analysis)**

**File: `src/main/webapp/WEB-INF/enrichment_analysis.py`** (NEW)
```python
import gseapy as gp
import pandas as pd

def run_enrichment(genes, organism='human', library='GO_Biological_Process_2021'):
    """
    Run enrichment analysis using gseapy
    """
    enr = gp.enrichr(
        gene_list=genes,
        gene_sets=library,
        organism=organism,
        outdir='enrichment_results'
    )

    return enr.results

def run_gsea(expression_matrix, gene_sets, phenotype_labels):
    """
    Run Gene Set Enrichment Analysis
    """
    gsea_results = gp.gsea(
        data=expression_matrix,
        gene_sets=gene_sets,
        cls=phenotype_labels,
        outdir='gsea_results'
    )

    return gsea_results

def plot_enrichment_dotplot(enrichment_df, top_n=20):
    """
    Create publication-ready dot plot
    """
    import matplotlib.pyplot as plt

    top_terms = enrichment_df.nsmallest(top_n, 'Adjusted P-value')

    fig, ax = plt.subplots(figsize=(10, 8))
    scatter = ax.scatter(
        top_terms['Combined Score'],
        range(len(top_terms)),
        s=top_terms['Odds Ratio'] * 10,
        c=-np.log10(top_terms['Adjusted P-value']),
        cmap='viridis'
    )

    ax.set_yticks(range(len(top_terms)))
    ax.set_yticklabels(top_terms['Term'])
    ax.set_xlabel('Combined Score')
    plt.colorbar(scatter, label='-log10(Adjusted P-value)')

    plt.tight_layout()
    plt.savefig('enrichment_dotplot.png', dpi=300)
```

#### **3.3 Frontend**

**File: `src/main/webapp/enrichment.jsp`** (NEW)
```jsp
<div class="enrichment-container">
    <h1>Pathway Enrichment Analysis</h1>

    <div class="enrichment-input">
        <h2>Input Gene List</h2>
        <textarea id="gene-input" rows="10" cols="50"
                  placeholder="Enter genes (one per line or comma-separated)">
        </textarea>

        <select id="organism-select">
            <option value="hsapiens">Human</option>
            <option value="mmusculus">Mouse</option>
        </select>

        <div class="library-selection">
            <h3>Select Libraries</h3>
            <label><input type="checkbox" value="GO_Biological_Process_2021" checked> GO Biological Process</label>
            <label><input type="checkbox" value="GO_Molecular_Function_2021"> GO Molecular Function</label>
            <label><input type="checkbox" value="GO_Cellular_Component_2021"> GO Cellular Component</label>
            <label><input type="checkbox" value="KEGG_2021_Human" checked> KEGG Pathways</label>
            <label><input type="checkbox" value="Reactome_2022"> Reactome</label>
            <label><input type="checkbox" value="WikiPathways_2021_Human"> WikiPathways</label>
            <label><input type="checkbox" value="MSigDB_Hallmark_2020"> MSigDB Hallmark</label>
        </div>

        <button id="run-enrichment">Run Enrichment</button>
    </div>

    <div class="enrichment-results" id="results-container">
        <h2>Results</h2>

        <div class="results-tabs">
            <button class="tab active" data-tab="table">Table</button>
            <button class="tab" data-tab="dotplot">Dot Plot</button>
            <button class="tab" data-tab="barplot">Bar Plot</button>
            <button class="tab" data-tab="network">Network</button>
        </div>

        <div id="table-view" class="tab-content active">
            <table id="enrichment-table">
                <thead>
                    <tr>
                        <th>Term</th>
                        <th>Database</th>
                        <th>P-value</th>
                        <th>Adjusted P-value</th>
                        <th>Odds Ratio</th>
                        <th>Combined Score</th>
                        <th>Genes</th>
                    </tr>
                </thead>
                <tbody></tbody>
            </table>
        </div>

        <div id="dotplot-view" class="tab-content">
            <div id="dotplot"></div>
        </div>

        <div id="barplot-view" class="tab-content">
            <div id="barplot"></div>
        </div>

        <div id="network-view" class="tab-content">
            <div id="enrichment-network"></div>
        </div>
    </div>

    <div class="export-section">
        <button id="export-csv">Export as CSV</button>
        <button id="export-png">Export Plot as PNG</button>
    </div>
</div>

<script>
document.getElementById('run-enrichment').addEventListener('click', async function() {
    const genes = document.getElementById('gene-input').value
        .split(/[,\n]/)
        .map(g => g.trim())
        .filter(g => g.length > 0);

    const organism = document.getElementById('organism-select').value;

    const libraries = Array.from(document.querySelectorAll('.library-selection input:checked'))
        .map(cb => cb.value);

    // Show loading
    document.getElementById('results-container').innerHTML = '<div class="loading">Running enrichment analysis...</div>';

    // Run analysis for each library
    for (const library of libraries) {
        const response = await fetch('/enrichment', {
            method: 'POST',
            headers: {'Content-Type': 'application/x-www-form-urlencoded'},
            body: `genes=${genes.join(',')}&organism=${organism}&library=${library}`
        });

        const results = await response.json();
        displayResults(results, library);
    }
});

function displayResults(results, library) {
    // Populate table
    const tbody = document.querySelector('#enrichment-table tbody');

    results.forEach(term => {
        const row = tbody.insertRow();
        row.innerHTML = `
            <td><a href="${term.link}" target="_blank">${term.term}</a></td>
            <td>${library}</td>
            <td>${term.pvalue.toExponential(2)}</td>
            <td>${term.adjusted_pvalue.toExponential(2)}</td>
            <td>${term.odds_ratio.toFixed(2)}</td>
            <td>${term.combined_score.toFixed(2)}</td>
            <td>${term.genes.join(', ')}</td>
        `;
    });

    // Create visualizations
    createDotPlot(results);
    createBarPlot(results);
    createEnrichmentNetwork(results);
}
</script>
```

**File: `src/main/webapp/CSS/enrichment.css`** (NEW)
```css
.enrichment-container {
    max-width: 1400px;
    margin: 0 auto;
    padding: 20px;
}

.enrichment-input {
    background: white;
    padding: 20px;
    border-radius: 8px;
    margin-bottom: 20px;
}

#gene-input {
    width: 100%;
    font-family: 'JetBrains Mono', monospace;
    padding: 10px;
    border: 1px solid #ddd;
    border-radius: 4px;
}

.library-selection {
    margin: 20px 0;
}

.library-selection label {
    display: block;
    margin: 8px 0;
}

#enrichment-table {
    width: 100%;
    border-collapse: collapse;
}

#enrichment-table th,
#enrichment-table td {
    padding: 12px;
    text-align: left;
    border-bottom: 1px solid #ddd;
}

#enrichment-table th {
    background: #f5f5f5;
    font-weight: 600;
}

.tab-content {
    display: none;
}

.tab-content.active {
    display: block;
}
```

#### **3.4 Dependencies**
```xml
<!-- Python -->
gseapy==1.0.4
pandas==1.5.3
matplotlib==3.7.1
```

---

## 4. Cell-Cell Communication Enhancement

### **Feature Description**
Enhanced visualization and analysis of cell-cell communication beyond existing CellPhoneDB integration.

### **Implementation Plan**

#### **4.1 Python Analysis**

**File: `src/main/webapp/WEB-INF/cellchat_analysis.py`** (NEW)
```python
import scanpy as sc
import pandas as pd
import cellphonedb
from cellphonedb.utils import db_utils
from cellphonedb.core.methods import cpdb_statistical_analysis_method

def run_cellphonedb_analysis(adata, cell_type_key='cell_type', output_path='cellphonedb_out'):
    """
    Run CellPhoneDB analysis
    """
    # Prepare counts and metadata
    counts = pd.DataFrame(
        adata.raw.X.toarray().T,
        index=adata.var_names,
        columns=adata.obs_names
    )

    meta = adata.obs[[cell_type_key]].copy()
    meta.columns = ['cell_type']

    # Run CellPhoneDB
    cpdb_results = cpdb_statistical_analysis_method.call(
        cpdb_file_path=db_utils.get_database_path(),
        meta_file_path=meta,
        counts_file_path=counts,
        output_path=output_path,
        iterations=1000,
        threshold=0.1
    )

    return cpdb_results

def create_interaction_network(cpdb_results, significance_threshold=0.05):
    """
    Create network graph of cell-cell interactions
    """
    import networkx as nx
    import plotly.graph_objects as go

    # Filter significant interactions
    sig_interactions = cpdb_results[cpdb_results['pvalue'] < significance_threshold]

    # Build graph
    G = nx.Graph()

    for _, row in sig_interactions.iterrows():
        source_ct = row['cell_type_a']
        target_ct = row['cell_type_b']
        interaction = row['interacting_pair']

        G.add_edge(source_ct, target_ct,
                  interaction=interaction,
                  weight=row['mean'])

    # Create plotly visualization
    pos = nx.spring_layout(G)

    edge_trace = []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_trace.append(
            go.Scatter(
                x=[x0, x1, None],
                y=[y0, y1, None],
                mode='lines',
                line=dict(width=G[edge[0]][edge[1]]['weight'], color='#888')
            )
        )

    node_trace = go.Scatter(
        x=[pos[node][0] for node in G.nodes()],
        y=[pos[node][1] for node in G.nodes()],
        mode='markers+text',
        text=list(G.nodes()),
        textposition='top center',
        marker=dict(size=20, color='lightblue')
    )

    fig = go.Figure(data=edge_trace + [node_trace])
    fig.update_layout(showlegend=False)

    return fig
```

#### **4.2 Interactive Visualization**

**File: `src/main/webapp/WEB-INF/cellcomm_dash.py`** (NEW)
```python
import dash
from dash import dcc, html, dash_table
from dash.dependencies import Input, Output
import plotly.graph_objects as go
import pandas as pd

app = dash.Dash(__name__, url_base_pathname='/cell-communication/')

@app.callback(
    Output('circos-plot', 'figure'),
    Input('dataset-selector', 'value')
)
def update_circos(dataset_id):
    # Load CellPhoneDB results
    means = pd.read_csv(f'/data/{dataset_id}/cellphonedb_out/means.txt', sep='\t')
    pvalues = pd.read_csv(f'/data/{dataset_id}/cellphonedb_out/pvalues.txt', sep='\t')

    # Create chord diagram
    fig = create_chord_diagram(means, pvalues)
    return fig

def create_chord_diagram(means, pvalues, threshold=0.05):
    """
    Create circular chord diagram for cell-cell communication
    """
    # Filter significant interactions
    sig_mask = pvalues.iloc[:, 11:] < threshold
    sig_means = means.iloc[:, 11:].where(sig_mask, 0)

    # Aggregate by cell type pairs
    cell_types = sig_means.columns
    interaction_matrix = pd.DataFrame(0, index=cell_types, columns=cell_types)

    for col in sig_means.columns:
        parts = col.split('|')
        if len(parts) == 2:
            ct1, ct2 = parts
            interaction_matrix.loc[ct1, ct2] = sig_means[col].sum()

    # Create Sankey diagram
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            label=list(cell_types),
            color='blue'
        ),
        link=dict(
            source=[],  # indices
            target=[],  # indices
            value=[]    # interaction strength
        )
    )])

    return fig
```

---

## 5. Trajectory & Pseudotime Analysis

### **Feature Description**
Infer developmental trajectories and pseudotime ordering of cells.

### **Implementation Plan**

#### **5.1 Python Backend**

**File: `src/main/webapp/WEB-INF/trajectory_analysis.py`** (NEW)
```python
import scanpy as sc
import scvelo as scv
import pandas as pd

def run_rna_velocity(adata, loom_path):
    """
    Run RNA velocity analysis using scVelo
    """
    # Load spliced/unspliced counts from loom
    ldata = scv.read(loom_path)

    # Merge with main adata
    adata = scv.utils.merge(adata, ldata)

    # Preprocess
    scv.pp.filter_and_normalize(adata)
    scv.pp.moments(adata)

    # Compute velocity
    scv.tl.velocity(adata, mode='stochastic')
    scv.tl.velocity_graph(adata)

    # Visualize
    scv.pl.velocity_embedding_stream(adata, basis='umap', save='velocity_stream.png')

    return adata

def run_paga_trajectory(adata, root_cell=None):
    """
    Run PAGA trajectory inference
    """
    # Compute neighbors if not already done
    if 'neighbors' not in adata.uns:
        sc.pp.neighbors(adata)

    # Run PAGA
    sc.tl.paga(adata, groups='leiden')

    # Compute trajectory
    if root_cell:
        adata.uns['iroot'] = root_cell
        sc.tl.dpt(adata)

    # Visualize
    sc.pl.paga(adata, save='_paga.png')

    return adata

def run_diffusion_pseudotime(adata, root_cell_type):
    """
    Compute diffusion pseudotime
    """
    # Set root
    root_cells = adata.obs[adata.obs['cell_type'] == root_cell_type].index
    adata.uns['iroot'] = np.flatnonzero(adata.obs_names == root_cells[0])[0]

    # Compute DPT
    sc.tl.diffmap(adata)
    sc.tl.dpt(adata)

    # Plot
    sc.pl.dpt_groups_pseudotime(adata, save='_dpt.png')

    return adata
```

**File: `src/main/webapp/WEB-INF/trajectory_dash.py`** (NEW)
```python
import dash
from dash import dcc, html
import plotly.graph_objects as go
import scanpy as sc

app = dash.Dash(__name__, url_base_pathname='/trajectory/')

@app.callback(
    Output('trajectory-plot', 'figure'),
    Input('dataset-selector', 'value')
)
def update_trajectory(dataset_id):
    adata = sc.read_h5ad(f'/data/{dataset_id}/trajectory_adata.h5ad')

    # Create 3D trajectory plot
    fig = go.Figure(data=[go.Scatter3d(
        x=adata.obsm['X_umap'][:, 0],
        y=adata.obsm['X_umap'][:, 1],
        z=adata.obs['dpt_pseudotime'],
        mode='markers',
        marker=dict(
            size=3,
            color=adata.obs['dpt_pseudotime'],
            colorscale='Viridis',
            showscale=True
        ),
        text=adata.obs['cell_type']
    )])

    fig.update_layout(
        title='Pseudotime Trajectory',
        scene=dict(
            xaxis_title='UMAP1',
            yaxis_title='UMAP2',
            zaxis_title='Pseudotime'
        )
    )

    return fig
```

---

## 6. Gene Regulatory Networks

### **Feature Description**
Visualize gene regulatory networks and co-expression modules.

### **Implementation Plan**

#### **6.1 Python Backend**

**File: `src/main/webapp/WEB-INF/grn_analysis.py`** (NEW)
```python
import scanpy as sc
import pandas as pd
import networkx as nx
from arboreto.algo import grnboost2
from arboreto.utils import load_tf_names

def infer_grn(adata, tf_names=None):
    """
    Infer gene regulatory network using GRNBoost2
    """
    # Prepare expression matrix
    expr_matrix = pd.DataFrame(
        adata.X.toarray(),
        index=adata.obs_names,
        columns=adata.var_names
    )

    # Load TF names if not provided
    if tf_names is None:
        tf_names = load_tf_names('resources/hs_hgnc_tfs.txt')

    # Run GRNBoost2
    network = grnboost2(expr_matrix, tf_names=tf_names)

    return network

def create_grn_visualization(network, top_n=100):
    """
    Create interactive network visualization
    """
    import plotly.graph_objects as go

    # Filter top edges
    top_edges = network.nlargest(top_n, 'importance')

    # Build NetworkX graph
    G = nx.from_pandas_edgelist(
        top_edges,
        source='TF',
        target='target',
        edge_attr='importance'
    )

    # Calculate layout
    pos = nx.spring_layout(G, k=0.5, iterations=50)

    # Create edges
    edge_trace = []
    for edge in G.edges(data=True):
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]

        edge_trace.append(
            go.Scatter(
                x=[x0, x1, None],
                y=[y0, y1, None],
                mode='lines',
                line=dict(width=edge[2]['importance']*2, color='#888'),
                hoverinfo='none'
            )
        )

    # Create nodes
    node_x = [pos[node][0] for node in G.nodes()]
    node_y = [pos[node][1] for node in G.nodes()]

    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        mode='markers+text',
        text=list(G.nodes()),
        textposition='top center',
        marker=dict(
            size=10,
            color='lightblue',
            line=dict(width=2, color='darkblue')
        )
    )

    fig = go.Figure(data=edge_trace + [node_trace])
    fig.update_layout(
        title='Gene Regulatory Network',
        showlegend=False,
        hovermode='closest'
    )

    return fig
```

---

## 7. Interactive Visualization Suite

### **Feature Description**
Comprehensive interactive visualization tools beyond basic UMAP.

### **Implementation Plan**

#### **7.1 Features to Implement**

1. **3D UMAP/tSNE viewer**
2. **Split violin plots**
3. **Dot plots (gene expression × cell type)**
4. **Stacked violin plots**
5. **Feature plots (gene overlay on UMAP)**
6. **Heatmaps with hierarchical clustering**
7. **Sankey diagrams for cell type proportions**
8. **Interactive correlation plots**

#### **7.2 Dash Application**

**File: `src/main/webapp/WEB-INF/visualization_suite.py`** (NEW)
```python
import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State
import plotly.graph_objects as go
import plotly.express as px
import scanpy as sc
import pandas as pd

app = dash.Dash(__name__, url_base_pathname='/viz/')

app.layout = html.Div([
    html.H1('Interactive Visualization Suite'),

    dcc.Tabs(id='viz-tabs', value='umap', children=[
        dcc.Tab(label='3D UMAP', value='umap'),
        dcc.Tab(label='Violin Plots', value='violin'),
        dcc.Tab(label='Dot Plot', value='dotplot'),
        dcc.Tab(label='Heatmap', value='heatmap'),
        dcc.Tab(label='Feature Plot', value='feature'),
    ]),

    html.Div(id='viz-content')
])

@app.callback(
    Output('viz-content', 'children'),
    Input('viz-tabs', 'value'),
    State('dataset-id', 'data')
)
def render_visualization(tab, dataset_id):
    adata = sc.read_h5ad(f'/data/{dataset_id}/adata.h5ad')

    if tab == 'umap':
        return create_3d_umap(adata)
    elif tab == 'violin':
        return create_violin_plot(adata)
    elif tab == 'dotplot':
        return create_dot_plot(adata)
    elif tab == 'heatmap':
        return create_heatmap(adata)
    elif tab == 'feature':
        return create_feature_plot(adata)

def create_3d_umap(adata):
    """
    Create 3D UMAP with rotation and zoom
    """
    # Compute 3D UMAP if not exists
    if adata.obsm['X_umap'].shape[1] == 2:
        sc.tl.umap(adata, n_components=3)

    fig = go.Figure(data=[go.Scatter3d(
        x=adata.obsm['X_umap'][:, 0],
        y=adata.obsm['X_umap'][:, 1],
        z=adata.obsm['X_umap'][:, 2],
        mode='markers',
        marker=dict(
            size=2,
            color=adata.obs['cell_type'].cat.codes,
            colorscale='Viridis',
            showscale=True
        ),
        text=adata.obs['cell_type']
    )])

    fig.update_layout(
        scene=dict(
            xaxis_title='UMAP1',
            yaxis_title='UMAP2',
            zaxis_title='UMAP3'
        ),
        height=800
    )

    return dcc.Graph(figure=fig)

def create_dot_plot(adata, genes=None, groupby='cell_type'):
    """
    Create dot plot showing gene expression across cell types
    """
    if genes is None:
        # Use top marker genes
        genes = get_top_markers(adata, groupby, n_genes=5)

    # Calculate mean expression and fraction expressing
    dot_data = []

    for ct in adata.obs[groupby].unique():
        subset = adata[adata.obs[groupby] == ct]

        for gene in genes:
            if gene in adata.var_names:
                expr = subset[:, gene].X.toarray().flatten()
                mean_expr = expr.mean()
                pct_expr = (expr > 0).sum() / len(expr) * 100

                dot_data.append({
                    'cell_type': ct,
                    'gene': gene,
                    'mean_expression': mean_expr,
                    'pct_expressed': pct_expr
                })

    df = pd.DataFrame(dot_data)

    fig = go.Figure()

    for gene in genes:
        gene_data = df[df['gene'] == gene]

        fig.add_trace(go.Scatter(
            x=gene_data['cell_type'],
            y=[gene] * len(gene_data),
            mode='markers',
            marker=dict(
                size=gene_data['pct_expressed'],
                color=gene_data['mean_expression'],
                colorscale='Reds',
                showscale=True,
                sizemode='diameter',
                sizemin=4
            ),
            name=gene
        ))

    fig.update_layout(
        title='Dot Plot: Gene Expression × Cell Type',
        xaxis_title='Cell Type',
        yaxis_title='Gene',
        height=600
    )

    return dcc.Graph(figure=fig)
```

---

## 8. Comparative Analysis Tools

### **Feature Description**
Compare datasets, conditions, species, or cell types.

### **Implementation Plan**

#### **8.1 Features**

1. **Disease vs Healthy comparison**
2. **Human vs Mouse ortholog mapping**
3. **Cross-dataset integration quality**
4. **Differential abundance testing**

#### **8.2 Python Backend**

**File: `src/main/webapp/WEB-INF/comparative_analysis.py`** (NEW)
```python
import scanpy as sc
import pandas as pd
from scipy import stats

def compare_disease_vs_healthy(adata, condition_key='disease_status'):
    """
    Compare gene expression between disease and healthy
    """
    # Split data
    healthy = adata[adata.obs[condition_key] == 'healthy']
    disease = adata[adata.obs[condition_key] == 'disease']

    # Perform differential expression
    sc.tl.rank_genes_groups(adata, groupby=condition_key, method='wilcoxon')

    # Get results
    result = adata.uns['rank_genes_groups']
    degs = pd.DataFrame({
        'gene': result['names']['disease'],
        'logfoldchange': result['logfoldchanges']['disease'],
        'pval': result['pvals']['disease'],
        'pval_adj': result['pvals_adj']['disease']
    })

    return degs

def differential_abundance_testing(adata, condition_key, cell_type_key):
    """
    Test for cell type proportion changes between conditions
    """
    from statsmodels.stats.multitest import multipletests

    # Calculate proportions
    props = adata.obs.groupby([condition_key, cell_type_key]).size().unstack(fill_value=0)
    props = props.div(props.sum(axis=1), axis=0)

    # Statistical test
    results = []

    for cell_type in props.columns:
        healthy_prop = props.loc['healthy', cell_type]
        disease_prop = props.loc['disease', cell_type]

        # Chi-square test
        contingency = [[healthy_counts[cell_type], disease_counts[cell_type]],
                      [healthy_total - healthy_counts[cell_type],
                       disease_total - disease_counts[cell_type]]]

        chi2, pval = stats.chi2_contingency(contingency)[:2]

        results.append({
            'cell_type': cell_type,
            'healthy_pct': healthy_prop * 100,
            'disease_pct': disease_prop * 100,
            'fold_change': disease_prop / healthy_prop,
            'pvalue': pval
        })

    df = pd.DataFrame(results)
    df['padj'] = multipletests(df['pvalue'], method='fdr_bh')[1]

    return df
```

---

## 9. Advanced Search & Filtering

### **Feature Description**
Enhanced search capabilities beyond gene names.

### **Implementation Plan**

#### **9.1 Features**

1. **Cell type marker search** - Find datasets containing specific cell types
2. **Disease/condition filtering** - Filter by pathology
3. **Tissue/anatomical location** - Filter by skin region
4. **Multi-gene co-expression** - Find cells expressing gene combinations
5. **Metadata search** - Search by author, publication, technique

#### **9.2 Backend**

**File: `src/main/java/Servlet/AdvancedSearchServlet.java`** (NEW)
```java
@WebServlet("/advanced-search")
public class AdvancedSearchServlet extends HttpServlet {
    protected void doPost(HttpServletRequest request, HttpServletResponse response) {
        // Parse search parameters
        String cellType = request.getParameter("cellType");
        String disease = request.getParameter("disease");
        String tissue = request.getParameter("tissue");
        String[] genes = request.getParameter("genes").split(",");
        String species = request.getParameter("species");
        String technique = request.getParameter("technique");

        // Build search query
        List<Dataset> results = performAdvancedSearch(
            cellType, disease, tissue, genes, species, technique
        );

        // Return JSON results
        Gson gson = new Gson();
        String json = gson.toJson(results);

        response.setContentType("application/json");
        response.getWriter().write(json);
    }

    private List<Dataset> performAdvancedSearch(
        String cellType, String disease, String tissue,
        String[] genes, String species, String technique
    ) {
        // Load metadata
        List<Dataset> allDatasets = DatasetLoader.loadFromExcel("WEB-INF/AllData.xlsx");

        // Filter by criteria
        return allDatasets.stream()
            .filter(d -> cellType == null || d.getCellTypes().contains(cellType))
            .filter(d -> disease == null || d.getDisease().equals(disease))
            .filter(d -> tissue == null || d.getTissue().equals(tissue))
            .filter(d -> species == null || d.getSpecies().equals(species))
            .filter(d -> hasGenes(d, genes))
            .collect(Collectors.toList());
    }
}
```

---

## 10. Cell Type Annotation & Prediction

### **Feature Description**
Automated cell type annotation using reference atlases.

### **Implementation Plan**

#### **10.1 Python Backend**

**File: `src/main/webapp/WEB-INF/cell_annotation.py`** (NEW)
```python
import scanpy as sc
import celltypist
from celltypist import models

def auto_annotate_celltypist(adata, model='Immune_All_Low.pkl'):
    """
    Automatically annotate cell types using CellTypist
    """
    # Download model if needed
    models.download_models(model=model)

    # Load model
    model = models.Model.load(model=model)

    # Predict
    predictions = celltypist.annotate(adata, model=model)

    # Add to adata
    adata.obs['celltypist_prediction'] = predictions.predicted_labels
    adata.obs['celltypist_confidence'] = predictions.probability

    return adata

def annotate_with_singler(adata, reference_path):
    """
    Annotate using SingleR (via rpy2)
    """
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

    pandas2ri.activate()

    # Load R packages
    ro.r('library(SingleR)')
    ro.r('library(celldex)')

    # Load reference
    ro.r(f'ref <- readRDS("{reference_path}")')

    # Convert adata to R
    expr_matrix = pandas2ri.py2rpy(pd.DataFrame(adata.X.toarray().T,
                                                 index=adata.var_names,
                                                 columns=adata.obs_names))

    # Run SingleR
    ro.globalenv['test_expr'] = expr_matrix
    predictions = ro.r('SingleR(test=test_expr, ref=ref, labels=ref$label.main)')

    # Add to adata
    adata.obs['singler_prediction'] = predictions.rx2('labels')

    return adata
```

---

## 11. REST API & Data Access

### **Feature Description**
Programmatic access to data and analyses.

### **Implementation Plan**

#### **11.1 API Endpoints**

**File: `src/main/java/API/DatasetAPI.java`** (NEW)
```java
@WebServlet("/api/v1/datasets")
public class DatasetAPI extends HttpServlet {

    @Override
    protected void doGet(HttpServletRequest request, HttpServletResponse response) {
        String action = request.getParameter("action");

        switch (action) {
            case "list":
                listDatasets(request, response);
                break;
            case "get":
                getDataset(request, response);
                break;
            case "search":
                searchDatasets(request, response);
                break;
            default:
                sendError(response, "Invalid action");
        }
    }

    private void listDatasets(HttpServletRequest request, HttpServletResponse response) {
        int page = Integer.parseInt(request.getParameter("page") != null ?
                                   request.getParameter("page") : "1");
        int pageSize = Integer.parseInt(request.getParameter("pageSize") != null ?
                                       request.getParameter("pageSize") : "50");

        List<Dataset> datasets = DatasetLoader.loadFromExcel("WEB-INF/AllData.xlsx");

        // Pagination
        int start = (page - 1) * pageSize;
        int end = Math.min(start + pageSize, datasets.size());
        List<Dataset> pageData = datasets.subList(start, end);

        // Build response
        Map<String, Object> responseData = new HashMap<>();
        responseData.put("datasets", pageData);
        responseData.put("total", datasets.size());
        responseData.put("page", page);
        responseData.put("pageSize", pageSize);

        sendJSON(response, responseData);
    }
}

@WebServlet("/api/v1/genes")
public class GeneAPI extends HttpServlet {
    // Endpoints for gene expression queries
}

@WebServlet("/api/v1/cells")
public class CellAPI extends HttpServlet {
    // Endpoints for cell metadata and expression
}
```

#### **11.2 API Documentation**

**File: `src/main/webapp/api-docs.jsp`** (NEW)
```jsp
<div class="api-docs">
    <h1>scSAID REST API Documentation</h1>

    <section>
        <h2>Authentication</h2>
        <p>API key required for all requests: <code>X-API-Key: your_key_here</code></p>
    </section>

    <section>
        <h2>Endpoints</h2>

        <div class="endpoint">
            <h3>GET /api/v1/datasets</h3>
            <p>List all datasets</p>

            <h4>Parameters:</h4>
            <ul>
                <li><code>page</code> (int): Page number (default: 1)</li>
                <li><code>pageSize</code> (int): Results per page (default: 50)</li>
                <li><code>species</code> (string): Filter by species</li>
                <li><code>disease</code> (string): Filter by disease</li>
            </ul>

            <h4>Example:</h4>
            <pre><code>curl -H "X-API-Key: your_key" \
  "https://scsaid.org/api/v1/datasets?page=1&pageSize=10"
</code></pre>

            <h4>Response:</h4>
            <pre><code>{
  "datasets": [
    {
      "SAID": "SAID001",
      "GSE": "GSE123456",
      "title": "Single-cell analysis of human skin",
      "species": "Homo sapiens",
      "disease": "healthy"
    }
  ],
  "total": 600,
  "page": 1,
  "pageSize": 10
}
</code></pre>
        </div>

        <!-- More endpoints -->
    </section>
</div>
```

---

## 12. User Workspace & Saved Analyses

### **Feature Description**
Allow users to save analyses, bookmark datasets, and share results.

### **Implementation Plan**

#### **12.1 Database Schema**

**File: `src/main/resources/schema.sql`** (NEW)
```sql
CREATE TABLE users (
    user_id INT PRIMARY KEY AUTO_INCREMENT,
    email VARCHAR(255) UNIQUE NOT NULL,
    password_hash VARCHAR(255) NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE saved_analyses (
    analysis_id INT PRIMARY KEY AUTO_INCREMENT,
    user_id INT,
    analysis_type VARCHAR(50),
    parameters TEXT,  -- JSON
    results_path VARCHAR(500),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (user_id) REFERENCES users(user_id)
);

CREATE TABLE bookmarks (
    bookmark_id INT PRIMARY KEY AUTO_INCREMENT,
    user_id INT,
    dataset_id VARCHAR(50),
    notes TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (user_id) REFERENCES users(user_id)
);

CREATE TABLE shared_analyses (
    share_id VARCHAR(36) PRIMARY KEY,  -- UUID
    analysis_id INT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    expires_at TIMESTAMP,
    FOREIGN KEY (analysis_id) REFERENCES saved_analyses(analysis_id)
);
```

#### **12.2 Backend**

**File: `src/main/java/Servlet/UserWorkspaceServlet.java`** (NEW)
```java
@WebServlet("/workspace")
public class UserWorkspaceServlet extends HttpServlet {
    protected void doGet(HttpServletRequest request, HttpServletResponse response) {
        HttpSession session = request.getSession();
        Integer userId = (Integer) session.getAttribute("userId");

        if (userId == null) {
            response.sendRedirect("login.jsp");
            return;
        }

        // Load user's saved analyses and bookmarks
        List<SavedAnalysis> analyses = getSavedAnalyses(userId);
        List<Bookmark> bookmarks = getBookmarks(userId);

        request.setAttribute("analyses", analyses);
        request.setAttribute("bookmarks", bookmarks);

        request.getRequestDispatcher("workspace.jsp").forward(request, response);
    }
}
```

---

## 13. Custom Dataset Upload & Analysis

### **Feature Description**
Allow researchers to upload their own datasets and analyze them using scSAID tools.

### **Implementation Plan**

#### **13.1 Upload Interface**

**File: `src/main/webapp/upload-dataset.jsp`** (NEW)
```jsp
<div class="upload-container">
    <h1>Upload Custom Dataset</h1>

    <form id="upload-form" enctype="multipart/form-data">
        <div class="form-group">
            <label>Dataset Format:</label>
            <select name="format">
                <option value="h5ad">H5AD (AnnData)</option>
                <option value="10x">10X Genomics (MEX)</option>
                <option value="loom">Loom</option>
                <option value="csv">CSV (genes × cells)</option>
            </select>
        </div>

        <div class="form-group">
            <label>Upload File:</label>
            <input type="file" name="dataset" accept=".h5ad,.h5,.loom,.csv,.mtx" />
        </div>

        <div class="form-group">
            <label>Metadata (optional):</label>
            <input type="file" name="metadata" accept=".csv,.tsv" />
        </div>

        <div class="form-group">
            <label>Dataset Name:</label>
            <input type="text" name="dataset_name" required />
        </div>

        <div class="form-group">
            <label>Species:</label>
            <select name="species">
                <option value="human">Human</option>
                <option value="mouse">Mouse</option>
            </select>
        </div>

        <button type="submit">Upload & Analyze</button>
    </form>

    <div id="upload-progress" style="display:none;">
        <div class="progress-bar">
            <div class="progress-fill"></div>
        </div>
        <p>Processing... This may take several minutes.</p>
    </div>
</div>
```

#### **13.2 Backend Processing**

**File: `src/main/java/Servlet/UploadDatasetServlet.java`** (NEW)
```java
@WebServlet("/upload-dataset")
@MultipartConfig(maxFileSize = 5368709120L) // 5GB max
public class UploadDatasetServlet extends HttpServlet {
    protected void doPost(HttpServletRequest request, HttpServletResponse response) {
        Part filePart = request.getPart("dataset");
        String format = request.getParameter("format");
        String datasetName = request.getParameter("dataset_name");

        // Save file
        String uploadId = UUID.randomUUID().toString();
        String uploadPath = "/data/user_uploads/" + uploadId + "/";
        Files.createDirectories(Paths.get(uploadPath));

        String filename = Paths.get(filePart.getSubmittedFileName()).getFileName().toString();
        filePart.write(uploadPath + filename);

        // Trigger processing pipeline
        ProcessBuilder pb = new ProcessBuilder(
            "python3",
            "WEB-INF/process_upload.py",
            "--input", uploadPath + filename,
            "--format", format,
            "--output", uploadPath
        );

        Process p = pb.start();

        // Return upload ID for tracking
        response.getWriter().write("{\"upload_id\": \"" + uploadId + "\"}");
    }
}
```

**File: `src/main/webapp/WEB-INF/process_upload.py`** (NEW)
```python
import scanpy as sc
import argparse

def process_uploaded_dataset(input_path, format, output_path):
    """
    Process uploaded dataset and run standard pipeline
    """
    # Load data
    if format == 'h5ad':
        adata = sc.read_h5ad(input_path)
    elif format == '10x':
        adata = sc.read_10x_mtx(input_path)
    elif format == 'loom':
        adata = sc.read_loom(input_path)
    elif format == 'csv':
        adata = sc.read_csv(input_path).T

    # QC
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

    # Normalization
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Feature selection
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)

    # Dimensionality reduction
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    # Clustering
    sc.tl.leiden(adata)

    # Marker genes
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

    # Save
    adata.write(output_path + '/processed_adata.h5ad')

    print(f"Processing complete. Output: {output_path}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--format', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    process_uploaded_dataset(args.input, args.format, args.output)
```

---

## 14. Publication-Ready Figure Generation

### **Feature Description**
Generate high-quality, publication-ready figures with customization options.

### **Implementation Plan**

#### **14.1 Python Backend**

**File: `src/main/webapp/WEB-INF/figure_generator.py`** (NEW)
```python
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

def generate_publication_umap(adata, color_by='cell_type', save_path='figures/'):
    """
    Create publication-quality UMAP
    """
    # Set publication style
    sns.set_style('white')
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.linewidth'] = 1.5

    fig, ax = plt.subplots(figsize=(8, 6))

    sc.pl.umap(
        adata,
        color=color_by,
        ax=ax,
        frameon=False,
        legend_fontsize=10,
        legend_fontweight='normal',
        show=False
    )

    plt.tight_layout()
    plt.savefig(save_path + 'umap.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(save_path + 'umap.png', dpi=300, bbox_inches='tight')

    return save_path + 'umap.pdf'

def generate_figure_panel(adata, plots=['umap', 'dotplot', 'violin'], genes=None):
    """
    Create multi-panel figure
    """
    fig = plt.figure(figsize=(18, 6))

    # UMAP
    ax1 = plt.subplot(1, 3, 1)
    sc.pl.umap(adata, color='cell_type', ax=ax1, show=False)
    ax1.set_title('A', fontweight='bold', loc='left', fontsize=16)

    # Dot plot
    ax2 = plt.subplot(1, 3, 2)
    sc.pl.dotplot(adata, genes, groupby='cell_type', ax=ax2, show=False)
    ax2.set_title('B', fontweight='bold', loc='left', fontsize=16)

    # Violin plot
    ax3 = plt.subplot(1, 3, 3)
    sc.pl.stacked_violin(adata, genes, groupby='cell_type', ax=ax3, show=False)
    ax3.set_title('C', fontweight='bold', loc='left', fontsize=16)

    plt.tight_layout()
    plt.savefig('figures/figure_panel.pdf', dpi=300)

    return 'figures/figure_panel.pdf'
```

---

## Implementation Priority & Timeline

### **Phase 1: Foundation (Months 1-2)**
1. Gene Information & External Links
2. REST API & Data Access
3. Advanced Search & Filtering

### **Phase 2: Core Analytics (Months 3-5)**
4. Pathway Enrichment Analysis
5. TF Activity Inference
6. Cell-Cell Communication Enhancement

### **Phase 3: Advanced Features (Months 6-8)**
7. Trajectory & Pseudotime Analysis
8. Gene Regulatory Networks
9. Interactive Visualization Suite

### **Phase 4: User Features (Months 9-10)**
10. Cell Type Annotation & Prediction
11. User Workspace & Saved Analyses
12. Custom Dataset Upload

### **Phase 5: Publication Tools (Months 11-12)**
13. Comparative Analysis Tools
14. Publication-Ready Figure Generation

---

## Technical Requirements Summary

### **Python Dependencies**
```
scanpy>=1.9.0
scvelo>=0.2.5
pyscenic>=0.12.1
gseapy>=1.0.4
celltypist>=1.3.0
cellphonedb>=3.0.0
dash>=2.8.0
plotly>=5.13.0
networkx>=2.8
arboreto>=0.1.6
```

### **Java Dependencies** (add to pom.xml)
```xml
<dependency>
    <groupId>org.apache.httpcomponents</groupId>
    <artifactId>httpclient</artifactId>
    <version>4.5.14</version>
</dependency>

<dependency>
    <groupId>mysql</groupId>
    <artifactId>mysql-connector-java</artifactId>
    <version>8.0.33</version>
</dependency>

<dependency>
    <groupId>com.auth0</groupId>
    <artifactId>java-jwt</artifactId>
    <version>4.3.0</version>
</dependency>
```

### **Database**
- MySQL 8.0+ or PostgreSQL 14+
- For user accounts, saved analyses, bookmarks

### **Storage Requirements**
- SCENIC databases: ~10GB
- Reference atlases: ~5GB per species
- User uploads: Scale as needed (recommend 1TB initial)

### **Compute Resources**
- RAM: 64GB minimum (128GB recommended for large integrations)
- CPU: 16+ cores for parallel processing
- GPU: Optional but recommended for deep learning features

---

## Estimated Development Effort

| Feature | Complexity | Est. Time | Priority |
|---------|-----------|-----------|----------|
| Gene Links | Low | 1 week | High |
| Enrichment Analysis | Medium | 2-3 weeks | High |
| TF Activity | High | 4-5 weeks | Medium |
| Cell Communication | Medium | 2-3 weeks | Medium |
| Trajectory Analysis | High | 4-5 weeks | Medium |
| GRN | High | 3-4 weeks | Low |
| Visualization Suite | Medium | 3-4 weeks | High |
| Comparative Analysis | Medium | 2-3 weeks | Medium |
| Advanced Search | Low | 1-2 weeks | High |
| Cell Annotation | Medium | 2 weeks | Medium |
| REST API | Medium | 2-3 weeks | High |
| User Workspace | High | 4-5 weeks | Low |
| Upload Feature | High | 3-4 weeks | Low |
| Figure Generator | Low | 1-2 weeks | Low |

**Total estimated development time: 8-12 months with 2-3 developers**

---

This comprehensive plan provides a roadmap to transform scSAID into a world-class single-cell analysis platform comparable to CellxGene, UCSC Cell Browser, and other leading tools in the field.
