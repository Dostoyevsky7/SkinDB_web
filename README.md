<p align="center">
  <img src="./src/main/webapp/images/scSAID_LOGO.png" alt="scSAID Logo" width="180"/>
</p>

<h1 align="center">scSAID</h1>

<p align="center">
  <strong>Single-Cell Skin & Appendages Integrated Database</strong>
</p>

<p align="center">
  <a href="https://skin-scsaid.com/">
    <img src="https://img.shields.io/badge/Live-skin--scsaid.com-e8927c?style=for-the-badge&logo=safari&logoColor=white" alt="Live Site"/>
  </a>
</p>

<p align="center">
  <a href="https://github.com/Dostoyevsky7/SkinDB_web/actions/workflows/deploy.yml">
    <img src="https://github.com/Dostoyevsky7/SkinDB_web/actions/workflows/deploy.yml/badge.svg" alt="Deploy Status"/>
  </a>
  <img src="https://img.shields.io/badge/cells-1M%2B-7c9eb8" alt="Cells"/>
  <img src="https://img.shields.io/badge/samples-600%2B-d4a574" alt="Samples"/>
  <img src="https://img.shields.io/badge/species-human%20%7C%20mouse-a3c4bc" alt="Species"/>
</p>

---

<p align="center">
  <em>Unifying the world's single-cell skin transcriptomics data in one place.</em>
</p>

## Overview

**scSAID** is a comprehensive web platform for exploring single-cell RNA sequencing (scRNA-seq) data from skin and appendage tissues. We aggregate, harmonize, and visualize data from over **1,000,000 cells** across **600+ samples**, spanning both human and mouse species in healthy and disease conditions.

Whether you're investigating keratinocyte differentiation, immune cell infiltration in psoriasis, or hair follicle stem cell dynamics — scSAID gives you the tools to explore it all.

---

## Features

| Feature | Description |
|---------|-------------|
| **Browse** | Explore datasets with advanced filtering, sorting, and pagination |
| **Gene Search** | Query expression patterns across cell types and conditions |
| **Interactive Visualization** | Embedded UMAP plots with real-time cell type exploration |
| **DEG Analysis** | Filter differentially expressed genes by p-value, fold change, and cell type |
| **CellPhoneDB** | Visualize cell-cell communication networks and receptor interactions |
| **Enrichment Analysis** | Gene set enrichment using MSigDB pathways |
| **Gene Set Scoring** | Score cells using Reactome pathways with AUCell, GSVA, or mean expression |
| **Data Export** | Download filtered results as Excel files |

---

## Tech Stack

```
┌─────────────────────────────────────────────────────────┐
│                      Frontend                           │
│   JSP  ·  CSS3  ·  JavaScript  ·  ECharts  ·  DataTables│
├─────────────────────────────────────────────────────────┤
│                      Backend                            │
│   Java Servlets  ·  Apache POI  ·  Python Dash          │
├─────────────────────────────────────────────────────────┤
│                   Infrastructure                        │
│   Maven  ·  Tomcat 9  ·  Nginx  ·  GitHub Actions       │
└─────────────────────────────────────────────────────────┘
```

---

## Quick Start

```bash
# Clone the repository
git clone https://github.com/Dostoyevsky7/SkinDB_web.git
cd SkinDB_web

# Make Maven wrapper executable
chmod +x mvnw

# Run locally
./mvnw tomcat7:run
```

Open **http://localhost:8080** in your browser.

---

## Project Structure

```
src/main/
├── java/
│   ├── Servlet/          # HTTP request handlers
│   ├── Utils/            # Data processing utilities
│   └── Entity/           # Data models
├── webapp/
│   ├── CSS/              # Design system & component styles
│   ├── JS/               # Interactive charts & micro-interactions
│   ├── images/           # Static assets
│   ├── WEB-INF/          # Configuration & data files
│   └── *.jsp             # Page templates
└── resources/
    └── mapping.json      # Field mappings
```

---

## Requirements

- Java 11+
- Maven 3.6+

---

## Production Deployment

### Dash App Integration

The enrichment and gene set scoring modules are embedded via iframes expecting same-origin paths:

- `/enrichment/` → Dash app on `localhost:8051`
- `/gene-scoring/` → Dash app on `localhost:8052`

Configure your nginx reverse proxy to forward these paths:

```nginx
map $http_upgrade $connection_upgrade {
    default upgrade;
    ''      close;
}

location /enrichment/ {
    proxy_pass http://127.0.0.1:8051/enrichment/;
    proxy_http_version 1.1;
    proxy_set_header Host $host;
    proxy_set_header Upgrade $http_upgrade;
    proxy_set_header Connection $connection_upgrade;
}

location /gene-scoring/ {
    proxy_pass http://127.0.0.1:8052/gene-scoring/;
    proxy_http_version 1.1;
    proxy_set_header Host $host;
    proxy_set_header Upgrade $http_upgrade;
    proxy_set_header Connection $connection_upgrade;
}
```

Use a process manager (systemd, supervisor, pm2) to keep Dash apps running.

---

## License

**Zhejiang University · ZJE**

---

<p align="center">
  <sub>Built with care for skin biology research</sub>
</p>
