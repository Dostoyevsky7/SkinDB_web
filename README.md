# scSAID: <u>S</u>ingle-<u>C</u>ell <u>S</u>kin & <u>A</u>ppendages <u>I</u>ntegrated <u>D</u>atabase <img src="./src/main/webapp/images/scSAID_LOGO.png" align="right" width="120"/>

> A comprehensive web platform for exploring single-cell RNA sequencing data from skin and appendage tissues.

> The database contains data from over **1,000,000 cells** across **600+ samples** from human and mouse skin tissues in both healthy and disease conditions.



## Features

-   **Browse** — Explore datasets with filtering and pagination
-   **Search & Integrate** — Select multiple datasets for integrated UMAP visualization
-   **Dataset Details** — View sample metadata, cell clustering, DEG results, and CellPhoneDB analysis
-   **Data Export** — Download filtered results as Excel files

## Tech Stack

| Layer    | Technology                |
|----------|---------------------------|
| Frontend | JSP, CSS3, JavaScript     |
| Backend  | Java Servlets             |
| Build    | Maven                     |
| Server   | Apache Tomcat 9           |
| Data     | Apache POI, Python (Dash) |

## Quick Start

``` bash
# Clone the repository
git clone https://github.com/Dostoyevsky7/SkinDB_web.git
cd SkinDB_web

# Make Maven wrapper executable
chmod +x mvnw

# Run locally
./mvnw tomcat7:run
```

Open [**http://localhost:8080**](http://localhost:8080){.uri} in your browser.

## Project Structure

```         
src/main/
├── java/
│   ├── Servlet/        # Request handlers
│   ├── Utils/          # Data processing utilities
│   └── Entity/         # Data models
├── webapp/
│   ├── CSS/            # Stylesheets
│   ├── images/         # Static assets
│   ├── WEB-INF/        # Config & data files
│   └── *.jsp           # Page templates
└── resources/
    └── mapping.json    # Field mappings
```

## Requirements

-   Java 11+
-   Maven 3.6+

## License

Zhejiang University · ZJE

------------------------------------------------------------------------

<p align="center">

<sub>Built for skin biology research</sub>

</p>
