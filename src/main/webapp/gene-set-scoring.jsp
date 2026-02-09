<%@ page contentType="text/html;charset=UTF-8" language="java" %>
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Gene Set Scoring - scSAID</title>

    <!-- Google Fonts -->
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Cormorant+Garamond:ital,wght@0,400;0,500;0,600;0,700;1,400&family=Source+Sans+3:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">

    <!-- Design System -->
    <link rel="stylesheet" href="CSS/design-system.css">
    <link rel="stylesheet" href="CSS/header.css">
    
    <style>
        .analysis-container {
            padding: var(--space-2xl) 0;
            background: var(--bg-body);
        }
        
        .analysis-header {
            text-align: center;
            max-width: 800px;
            margin: 0 auto var(--space-3xl);
            padding: 0 var(--space-lg);
        }
        
        .analysis-title {
            font-family: var(--font-display);
            font-size: 2.5rem;
            font-weight: 600;
            color: var(--text-primary);
            margin-bottom: var(--space-md);
        }
        
        .analysis-subtitle {
            font-size: 1.2rem;
            color: var(--text-secondary);
            margin-bottom: var(--space-xl);
        }
        
        .iframe-container {
            width: 100%;
            max-width: 1200px;
            height: 80vh;
            margin: 0 auto;
            border-radius: var(--radius-lg);
            overflow: hidden;
            box-shadow: var(--shadow-xl);
            border: 1px solid var(--border-light);
        }
        
        .iframe-container iframe {
            width: 100%;
            height: 100%;
            border: none;
        }
        
        .analysis-description {
            max-width: 800px;
            margin: 0 auto var(--space-xl);
            padding: 0 var(--space-lg);
            text-align: center;
            color: var(--text-secondary);
            line-height: 1.7;
        }
    </style>
</head>
<body>
    <!-- Header -->
    <header class="main-header">
        <div class="header-content">
            <div class="logo">
                <a href="index.jsp">scSAID</a>
            </div>
            <nav class="main-nav">
                <a href="index.jsp">Home</a>
                <a href="browse.jsp">Browse</a>
                <a href="search.jsp">Search</a>
                <a href="gene-search.jsp">Gene Search</a>
                <a href="enrichment-analysis.jsp" class="main-nav__link">Enrichment Analysis</a>
                <a href="gene-set-scoring.jsp" class="main-nav__link main-nav__link--active">Gene Set Scoring</a>
                <a href="download.jsp">Download</a>
            </nav>
        </div>
    </header>

    <main class="analysis-container">
        <section class="analysis-header">
            <h1 class="analysis-title">Gene Set Scoring Analysis</h1>
            <p class="analysis-subtitle">Score cell types using Reactome pathway gene sets</p>
            <p class="analysis-description">
                Perform gene set scoring analysis using multiple methods: AUCell, GSVA, and Mean Expression. 
                This tool allows you to evaluate pathway activity and functional characteristics across different cell types.
            </p>
        </section>

        <section class="iframe-container">
            <iframe 
                src="/gene-scoring/" 
                title="Gene Set Scoring Tool"
                allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture"
                allowfullscreen>
            </iframe>
        </section>
    </main>

    <!-- Footer -->
    <footer class="main-footer">
        <div class="footer-content">
            <p>&copy; 2024 scSAID - Single-Cell Skin & Appendages Integrated Database</p>
            <p>Zhejiang University · ZJE</p>
        </div>
    </footer>
</body>
</html>