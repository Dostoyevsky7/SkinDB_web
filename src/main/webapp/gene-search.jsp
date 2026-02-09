<%@ page language="java"
         contentType="text/html; charset=UTF-8"
         pageEncoding="UTF-8" %>
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Gene Search - scSAID</title>

    <!-- Google Fonts -->
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Cormorant+Garamond:wght@400;500;600;700&family=Source+Sans+3:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">

    <!-- Design System -->
    <link rel="stylesheet" href="CSS/design-system.css">
    <link rel="stylesheet" href="CSS/header.css">
    <link rel="stylesheet" href="CSS/construction-modal-simple.css">

    <style>
        body {
            background-color: #faf8f5;
        }

        .gene-search-page {
            min-height: 100vh;
            padding-top: 72px;
        }

        /* Hero Section */
        .search-hero {
            background: #1a2332;
            padding: 5rem 0;
            text-align: center;
        }

        .search-hero__content {
            max-width: 800px;
            margin: 0 auto;
            padding: 0 2rem;
        }

        .search-hero__eyebrow {
            display: inline-block;
            font-size: 0.75rem;
            font-weight: 700;
            letter-spacing: 0.15em;
            text-transform: uppercase;
            color: #d4a574;
            margin-bottom: 1.5rem;
        }

        .search-hero__title {
            font-family: 'Cormorant Garamond', Georgia, serif;
            font-size: clamp(2.5rem, 5vw, 3.5rem);
            font-weight: 500;
            color: #ffffff;
            margin: 0 0 1rem;
        }

        .search-hero__description {
            font-size: 1.15rem;
            color: rgba(255, 255, 255, 0.7);
            margin: 0 0 2.5rem;
            line-height: 1.7;
        }

        /* Search Box */
        .search-box {
            display: flex;
            max-width: 600px;
            margin: 0 auto;
            background: #ffffff;
            border-radius: 12px;
            overflow: hidden;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.15);
        }

        .search-box__input {
            flex: 1;
            padding: 1.25rem 1.5rem;
            font-family: 'Source Sans 3', sans-serif;
            font-size: 1.1rem;
            border: none;
            outline: none;
            background: transparent;
        }

        .search-box__input::placeholder {
            color: #9ca3af;
        }

        .search-box__btn {
            padding: 1.25rem 2rem;
            background: #e8927c;
            color: #ffffff;
            font-family: 'Source Sans 3', sans-serif;
            font-size: 1rem;
            font-weight: 600;
            border: none;
            cursor: pointer;
            transition: background 0.2s ease;
            display: flex;
            align-items: center;
            gap: 0.5rem;
        }

        .search-box__btn:hover {
            background: #d4755d;
        }

        .search-box__btn:disabled {
            background: #ccc;
            cursor: not-allowed;
        }

        /* Quick Examples */
        .search-examples {
            margin-top: 1.5rem;
            display: flex;
            flex-wrap: wrap;
            justify-content: center;
            gap: 0.5rem;
        }

        .search-examples__label {
            color: rgba(255, 255, 255, 0.5);
            font-size: 0.85rem;
            margin-right: 0.5rem;
        }

        .search-examples__btn {
            padding: 0.35rem 0.75rem;
            font-family: 'JetBrains Mono', monospace;
            font-size: 0.8rem;
            color: rgba(255, 255, 255, 0.7);
            background: rgba(255, 255, 255, 0.1);
            border: 1px solid rgba(255, 255, 255, 0.2);
            border-radius: 4px;
            cursor: pointer;
            transition: all 0.15s ease;
        }

        .search-examples__btn:hover {
            color: #ffffff;
            background: rgba(255, 255, 255, 0.15);
            border-color: rgba(255, 255, 255, 0.3);
        }

        /* Results Section */
        .results-section {
            max-width: 1400px;
            margin: 0 auto;
            padding: 3rem 2rem 4rem;
        }

        /* Results Header */
        .results-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 1.5rem;
        }

        .results-header__title {
            font-family: 'Cormorant Garamond', Georgia, serif;
            font-size: 1.75rem;
            font-weight: 500;
            color: #1a2332;
            margin: 0;
        }

        .results-header__count {
            font-size: 0.95rem;
            color: #5a6473;
        }

        .results-header__count strong {
            color: #e8927c;
            font-weight: 600;
        }

        /* Results Card */
        .results-card {
            background: #ffffff;
            border-radius: 16px;
            box-shadow: 0 4px 12px rgba(26, 35, 50, 0.08);
            overflow: hidden;
        }

        /* Results Table */
        .results-table {
            width: 100%;
            border-collapse: collapse;
            font-family: 'Source Sans 3', sans-serif;
            font-size: 0.9rem;
        }

        .results-table thead {
            background: #1a2332;
        }

        .results-table th {
            padding: 1rem 1.25rem;
            text-align: left;
            font-weight: 600;
            font-size: 0.75rem;
            text-transform: uppercase;
            letter-spacing: 0.08em;
            color: #ffffff;
            white-space: nowrap;
        }

        .results-table td {
            padding: 1rem 1.25rem;
            border-bottom: 1px solid #e5e0d8;
            color: #5a6473;
            vertical-align: middle;
        }

        .results-table tbody tr {
            transition: background 0.15s ease;
        }

        .results-table tbody tr:hover {
            background: #faf8f5;
        }

        .results-table tbody tr:last-child td {
            border-bottom: none;
        }

        /* Gene cell */
        .cell-gene {
            font-family: 'JetBrains Mono', monospace;
            font-size: 0.9rem;
            font-weight: 600;
            color: #1a2332;
        }

        .cell-gene mark {
            background: rgba(232, 146, 124, 0.25);
            color: #1a2332;
            padding: 0.1em 0.2em;
            border-radius: 2px;
        }

        /* Gene link */
        .gene-link {
            color: #1a2332;
            text-decoration: none;
            border-bottom: 2px solid transparent;
            transition: all 0.2s ease;
            display: inline-flex;
            align-items: center;
            gap: 0.3rem;
        }

        .gene-link:hover {
            color: #e8927c;
            border-bottom-color: #e8927c;
        }

        .gene-link-icon {
            opacity: 0;
            transition: opacity 0.2s ease;
            font-size: 0.8rem;
        }

        .gene-link:hover .gene-link-icon {
            opacity: 1;
        }

        /* Dataset link */
        .cell-link {
            color: #e8927c;
            font-weight: 600;
            text-decoration: none;
            transition: color 0.15s ease;
        }

        .cell-link:hover {
            color: #d4755d;
            text-decoration: underline;
        }

        /* Expression badges */
        .expression-badge {
            display: inline-block;
            padding: 0.2rem 0.5rem;
            font-size: 0.75rem;
            font-weight: 600;
            border-radius: 4px;
            font-family: 'JetBrains Mono', monospace;
        }

        .expression-badge--up {
            background: rgba(34, 197, 94, 0.15);
            color: #16a34a;
        }

        .expression-badge--down {
            background: rgba(239, 68, 68, 0.15);
            color: #dc2626;
        }

        /* P-value cell */
        .cell-pval {
            font-family: 'JetBrains Mono', monospace;
            font-size: 0.85rem;
        }

        /* Group badge */
        .group-badge {
            display: inline-block;
            padding: 0.25rem 0.6rem;
            font-size: 0.75rem;
            font-weight: 500;
            background: rgba(26, 35, 50, 0.08);
            color: #5a6473;
            border-radius: 4px;
        }

        /* Empty State */
        .empty-state {
            padding: 4rem 2rem;
            text-align: center;
        }

        .empty-state__icon {
            width: 80px;
            height: 80px;
            margin: 0 auto 1.5rem;
            color: #d1d5db;
        }

        .empty-state__title {
            font-family: 'Cormorant Garamond', Georgia, serif;
            font-size: 1.5rem;
            font-weight: 500;
            color: #1a2332;
            margin: 0 0 0.5rem;
        }

        .empty-state__text {
            font-size: 1rem;
            color: #5a6473;
            margin: 0;
        }

        /* Loading State */
        .loading-state {
            padding: 4rem 2rem;
            text-align: center;
            display: none;
        }

        .loading-state.active {
            display: block;
        }

        .loading-spinner {
            width: 48px;
            height: 48px;
            border: 3px solid #e5e0d8;
            border-top-color: #e8927c;
            border-radius: 50%;
            animation: spin 1s linear infinite;
            margin: 0 auto 1rem;
        }

        @keyframes spin {
            to { transform: rotate(360deg); }
        }

        .loading-state__text {
            font-size: 1rem;
            color: #5a6473;
        }

        /* Checkbox column */
        .cell-checkbox {
            text-align: center;
            width: 40px;
        }

        .gene-checkbox {
            width: 18px;
            height: 18px;
            cursor: pointer;
            accent-color: #e8927c;
        }

        /* Visualization Panel */
        .viz-panel {
            margin-top: 2rem;
            background: #ffffff;
            border-radius: 16px;
            box-shadow: 0 4px 12px rgba(26, 35, 50, 0.08);
            overflow: hidden;
            display: none;
        }

        .viz-panel.active {
            display: block;
        }

        .viz-panel__header {
            padding: 1.5rem 2rem;
            background: linear-gradient(135deg, #1a2332 0%, #2a3342 100%);
            display: flex;
            justify-content: space-between;
            align-items: center;
            border-bottom: 2px solid #e8927c;
        }

        .viz-panel__title {
            font-family: 'Cormorant Garamond', Georgia, serif;
            font-size: 1.5rem;
            font-weight: 500;
            color: #ffffff;
            margin: 0;
        }

        .viz-panel__count {
            font-size: 0.9rem;
            color: rgba(255, 255, 255, 0.7);
            font-family: 'Source Sans 3', sans-serif;
        }

        .viz-panel__controls {
            padding: 1.5rem 2rem;
            background: #faf8f5;
            border-bottom: 1px solid #e5e0d8;
            display: flex;
            justify-content: space-between;
            align-items: center;
            flex-wrap: wrap;
            gap: 1rem;
        }

        .viz-panel__selected {
            font-size: 0.9rem;
            color: #5a6473;
            flex: 1;
        }

        .viz-panel__selected strong {
            color: #1a2332;
            font-weight: 600;
        }

        .viz-panel__genes {
            display: inline;
            font-family: 'JetBrains Mono', monospace;
            font-size: 0.85rem;
            color: #e8927c;
        }

        .viz-panel__actions {
            display: flex;
            gap: 0.75rem;
        }

        .viz-btn {
            padding: 0.65rem 1.5rem;
            font-family: 'Source Sans 3', sans-serif;
            font-size: 0.9rem;
            font-weight: 600;
            border: none;
            border-radius: 8px;
            cursor: pointer;
            transition: all 0.2s ease;
            display: inline-flex;
            align-items: center;
            gap: 0.5rem;
        }

        .viz-btn--primary {
            background: #e8927c;
            color: #ffffff;
        }

        .viz-btn--primary:hover {
            background: #d4755d;
            transform: translateY(-1px);
            box-shadow: 0 4px 12px rgba(232, 146, 124, 0.3);
        }

        .viz-btn--primary:disabled {
            background: #ccc;
            cursor: not-allowed;
            transform: none;
        }

        .viz-btn--secondary {
            background: #ffffff;
            color: #5a6473;
            border: 1px solid #e5e0d8;
        }

        .viz-btn--secondary:hover {
            background: #faf8f5;
            border-color: #d4a574;
            color: #1a2332;
        }

        .viz-panel__iframe-container {
            position: relative;
            width: 100%;
            height: 900px;
            background: #faf8f5;
        }

        .viz-panel__iframe {
            width: 100%;
            height: 100%;
            border: none;
        }

        .viz-panel__loading {
            position: absolute;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            background: rgba(250, 248, 245, 0.95);
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            z-index: 10;
        }

        .viz-panel__loading.hidden {
            display: none;
        }

        /* Responsive */
        @media (max-width: 768px) {
            .search-hero {
                padding: 3rem 0;
            }

            .search-box {
                flex-direction: column;
            }

            .search-box__btn {
                justify-content: center;
            }

            .results-header {
                flex-direction: column;
                align-items: flex-start;
                gap: 0.5rem;
            }

            .results-table th,
            .results-table td {
                padding: 0.75rem 1rem;
            }

            .viz-panel__controls {
                flex-direction: column;
                align-items: flex-start;
            }

            .viz-panel__actions {
                width: 100%;
            }

            .viz-btn {
                flex: 1;
                justify-content: center;
            }

            .viz-panel__iframe-container {
                height: 700px;
            }
        }
    </style>
</head>
<body>

<!-- Header -->
<header class="site-header">
    <div class="container">
        <a href="index.jsp" class="site-logo">scSAID</a>
        <nav class="main-nav">
            <a href="index.jsp" class="main-nav__link">Home</a>
            <a href="browse.jsp" class="main-nav__link">Browse</a>
            <a href="search.jsp" class="main-nav__link">Search</a>
            <a href="gene-search.jsp" class="main-nav__link main-nav__link--active">Gene Search</a>
            <a href="download.jsp" class="main-nav__link">Download</a>
        </nav>
        <a href="https://zje.zju.edu.cn/zje/main.htm" target="_blank">
            <img src="images/ZJE_Logo.png" alt="ZJE - Zhejiang University" class="university-logo">
        </a>
    </div>
</header>

<main class="gene-search-page">
    <!-- Hero Section -->
    <section class="search-hero">
        <div class="search-hero__content">
            <span class="search-hero__eyebrow">Gene Explorer</span>
            <h1 class="search-hero__title">Search Genes Across All Datasets</h1>
            <p class="search-hero__description">
                Query differentially expressed genes across our entire collection of single-cell RNA sequencing datasets.
            </p>

            <div class="search-box">
                <input type="text" id="gene-query" class="search-box__input"
                       placeholder="Enter gene name (e.g., KRT14, COL1A1, ACTA2)">
                <button id="search-btn" class="search-box__btn">
                    <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                        <circle cx="11" cy="11" r="8"></circle>
                        <path d="M21 21l-4.35-4.35"></path>
                    </svg>
                    Search
                </button>
            </div>

            <div class="search-examples">
                <span class="search-examples__label">Try:</span>
                <button class="search-examples__btn" data-gene="KRT14">KRT14</button>
                <button class="search-examples__btn" data-gene="COL1A1">COL1A1</button>
                <button class="search-examples__btn" data-gene="ACTA2">ACTA2</button>
                <button class="search-examples__btn" data-gene="PECAM1">PECAM1</button>
                <button class="search-examples__btn" data-gene="CD3E">CD3E</button>
            </div>
        </div>
    </section>

    <!-- Results Section -->
    <section class="results-section">
        <!-- Loading State -->
        <div id="loading-state" class="loading-state">
            <div class="loading-spinner"></div>
            <p class="loading-state__text">Searching across all datasets...</p>
        </div>

        <!-- Results Container -->
        <div id="results-container" style="display: none;">
            <div class="results-header">
                <h2 class="results-header__title">Search Results</h2>
                <span id="results-count" class="results-header__count"></span>
            </div>

            <div class="results-card">
                <div style="overflow-x: auto;">
                    <table class="results-table">
                        <thead>
                            <tr>
                                <th class="cell-checkbox">
                                    <input type="checkbox" id="select-all" class="gene-checkbox" title="Select all">
                                </th>
                                <th>Gene</th>
                                <th>Dataset</th>
                                <th>GSE</th>
                                <th>Cell Type</th>
                                <th>Log2 FC</th>
                                <th>Adj. P-value</th>
                                <th>Action</th>
                            </tr>
                        </thead>
                        <tbody id="results-body">
                        </tbody>
                    </table>
                </div>
            </div>

            <!-- Visualization Panel -->
            <div id="viz-panel" class="viz-panel">
                <div class="viz-panel__header">
                    <h3 class="viz-panel__title">Gene Expression on Integrated UMAP</h3>
                    <span id="viz-gene-count" class="viz-panel__count"></span>
                </div>

                <div class="viz-panel__controls">
                    <div class="viz-panel__selected">
                        <strong>Selected Genes:</strong>
                        <span id="selected-genes-display" class="viz-panel__genes">None</span>
                    </div>

                    <div class="viz-panel__actions">
                        <button id="clear-selection-btn" class="viz-btn viz-btn--secondary">
                            <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                <line x1="18" y1="6" x2="6" y2="18"></line>
                                <line x1="6" y1="6" x2="18" y2="18"></line>
                            </svg>
                            Clear Selection
                        </button>
                        <button id="visualize-btn" class="viz-btn viz-btn--primary" disabled>
                            <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                <circle cx="12" cy="12" r="10"></circle>
                                <circle cx="12" cy="12" r="3"></circle>
                            </svg>
                            Visualize on UMAP
                        </button>
                    </div>
                </div>

                <div id="viz-iframe-container" class="viz-panel__iframe-container" style="display: none;">
                    <div id="viz-loading" class="viz-panel__loading">
                        <div class="loading-spinner"></div>
                        <p class="loading-state__text">Loading visualization...</p>
                    </div>
                    <iframe id="viz-iframe" class="viz-panel__iframe" src=""></iframe>
                </div>
            </div>
        </div>

        <!-- Empty/Initial State -->
        <div id="empty-state" class="results-card">
            <div class="empty-state">
                <svg class="empty-state__icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5">
                    <path d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" stroke-linecap="round" stroke-linejoin="round"/>
                </svg>
                <h3 class="empty-state__title">Enter a gene name to search</h3>
                <p class="empty-state__text">Search results will appear here</p>
            </div>
        </div>

        <!-- No Results State -->
        <div id="no-results-state" class="results-card" style="display: none;">
            <div class="empty-state">
                <svg class="empty-state__icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5">
                    <path d="M9.172 9.172a4 4 0 015.656 0M9 9l6 6m-6 0l6-6" stroke-linecap="round" stroke-linejoin="round"/>
                    <circle cx="12" cy="12" r="10" stroke-linecap="round" stroke-linejoin="round"/>
                </svg>
                <h3 class="empty-state__title">No results found</h3>
                <p class="empty-state__text">Try a different gene name or check your spelling</p>
            </div>
        </div>
    </section>
</main>

<script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
<script>
$(document).ready(function() {
    const $input = $('#gene-query');
    const $searchBtn = $('#search-btn');
    const $loadingState = $('#loading-state');
    const $resultsContainer = $('#results-container');
    const $emptyState = $('#empty-state');
    const $noResultsState = $('#no-results-state');
    const $resultsBody = $('#results-body');
    const $resultsCount = $('#results-count');
    const $selectAll = $('#select-all');

    // Visualization elements
    const $vizPanel = $('#viz-panel');
    const $vizBtn = $('#visualize-btn');
    const $clearSelectionBtn = $('#clear-selection-btn');
    const $selectedGenesDisplay = $('#selected-genes-display');
    const $vizGeneCount = $('#viz-gene-count');
    const $vizIframeContainer = $('#viz-iframe-container');
    const $vizIframe = $('#viz-iframe');
    const $vizLoading = $('#viz-loading');

    let currentQuery = '';
    let selectedGenes = new Set();

    function performSearch(query) {
        if (!query || query.trim() === '') return;

        query = query.trim();
        currentQuery = query.toLowerCase();

        // Reset selection
        selectedGenes.clear();
        updateSelectionUI();

        // Show loading
        $emptyState.hide();
        $noResultsState.hide();
        $resultsContainer.hide();
        $loadingState.addClass('active');
        $searchBtn.prop('disabled', true);

        $.ajax({
            url: 'gene-search',
            type: 'GET',
            data: { q: query },
            success: function(response) {
                $loadingState.removeClass('active');
                $searchBtn.prop('disabled', false);

                if (response.results && response.results.length > 0) {
                    renderResults(response.results, response.count);
                    $resultsContainer.show();
                } else {
                    $noResultsState.show();
                }
            },
            error: function(xhr, status, error) {
                $loadingState.removeClass('active');
                $searchBtn.prop('disabled', false);
                alert('Search error: ' + error);
                $emptyState.show();
            }
        });
    }

    function renderResults(results, count) {
        $resultsCount.html('Found <strong>' + count + '</strong> results' + (count >= 500 ? ' (showing top 500)' : ''));

        // Get unique genes for selection
        const uniqueGenes = new Set();
        results.forEach(row => uniqueGenes.add(row.gene));

        let html = '';
        let processedGenes = new Set();

        results.forEach(function(row) {
            const logfc = parseFloat(row.logfc);
            const isUp = logfc > 0;
            const fcClass = isUp ? 'expression-badge--up' : 'expression-badge--down';
            const fcSign = isUp ? '+' : '';

            // Highlight matching text
            const geneName = highlightMatch(row.gene, currentQuery);

            // Add checkbox only for first occurrence of each gene
            let checkboxHtml = '';
            if (!processedGenes.has(row.gene)) {
                checkboxHtml = '<input type="checkbox" class="gene-checkbox gene-select" data-gene="' + escapeHtml(row.gene) + '">';
                processedGenes.add(row.gene);
            }

            html += '<tr>';
            html += '<td class="cell-checkbox">' + checkboxHtml + '</td>';
            html += '<td class="cell-gene"><a href="gene-details?gene=' + encodeURIComponent(row.gene) + '&species=human" class="gene-link">' + geneName + '<span class="gene-link-icon">→</span></a></td>';
            html += '<td><a href="details.jsp?said=' + encodeURIComponent(row.said) + '" class="cell-link">' + row.said + '</a></td>';
            html += '<td>' + row.gse + '</td>';
            html += '<td><span class="group-badge">' + escapeHtml(row.group) + '</span></td>';
            html += '<td><span class="expression-badge ' + fcClass + '">' + fcSign + row.logfc + '</span></td>';
            html += '<td class="cell-pval">' + row.pval + '</td>';
            html += '<td><a href="details.jsp?said=' + encodeURIComponent(row.said) + '#DEG" class="cell-link">Dataset</a> | <a href="gene-details?gene=' + encodeURIComponent(row.gene) + '&species=human" class="cell-link">Gene Info</a></td>';
            html += '</tr>';
        });

        $resultsBody.html(html);

        // Attach checkbox handlers
        $('.gene-select').on('change', function() {
            const gene = $(this).data('gene');
            if ($(this).prop('checked')) {
                selectedGenes.add(gene);
            } else {
                selectedGenes.delete(gene);
            }
            updateSelectionUI();
        });

        // Update select all state
        $selectAll.prop('checked', false);
    }

    function updateSelectionUI() {
        const count = selectedGenes.size;

        if (count === 0) {
            $selectedGenesDisplay.text('None');
            $vizBtn.prop('disabled', true);
            $vizGeneCount.text('');
        } else {
            const genesList = Array.from(selectedGenes).join(', ');
            $selectedGenesDisplay.text(genesList);
            $vizBtn.prop('disabled', false);
            $vizGeneCount.text(count + ' gene' + (count > 1 ? 's' : '') + ' selected');
        }

        // Show/hide viz panel based on selection
        if (count > 0) {
            $vizPanel.addClass('active');
        }
    }

    function visualizeGenes() {
        if (selectedGenes.size === 0) return;

        const genesList = Array.from(selectedGenes).join(',');
        const vizUrl = 'http://localhost:8053/gene-viz/?genes=' + encodeURIComponent(genesList);

        // Show iframe container and loading
        $vizIframeContainer.show();
        $vizLoading.removeClass('hidden');
        $vizBtn.prop('disabled', true);

        // Set iframe source
        $vizIframe.attr('src', vizUrl);

        // Hide loading after iframe loads
        $vizIframe.on('load', function() {
            setTimeout(function() {
                $vizLoading.addClass('hidden');
                $vizBtn.prop('disabled', false);
            }, 500);
        });

        // Scroll to visualization
        setTimeout(function() {
            $vizPanel[0].scrollIntoView({ behavior: 'smooth', block: 'nearest' });
        }, 100);
    }

    function clearSelection() {
        selectedGenes.clear();
        $('.gene-select').prop('checked', false);
        $selectAll.prop('checked', false);
        updateSelectionUI();
        $vizIframeContainer.hide();
        $vizIframe.attr('src', '');
    }

    function highlightMatch(text, query) {
        const regex = new RegExp('(' + escapeRegex(query) + ')', 'gi');
        return escapeHtml(text).replace(regex, '<mark>$1</mark>');
    }

    function escapeHtml(text) {
        const div = document.createElement('div');
        div.textContent = text;
        return div.innerHTML;
    }

    function escapeRegex(string) {
        return string.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
    }

    // Search button click
    $searchBtn.on('click', function() {
        performSearch($input.val());
    });

    // Enter key
    $input.on('keypress', function(e) {
        if (e.which === 13) {
            performSearch($input.val());
        }
    });

    // Example buttons
    $('.search-examples__btn').on('click', function() {
        const gene = $(this).data('gene');
        $input.val(gene);
        performSearch(gene);
    });

    // Select all checkbox
    $selectAll.on('change', function() {
        const isChecked = $(this).prop('checked');
        $('.gene-select').each(function() {
            $(this).prop('checked', isChecked);
            const gene = $(this).data('gene');
            if (isChecked) {
                selectedGenes.add(gene);
            } else {
                selectedGenes.delete(gene);
            }
        });
        updateSelectionUI();
    });

    // Visualize button
    $vizBtn.on('click', visualizeGenes);

    // Clear selection button
    $clearSelectionBtn.on('click', clearSelection);
});
</script>

<!-- Under Construction Modal Script -->
<script src="JS/construction-modal-simple.js"></script>

</body>
</html>
