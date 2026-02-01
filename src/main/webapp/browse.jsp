<%@ page language="java"
         contentType="text/html; charset=UTF-8"
         pageEncoding="UTF-8" %>
<%@ page import="java.io.FileInputStream, java.io.InputStream, org.apache.poi.ss.usermodel.*" %>
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Browse Datasets - scSAID</title>

    <!-- Google Fonts -->
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Cormorant+Garamond:wght@400;500;600;700&family=Source+Sans+3:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">

    <!-- Design System -->
    <link rel="stylesheet" href="CSS/design-system.css">
    <link rel="stylesheet" href="CSS/header.css">
    <link rel="stylesheet" href="CSS/animations.css">

    <style>
        /* ==========================================================================
           Browse Page Specific Styles
           ========================================================================== */

        body {
            background-color: #faf8f5;
        }

        /* Page Layout */
        .browse-page {
            min-height: 100vh;
            padding-top: 72px;
        }

        /* Page Header */
        .page-header {
            background: #1a2332;
            padding: 4rem 0;
            margin-bottom: 3rem;
        }

        .page-header__content {
            max-width: 1400px;
            margin: 0 auto;
            padding: 0 2rem;
        }

        .page-header__eyebrow {
            display: inline-block;
            font-size: 0.75rem;
            font-weight: 700;
            letter-spacing: 0.15em;
            text-transform: uppercase;
            color: #d4a574;
            margin-bottom: 1rem;
        }

        .page-header__title {
            font-family: 'Cormorant Garamond', Georgia, serif;
            font-size: clamp(2rem, 4vw, 3rem);
            font-weight: 500;
            color: #ffffff;
            margin: 0 0 1rem;
        }

        .page-header__description {
            font-size: 1.1rem;
            color: rgba(255, 255, 255, 0.7);
            max-width: 600px;
            margin: 0;
        }

        /* Main Content */
        .browse-content {
            max-width: 1400px;
            margin: 0 auto;
            padding: 0 2rem 4rem;
        }

        /* Table Card */
        .table-card {
            background: #ffffff;
            border-radius: 16px;
            box-shadow: 0 4px 12px rgba(26, 35, 50, 0.08);
            overflow: hidden;
        }

        .table-card__header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            padding: 1.5rem 2rem;
            border-bottom: 1px solid #e5e0d8;
        }

        .table-card__title {
            font-family: 'Cormorant Garamond', Georgia, serif;
            font-size: 1.5rem;
            font-weight: 500;
            color: #1a2332;
            margin: 0;
        }

        .table-card__actions {
            display: flex;
            gap: 1rem;
            align-items: center;
        }

        /* Enhanced Table */
        .data-table-wrapper {
            overflow-x: auto;
        }

        .browse-table {
            width: 100%;
            border-collapse: collapse;
            font-family: 'Source Sans 3', sans-serif;
            font-size: 0.9rem;
        }

        .browse-table thead {
            background: #1a2332;
        }

        .browse-table th {
            padding: 1rem 1.25rem;
            text-align: left;
            font-weight: 600;
            font-size: 0.75rem;
            text-transform: uppercase;
            letter-spacing: 0.08em;
            color: #ffffff;
            white-space: nowrap;
        }

        .browse-table th:first-child {
            padding-left: 2rem;
        }

        .browse-table th:last-child {
            padding-right: 2rem;
        }

        .browse-table td {
            padding: 1rem 1.25rem;
            border-bottom: 1px solid #e5e0d8;
            color: #5a6473;
            vertical-align: middle;
        }

        .browse-table td:first-child {
            padding-left: 2rem;
        }

        .browse-table td:last-child {
            padding-right: 2rem;
        }

        .browse-table tbody tr {
            transition: all 0.15s ease;
        }

        .browse-table tbody tr:hover {
            background-color: #faf8f5;
        }

        .browse-table tbody tr:last-child td {
            border-bottom: none;
        }

        /* Cell styling */
        .browse-table .cell-id {
            font-family: 'JetBrains Mono', monospace;
            font-size: 0.85rem;
            color: #1a2332;
            font-weight: 500;
        }

        .browse-table .cell-link {
            color: #e8927c;
            font-weight: 600;
            text-decoration: none;
            display: inline-flex;
            align-items: center;
            gap: 0.5rem;
            transition: all 0.15s ease;
        }

        .browse-table .cell-link:hover {
            color: #d4755d;
        }

        .browse-table .cell-link svg {
            width: 16px;
            height: 16px;
            transition: transform 0.15s ease;
        }

        .browse-table .cell-link:hover svg {
            transform: translateX(3px);
        }

        /* Species Badge */
        .species-badge {
            display: inline-flex;
            align-items: center;
            padding: 0.25rem 0.75rem;
            font-size: 0.75rem;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            border-radius: 20px;
        }

        .species-badge--human {
            background: rgba(232, 146, 124, 0.15);
            color: #d4755d;
        }

        .species-badge--mouse {
            background: rgba(212, 165, 116, 0.2);
            color: #b8864a;
        }

        /* Checkbox */
        .browse-table input[type="checkbox"] {
            width: 18px;
            height: 18px;
            accent-color: #e8927c;
            cursor: pointer;
        }

        /* Selected row */
        .browse-table tbody tr.selected-row {
            background-color: rgba(232, 146, 124, 0.08);
        }

        .browse-table tbody tr.selected-row td:first-child {
            box-shadow: inset 3px 0 0 #e8927c;
        }

        /* Pagination */
        .table-card__footer {
            display: flex;
            justify-content: center;
            align-items: center;
            padding: 1.5rem 2rem;
            border-top: 1px solid #e5e0d8;
            background: #faf8f5;
        }

        .pagination {
            display: flex;
            align-items: center;
            gap: 0.5rem;
            margin: 0;
        }

        .pagination__btn {
            display: inline-flex;
            align-items: center;
            justify-content: center;
            min-width: 40px;
            height: 40px;
            padding: 0 1rem;
            font-family: 'Source Sans 3', sans-serif;
            font-size: 0.9rem;
            font-weight: 500;
            color: #5a6473;
            background: #ffffff;
            border: 1px solid #e5e0d8;
            border-radius: 8px;
            cursor: pointer;
            text-decoration: none;
            transition: all 0.15s ease;
        }

        .pagination__btn:hover:not(.pagination__btn--disabled) {
            color: #e8927c;
            border-color: #e8927c;
        }

        .pagination__btn--disabled {
            opacity: 0.4;
            cursor: not-allowed;
        }

        .pagination__input-group {
            display: flex;
            align-items: center;
            gap: 0.5rem;
            margin: 0 0.5rem;
        }

        .pagination__label {
            font-size: 0.9rem;
            color: #5a6473;
        }

        .pagination__input {
            width: 60px;
            height: 40px;
            text-align: center;
            font-family: 'Source Sans 3', sans-serif;
            font-size: 0.9rem;
            border: 1px solid #e5e0d8;
            border-radius: 8px;
            outline: none;
            transition: all 0.15s ease;
        }

        .pagination__input:focus {
            border-color: #e8927c;
            box-shadow: 0 0 0 3px rgba(232, 146, 124, 0.15);
        }

        /* UMAP Result Container */
        .umap-container {
            margin-top: 2rem;
            padding: 2rem;
            background: #ffffff;
            border-radius: 16px;
            box-shadow: 0 4px 12px rgba(26, 35, 50, 0.08);
            text-align: center;
        }

        .umap-container__title {
            font-family: 'Cormorant Garamond', Georgia, serif;
            font-size: 1.25rem;
            font-weight: 500;
            color: #1a2332;
            margin-bottom: 1rem;
        }

        .umap-container img {
            max-width: 100%;
            border-radius: 8px;
            border: 1px solid #e5e0d8;
        }

        /* Loading Indicator */
        .loading-indicator {
            display: none;
            flex-direction: column;
            align-items: center;
            gap: 1rem;
            padding: 3rem;
        }

        .loading-indicator__spinner {
            width: 48px;
            height: 48px;
            border: 3px solid #e5e0d8;
            border-top-color: #e8927c;
            border-radius: 50%;
            animation: spin 1s linear infinite;
        }

        @keyframes spin {
            to { transform: rotate(360deg); }
        }

        .loading-indicator__text {
            font-size: 0.95rem;
            color: #5a6473;
        }

        /* Error Message */
        .error-message {
            padding: 1rem 1.5rem;
            background: rgba(220, 53, 69, 0.1);
            color: #dc3545;
            border-radius: 8px;
            margin: 2rem;
            text-align: center;
        }

        /* Filter Bar */
        .filter-bar {
            display: flex;
            flex-wrap: wrap;
            gap: 1rem;
            padding: 1.5rem 2rem;
            background: #f5f3f0;
            border-bottom: 1px solid #e5e0d8;
            align-items: flex-end;
        }

        .filter-group {
            display: flex;
            flex-direction: column;
            gap: 0.4rem;
        }

        .filter-group__label {
            font-size: 0.75rem;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            color: #5a6473;
        }

        .filter-group__select {
            min-width: 160px;
            padding: 0.65rem 2.5rem 0.65rem 1rem;
            font-family: 'Source Sans 3', sans-serif;
            font-size: 0.9rem;
            color: #1a2332;
            background: #ffffff url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='12' height='12' viewBox='0 0 24 24' fill='none' stroke='%235a6473' stroke-width='2'%3E%3Cpath d='M6 9l6 6 6-6'/%3E%3C/svg%3E") no-repeat right 1rem center;
            border: 1px solid #e5e0d8;
            border-radius: 8px;
            cursor: pointer;
            appearance: none;
            -webkit-appearance: none;
            transition: all 0.15s ease;
        }

        .filter-group__select:hover {
            border-color: #d4a574;
        }

        .filter-group__select:focus {
            outline: none;
            border-color: #e8927c;
            box-shadow: 0 0 0 3px rgba(232, 146, 124, 0.15);
        }

        .filter-bar__actions {
            display: flex;
            gap: 0.75rem;
            margin-left: auto;
            align-items: flex-end;
        }

        .filter-bar__btn {
            padding: 0.65rem 1.25rem;
            font-family: 'Source Sans 3', sans-serif;
            font-size: 0.85rem;
            font-weight: 600;
            border-radius: 8px;
            cursor: pointer;
            transition: all 0.15s ease;
        }

        .filter-bar__btn--clear {
            background: transparent;
            color: #5a6473;
            border: 1px solid #e5e0d8;
        }

        .filter-bar__btn--clear:hover {
            color: #1a2332;
            border-color: #5a6473;
        }

        .filter-count {
            font-size: 0.9rem;
            color: #5a6473;
            padding: 0.65rem 0;
        }

        .filter-count strong {
            color: #e8927c;
            font-weight: 600;
        }

        /* Hidden row (filtered out) */
        .browse-table tbody tr.filtered-out {
            display: none;
        }

        /* Responsive */
        @media (max-width: 768px) {
            .page-header {
                padding: 3rem 0;
            }

            .browse-content {
                padding: 0 1rem 3rem;
            }

            .table-card__header {
                flex-direction: column;
                gap: 1rem;
                align-items: flex-start;
            }

            .browse-table th,
            .browse-table td {
                padding: 0.75rem 1rem;
            }

            .browse-table th:first-child,
            .browse-table td:first-child {
                padding-left: 1rem;
            }

            .filter-bar {
                padding: 1rem;
            }

            .filter-group__select {
                min-width: 140px;
            }

            .filter-bar__actions {
                width: 100%;
                margin-left: 0;
                margin-top: 0.5rem;
            }
        }
    </style>
    <script src="JS/micro-interactions.js"></script>
</head>
<body class="content-fade-in">

<!-- Header -->
<header class="site-header">
    <div class="container">
        <a href="index.jsp" class="site-logo">scSAID</a>
        <nav class="main-nav">
            <a href="index.jsp" class="main-nav__link">Home</a>
            <a href="browse.jsp" class="main-nav__link main-nav__link--active">Browse</a>
            <a href="search.jsp" class="main-nav__link">Search</a>
            <a href="gene-search.jsp" class="main-nav__link">Gene Search</a>
            <a href="download.jsp" class="main-nav__link">Download</a>
        </nav>
        <div class="header-icons">
            <a href="https://github.com/Dostoyevsky7/SkinDB_web" target="_blank" class="header-icon-link" title="View on GitHub">
                <svg class="github-icon" viewBox="0 0 24 24" fill="currentColor">
                    <path d="M12 0c-6.626 0-12 5.373-12 12 0 5.302 3.438 9.8 8.207 11.387.599.111.793-.261.793-.577v-2.234c-3.338.726-4.033-1.416-4.033-1.416-.546-1.387-1.333-1.756-1.333-1.756-1.089-.745.083-.729.083-.729 1.205.084 1.839 1.237 1.839 1.237 1.07 1.834 2.807 1.304 3.492.997.107-.775.418-1.305.762-1.604-2.665-.305-5.467-1.334-5.467-5.931 0-1.311.469-2.381 1.236-3.221-.124-.303-.535-1.524.117-3.176 0 0 1.008-.322 3.301 1.23.957-.266 1.983-.399 3.003-.404 1.02.005 2.047.138 3.006.404 2.291-1.552 3.297-1.23 3.297-1.23.653 1.653.242 2.874.118 3.176.77.84 1.235 1.911 1.235 3.221 0 4.609-2.807 5.624-5.479 5.921.43.372.823 1.102.823 2.222v3.293c0 .319.192.694.801.576 4.765-1.589 8.199-6.086 8.199-11.386 0-6.627-5.373-12-12-12z"/>
                </svg>
            </a>
            <a href="https://zje.zju.edu.cn/zje/main.htm" target="_blank" class="header-icon-link" title="ZJE - Zhejiang University">
                <img src="images/ZJE_Logo.png" alt="ZJE - Zhejiang University" class="university-logo">
            </a>
        </div>
    </div>
</header>

<main class="browse-page">
    <!-- Page Header -->
    <div class="page-header">
        <div class="page-header__content">
            <span class="page-header__eyebrow">Data Explorer</span>
            <h1 class="page-header__title">Browse Datasets</h1>
            <p class="page-header__description">
                Explore our comprehensive collection of single-cell RNA sequencing datasets from skin and appendage tissues.
            </p>
        </div>
    </div>

    <!-- Main Content -->
    <div class="browse-content">
        <%
            String excelPath = application.getRealPath("/WEB-INF/BrowseShow.xlsx");
            InputStream input = null;
            Workbook workbook = null;
            try {
                input = new FileInputStream(excelPath);
                workbook = WorkbookFactory.create(input);
                Sheet sheet = workbook.getSheetAt(0);

                int rowsPerPage = 10;
                int totalRows = sheet.getLastRowNum();
                int totalPages = (int) Math.ceil((double) totalRows / rowsPerPage);

                String pageParam = request.getParameter("page");
                int pageNum = 1;
                try { pageNum = Integer.parseInt(pageParam); } catch(Exception ignore){}
                if (pageNum < 1) pageNum = 1;
                if (pageNum > totalPages) pageNum = totalPages;

                int startRow = (pageNum - 1) * rowsPerPage + 1;
                int endRow = Math.min(startRow + rowsPerPage - 1, totalRows);
        %>

        <div class="table-card" data-panel-enter>
            <div class="table-card__header">
                <h2 class="table-card__title">Dataset Preview</h2>
                <div class="table-card__actions">
                    <button id="integrate-button" class="btn btn--primary" data-btn-morph>
                        <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" style="margin-right: 8px;">
                            <circle cx="12" cy="12" r="10"></circle>
                            <path d="M12 6v12M6 12h12"></path>
                        </svg>
                        Generate Integrated UMAP
                    </button>
                </div>
            </div>

            <!-- Filter Bar -->
            <div class="filter-bar">
                <div class="filter-group">
                    <label class="filter-group__label">Species</label>
                    <select id="filter-species" class="filter-group__select">
                        <option value="">All Species</option>
                        <option value="human">Human</option>
                        <option value="mouse">Mouse</option>
                    </select>
                </div>

                <div class="filter-group">
                    <label class="filter-group__label">Disease Status</label>
                    <select id="filter-disease" class="filter-group__select">
                        <option value="">All Status</option>
                        <option value="healthy">Healthy</option>
                        <option value="disease">Disease</option>
                    </select>
                </div>

                <div class="filter-group">
                    <label class="filter-group__label">Tissue</label>
                    <select id="filter-tissue" class="filter-group__select">
                        <option value="">All Tissues</option>
                    </select>
                </div>

                <div class="filter-bar__actions">
                    <span id="filter-count" class="filter-count"></span>
                    <button id="clear-filters" class="filter-bar__btn filter-bar__btn--clear">Clear Filters</button>
                </div>
            </div>

            <div class="data-table-wrapper">
                <table class="browse-table">
                    <thead>
                    <tr>
                        <th><input type="checkbox" id="select-all" title="Select all"></th>
                        <th>SAID</th>
                        <th>GSE</th>
                        <th>GSM</th>
                        <th>Species</th>
                        <th>Disease</th>
                        <th>Tissue</th>
                        <th>Details</th>
                    </tr>
                    </thead>
                    <tbody data-stagger-group data-stagger-type="fade-up">
                    <%
                        for (int r = startRow; r <= endRow; r++) {
                            Row row = sheet.getRow(r);
                            if (row == null) continue;
                            String said_display = row.getCell(0, Row.MissingCellPolicy.CREATE_NULL_AS_BLANK).toString();
                            String gse = row.getCell(1, Row.MissingCellPolicy.CREATE_NULL_AS_BLANK).toString();
                            String gsm_value = row.getCell(2, Row.MissingCellPolicy.CREATE_NULL_AS_BLANK).toString();
                            String species = row.getCell(3, Row.MissingCellPolicy.CREATE_NULL_AS_BLANK).toString();
                            String disease = row.getCell(4, Row.MissingCellPolicy.CREATE_NULL_AS_BLANK).toString();
                            String tissue = row.getCell(5, Row.MissingCellPolicy.CREATE_NULL_AS_BLANK).toString();
                            String speciesLower = species.toLowerCase();
                            String diseaseLower = disease.toLowerCase();
                            String tissueLower = tissue.toLowerCase().trim();
                    %>
                    <tr data-species="<%= speciesLower %>" data-disease="<%= diseaseLower %>" data-tissue="<%= tissueLower %>" data-stagger-item>
                        <td><input type="checkbox" name="dataset_checkbox" value="<%= gsm_value %>"></td>
                        <td class="cell-id"><%= said_display %></td>
                        <td><%= gse %></td>
                        <td><%= gsm_value %></td>
                        <td>
                            <span class="species-badge <%= speciesLower.contains("human") ? "species-badge--human" : "species-badge--mouse" %>">
                                <%= species %>
                            </span>
                        </td>
                        <td><%= disease %></td>
                        <td><%= tissue %></td>
                        <td>
                            <a href="details.jsp?said=<%= java.net.URLEncoder.encode(said_display, "UTF-8") %>" class="cell-link">
                                View
                                <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                                    <path d="M5 12h14M12 5l7 7-7 7"></path>
                                </svg>
                            </a>
                        </td>
                    </tr>
                    <% } %>
                    </tbody>
                </table>
            </div>

            <div class="table-card__footer">
                <div class="pagination">
                    <% if (pageNum > 1) { %>
                    <a href="?page=<%= pageNum - 1 %>" class="pagination__btn">
                        <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <path d="M15 18l-6-6 6-6"></path>
                        </svg>
                        Previous
                    </a>
                    <% } else { %>
                    <span class="pagination__btn pagination__btn--disabled">
                        <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <path d="M15 18l-6-6 6-6"></path>
                        </svg>
                        Previous
                    </span>
                    <% } %>

                    <div class="pagination__input-group">
                        <span class="pagination__label">Page</span>
                        <form method="get" style="display: inline-flex; align-items: center; gap: 0.5rem;">
                            <input type="number" name="page" min="1" max="<%= totalPages %>" value="<%= pageNum %>" class="pagination__input">
                            <span class="pagination__label">of <%= totalPages %></span>
                            <button type="submit" class="pagination__btn">Go</button>
                        </form>
                    </div>

                    <% if (pageNum < totalPages) { %>
                    <a href="?page=<%= pageNum + 1 %>" class="pagination__btn">
                        Next
                        <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <path d="M9 18l6-6-6-6"></path>
                        </svg>
                    </a>
                    <% } else { %>
                    <span class="pagination__btn pagination__btn--disabled">
                        Next
                        <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                            <path d="M9 18l6-6-6-6"></path>
                        </svg>
                    </span>
                    <% } %>
                </div>
            </div>
        </div>

        <!-- UMAP Result Container -->
        <div id="umap-result-container" class="umap-container" style="display: none;">
            <div id="loading-indicator" class="loading-indicator">
                <div class="loading-indicator__spinner"></div>
                <p class="loading-indicator__text">Generating UMAP... This may take a few moments.</p>
            </div>
            <img id="umap-image" src="" alt="Integrated UMAP plot" style="display: none;">
        </div>

        <%
            } catch (Exception e) {
        %>
        <div class="error-message">
            <strong>Error loading data:</strong> <%= e.getMessage() %>
        </div>
        <%
            } finally {
                if (workbook != null) try { workbook.close(); } catch(Exception ignore){}
                if (input != null) try { input.close(); } catch(Exception ignore){}
            }
        %>
    </div>
</main>

<script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
<script>
$(document).ready(function() {
    // ===== Filter Functionality =====
    const $filterSpecies = $('#filter-species');
    const $filterDisease = $('#filter-disease');
    const $filterTissue = $('#filter-tissue');
    const $filterCount = $('#filter-count');
    const $clearFilters = $('#clear-filters');
    const $tableRows = $('.browse-table tbody tr');

    // Populate tissue options from table data
    const tissues = new Set();
    $tableRows.each(function() {
        const tissue = $(this).data('tissue');
        if (tissue) tissues.add(tissue);
    });

    // Sort and add tissue options
    Array.from(tissues).sort().forEach(function(tissue) {
        if (tissue.trim()) {
            $filterTissue.append('<option value="' + tissue + '">' + capitalizeFirst(tissue) + '</option>');
        }
    });

    function capitalizeFirst(str) {
        return str.charAt(0).toUpperCase() + str.slice(1);
    }

    function applyFilters() {
        const speciesVal = $filterSpecies.val().toLowerCase();
        const diseaseVal = $filterDisease.val().toLowerCase();
        const tissueVal = $filterTissue.val().toLowerCase();

        let visibleCount = 0;
        const totalCount = $tableRows.length;

        $tableRows.each(function() {
            const $row = $(this);
            const rowSpecies = ($row.data('species') || '').toLowerCase();
            const rowDisease = ($row.data('disease') || '').toLowerCase();
            const rowTissue = ($row.data('tissue') || '').toLowerCase();

            let show = true;

            // Species filter
            if (speciesVal && !rowSpecies.includes(speciesVal)) {
                show = false;
            }

            // Disease filter
            if (diseaseVal) {
                if (diseaseVal === 'healthy' && !rowDisease.includes('healthy') && !rowDisease.includes('normal')) {
                    show = false;
                } else if (diseaseVal === 'disease' && (rowDisease.includes('healthy') || rowDisease.includes('normal') || rowDisease === '')) {
                    show = false;
                }
            }

            // Tissue filter
            if (tissueVal && rowTissue !== tissueVal) {
                show = false;
            }

            $row.toggleClass('filtered-out', !show);
            if (show) visibleCount++;
        });

        // Update count
        if (speciesVal || diseaseVal || tissueVal) {
            $filterCount.html('Showing <strong>' + visibleCount + '</strong> of ' + totalCount);
        } else {
            $filterCount.html('');
        }
    }

    // Filter event listeners
    $filterSpecies.on('change', applyFilters);
    $filterDisease.on('change', applyFilters);
    $filterTissue.on('change', applyFilters);

    // Clear filters
    $clearFilters.on('click', function() {
        $filterSpecies.val('');
        $filterDisease.val('');
        $filterTissue.val('');
        applyFilters();
    });

    // ===== Row Selection =====
    // Row selection highlighting
    document.querySelector('#select-all').closest('table').addEventListener('change', function (e) {
        const cb = e.target;
        if (cb.type !== 'checkbox' || cb.name !== 'dataset_checkbox') return;
        const tr = cb.closest('tr');
        tr.classList.toggle('selected-row', cb.checked);
    });

    // Select all checkbox (only select visible rows)
    $('#select-all').on('click', function() {
        const flag = this.checked;
        $('input[name="dataset_checkbox"]').each(function () {
            const $row = $(this).closest('tr');
            if (!$row.hasClass('filtered-out')) {
                this.checked = flag;
                $row.toggleClass('selected-row', flag);
            }
        });
    });

    // Integration button
    $('#integrate-button').on('click', function() {
        let selectedSaids = [];
        $('input[name="dataset_checkbox"]:checked').each(function() {
            selectedSaids.push($(this).val());
        });

        if (selectedSaids.length < 2) {
            alert('Please select at least two datasets to integrate.');
            return;
        }

        $('#umap-result-container').show();
        $('#loading-indicator').css('display', 'flex');
        $('#umap-image').hide();

        $.ajax({
            url: 'integrate',
            type: 'POST',
            data: { 'saids[]': selectedSaids },
            success: function(response) {
                $('#loading-indicator').hide();
                if (response.redirectUrl) {
                    window.open(response.redirectUrl, '_blank');
                } else if (response.error) {
                    alert('Error: ' + response.error);
                } else {
                    alert('An unknown error occurred.');
                }
            },
            error: function(xhr, status, error) {
                $('#loading-indicator').hide();
                alert('AJAX error: ' + error + '\n' + xhr.responseText);
            }
        });
    });
});
</script>

</body>
</html>
