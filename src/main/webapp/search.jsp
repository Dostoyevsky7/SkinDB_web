<%@ page import="java.io.InputStream" %>
<%@ page import="org.apache.poi.ss.usermodel.Workbook" %>
<%@ page import="java.io.FileInputStream" %>
<%@ page import="org.apache.poi.ss.usermodel.WorkbookFactory" %>
<%@ page import="org.apache.poi.ss.usermodel.Sheet" %>
<%@ page import="org.apache.poi.ss.usermodel.Row" %>
<%@ page contentType="text/html;charset=UTF-8" %>
<%@ taglib uri="http://java.sun.com/jsp/jstl/core" prefix="c" %>
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Search & Integrate - scSAID</title>

    <!-- Google Fonts -->
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Cormorant+Garamond:wght@400;500;600;700&family=Source+Sans+3:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">

    <!-- Design System -->
    <link rel="stylesheet" href="CSS/design-system.css">
    <link rel="stylesheet" href="CSS/header.css">

    <style>
        /* ==========================================================================
           Search Page Specific Styles
           ========================================================================== */

        body {
            background-color: #faf8f5;
        }

        /* Page Layout */
        .search-page {
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
        .search-content {
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

        .table-card__subtitle {
            font-size: 0.9rem;
            color: #5a6473;
            margin-top: 0.25rem;
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

        .search-table {
            width: 100%;
            border-collapse: collapse;
            font-family: 'Source Sans 3', sans-serif;
            font-size: 0.9rem;
        }

        .search-table thead {
            background: #1a2332;
        }

        .search-table th {
            padding: 1rem 1.25rem;
            text-align: left;
            font-weight: 600;
            font-size: 0.75rem;
            text-transform: uppercase;
            letter-spacing: 0.08em;
            color: #ffffff;
            white-space: nowrap;
        }

        .search-table th:first-child {
            padding-left: 2rem;
        }

        .search-table th:last-child {
            padding-right: 2rem;
        }

        .search-table td {
            padding: 1rem 1.25rem;
            border-bottom: 1px solid #e5e0d8;
            color: #5a6473;
            vertical-align: middle;
        }

        .search-table td:first-child {
            padding-left: 2rem;
        }

        .search-table td:last-child {
            padding-right: 2rem;
        }

        .search-table tbody tr {
            transition: all 0.15s ease;
        }

        .search-table tbody tr:hover {
            background-color: #faf8f5;
        }

        .search-table tbody tr:last-child td {
            border-bottom: none;
        }

        /* Cell styling */
        .search-table .cell-id {
            font-family: 'JetBrains Mono', monospace;
            font-size: 0.85rem;
            color: #1a2332;
            font-weight: 500;
        }

        .search-table .cell-link {
            color: #e8927c;
            font-weight: 600;
            text-decoration: none;
            display: inline-flex;
            align-items: center;
            gap: 0.5rem;
            transition: all 0.15s ease;
        }

        .search-table .cell-link:hover {
            color: #d4755d;
        }

        .search-table .cell-link svg {
            width: 16px;
            height: 16px;
            transition: transform 0.15s ease;
        }

        .search-table .cell-link:hover svg {
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
        .search-table input[type="checkbox"] {
            width: 18px;
            height: 18px;
            accent-color: #e8927c;
            cursor: pointer;
        }

        /* Selected row */
        .search-table tbody tr.selected-row {
            background-color: rgba(232, 146, 124, 0.08);
        }

        .search-table tbody tr.selected-row td:first-child {
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

        /* Responsive */
        @media (max-width: 768px) {
            .page-header {
                padding: 3rem 0;
            }

            .search-content {
                padding: 0 1rem 3rem;
            }

            .table-card__header {
                flex-direction: column;
                gap: 1rem;
                align-items: flex-start;
            }

            .search-table th,
            .search-table td {
                padding: 0.75rem 1rem;
            }

            .search-table th:first-child,
            .search-table td:first-child {
                padding-left: 1rem;
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
            <a href="search.jsp" class="main-nav__link main-nav__link--active">Search</a>
            <a href="#" class="main-nav__link">Help</a>
            <a href="#" class="main-nav__link">Download</a>
        </nav>
        <a href="https://zje.zju.edu.cn/zje/main.htm" target="_blank">
            <img src="images/ZJE_Logo.png"
                 alt="ZJE - Zhejiang University" class="university-logo">
        </a>
    </div>
</header>

<main class="search-page">
    <!-- Page Header -->
    <div class="page-header">
        <div class="page-header__content">
            <span class="page-header__eyebrow">Data Integration</span>
            <h1 class="page-header__title">Search & Integrate</h1>
            <p class="page-header__description">
                Select multiple datasets to generate integrated UMAP visualizations for cross-sample analysis.
            </p>
        </div>
    </div>

    <!-- Main Content -->
    <div class="search-content">
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

        <div class="table-card">
            <div class="table-card__header">
                <div>
                    <h2 class="table-card__title">Dataset Selection</h2>
                    <p class="table-card__subtitle">Select at least 2 datasets to generate an integrated UMAP</p>
                </div>
                <div class="table-card__actions">
                    <button id="integrate-button" class="btn btn--primary">
                        <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" style="margin-right: 8px;">
                            <circle cx="12" cy="12" r="10"></circle>
                            <path d="M12 6v12M6 12h12"></path>
                        </svg>
                        Generate Integrated UMAP
                    </button>
                </div>
            </div>

            <div class="data-table-wrapper">
                <table class="search-table">
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
                    <tbody>
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
                    %>
                    <tr>
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
    // Row selection highlighting
    document.querySelector('#select-all').closest('table').addEventListener('change', function (e) {
        const cb = e.target;
        if (cb.type !== 'checkbox' || cb.name !== 'dataset_checkbox') return;
        const tr = cb.closest('tr');
        tr.classList.toggle('selected-row', cb.checked);
    });

    // Select all checkbox
    $('#select-all').on('click', function() {
        const flag = this.checked;
        $('input[name="dataset_checkbox"]').each(function () {
            this.checked = flag;
            $(this).closest('tr').toggleClass('selected-row', flag);
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
