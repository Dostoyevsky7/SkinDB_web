<%@ page language="java"
         contentType="text/html; charset=UTF-8"
         pageEncoding="UTF-8" %>
<%@ page import="
    java.io.*,
    java.util.*,
    java.nio.file.Paths,
    org.apache.poi.ss.usermodel.*" %>
<%!
    // =========================================================================
    // CONFIGURATION: PLEASE SET THESE PATHS
    // =========================================================================
    // 1. Set the absolute base path to your main data directory
    private static final String DATA_BASE_PATH = "D:\\web\\web_application\\scrna_website_test\\src\\main\\webapp\\SkinDB_New";
    // 2. Set the command for your python executable
    // WARNING: If you are using a virtual environment (like dash_env),
    // you MUST set this to the absolute path of the python.exe file.
    // Example: private static final String PYTHON_COMMAND = "D:\\path\\to\\your\\virtual\\env\\Scripts\\python.exe";
    private static final String PYTHON_COMMAND = "C:\\Anaconda3\\envs\\dash_env\\python.exe";
    // =========================================================================


    // Helper method to build the path to the cpdb_out directory for a given sample
    private String getCpdbPath(String said, String gse, String gsm) {
        // This logic assumes a path structure like /SkinDB_New/10X/human/GSE.../GSM.../
        // You MUST adjust this to match your actual directory structure.
        if (gse == null || gsm == null || gse.isEmpty() || gsm.isEmpty()) {
            return null;
        }
        // This is an example, please update it to match your structure.
        return Paths.get(DATA_BASE_PATH, "10X", "human", gse, gsm, "cpdb_out").toString();
    }
%>
<%
    // =========================================================================
    // SECTION A: NEW BACKEND LOGIC FOR HANDLING AJAX REQUESTS FOR CELLPHONEDB
    // =========================================================================
    String action = request.getParameter("action");
    if (action != null) {
        response.setContentType("application/json");
        response.setCharacterEncoding("UTF-8");
        PrintWriter jsonOut = response.getWriter();
        String saidParamForAction = request.getParameter("said");
        String gseForAction = request.getParameter("gse");
        String gsmForAction = request.getParameter("gsm");

        String cpdbPath = getCpdbPath(saidParamForAction, gseForAction, gsmForAction);

        if (cpdbPath == null || !new File(cpdbPath).exists()) {
            jsonOut.print("{\"error\": \"CPDB data path not found for the given sample. Path was: " + (cpdbPath == null ? "null" : cpdbPath.replace("\\", "\\\\")) + "\"}");
            jsonOut.flush();
            return;
        }

        try {
            if ("get_cell_types".equals(action)) {
                // 调用 Python 脚本 --list 模式
                String pythonScriptPath = Paths.get(cpdbPath, "plot_cpdb_receiver_top15.py").toString();
                ProcessBuilder pb = new ProcessBuilder(PYTHON_COMMAND, pythonScriptPath, "--list");
                pb.directory(new File(cpdbPath));
                pb.redirectErrorStream(true);
                try {
                    Process p = pb.start();
                    BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream(), "UTF-8"));
                    StringBuilder sb = new StringBuilder();
                    String line;
                    while ((line = br.readLine()) != null) {
                        sb.append(line);
                    }
                    br.close();
                    int exitCode = p.waitFor();
                    if (exitCode == 0) {
                        jsonOut.print(sb.toString());
                        // 直接输出 Python 的 JSON
                    } else {
                        jsonOut.print("{\"error\": \"Python script exited with code " + exitCode + ". Output: " + sb.toString().replace("\"", "'") + "\"}");
                    }
                } catch (Exception e) {
                    jsonOut.print("{\"error\": \"Exception: " + e.getMessage().replace("\"", "'") + "\"}");
                }
                return;
            }
            else if ("generate_plot".equals(action)) {
                String plotType = request.getParameter("plot_type");
                String scriptName = "";
                List<String> command = new ArrayList<String>();
                command.add(PYTHON_COMMAND);
                String outputFileName = "";
                // 定义输出文件名变量

                if ("summary".equals(plotType)) {
                    scriptName = "plot_cpdb_sum_sig.py";
                    outputFileName = "sum_sig_heatmap.png"; // 使用您本地运行的硬编码文件名
                } else if ("receiver".equals(plotType)) {
                    scriptName = "plot_cpdb_receiver_top15.py";
                    String cellType = request.getParameter("cell_type");
                    if (cellType == null || cellType.isEmpty()) {
                        jsonOut.print("{\"error\": \"Cell type is required for this plot.\"}");
                        return;
                    }
                    command.add(new File(cpdbPath, scriptName).getAbsolutePath());
                    command.add(cellType);
                    outputFileName = "top15_receiver_" + cellType.replace(" ", "_").replace("/", "_") + ".png";
                } else {
                    jsonOut.print("{\"error\": \"Invalid plot type specified.\"}");
                    return;
                }

                if (command.size() == 1) { // For scripts that need no extra args, like summary plot
                    command.add(new File(cpdbPath, scriptName).getAbsolutePath());
                }

                // Execute the script
                ProcessBuilder pb = new ProcessBuilder(command);
                pb.directory(new File(cpdbPath)); // Execute script in its directory
                // 修改：不捕获脚本输出，只等待进程结束，以防止JSON响应被污染
                Process process = pb.start();
                InputStream is = null, es = null;
                try {
                    is = process.getInputStream();
                    es = process.getErrorStream();
                    while (is.read() != -1) ;
                    while (es.read() != -1) ;
                } finally {
                    if (is != null) try { is.close(); } catch (IOException ignore) {}
                    if (es != null) try { es.close(); } catch (IOException ignore) {}
                }

                int exitCode = process.waitFor();
                if (exitCode == 0) {
                    // Check if the output file was generated
                    File outputFile = new File(cpdbPath, outputFileName);
                    if (outputFile.exists()) {
                        String imageUrl = request.getContextPath() + "/SkinDB_New/10X/human/" + gseForAction + "/" + gsmForAction + "/cpdb_out/" + outputFileName;
                        jsonOut.print("{\"imageUrl\": \"" + imageUrl + "\"}");
                    } else {
                        // 修改: 移除 scriptOutput 变量，因为我们已不再捕获它
                        jsonOut.print("{\"error\": \"Plot generated successfully but output file not found: " + outputFileName + ".\"}");
                    }
                } else {
                    // 修改: 移除 scriptOutput 变量，因为我们已不再捕获它
                    jsonOut.print("{\"error\": \"Failed to generate plot. Exit code: " + exitCode + ".\"}");
                }
            }
        } catch (Exception e) {
            jsonOut.print("{\"error\": \"An exception occurred: " + e.getMessage().replace("\"", "'") + "\"}");
        } finally {
            jsonOut.flush();
        }
        return; // End execution here, do not render the HTML page
    }

    // =========================================================================
    // SECTION B: EXISTING JSP LOGIC FOR PAGE DISPLAY
    // =========================================================================
    // 1) 获取 URL 中的 said 参数
    String saidParam = request.getParameter("said");
    if (saidParam == null || saidParam.trim().isEmpty()) {
        out.println("<h2 style='color:red;'>Error: no SAID specified.</h2>");
        return;
    }

    // 2) 声明变量
    String saidVal = "";
    String gseVal = "";
    String gsmVal = "";
    String speciesVal = "";
    String diseaseVal = "";
    String tissueVal = "";
    String charVal = "";
    String titleVal = "";
    String summaryVal = "";
    String designVal = "";

    DataFormatter fmt = new DataFormatter();
    Workbook workbook = null;
    InputStream fis = null;
    try {
        // 3) 打开 Excel
        String excelPath = application.getRealPath("/WEB-INF/AllData.xlsx");
        fis = new FileInputStream(excelPath);
        workbook = WorkbookFactory.create(fis);
        Sheet sheet = workbook.getSheetAt(0);
        // 4) 构建表头映射
        Row headerRow = sheet.getRow(0);
        Map<String,Integer> colIdx = new HashMap<String,Integer>();
        for (Cell c : headerRow) {
            colIdx.put(fmt.formatCellValue(c).trim(), c.getColumnIndex());
        }

        // 5) 各列索引
        Integer idxSAID    = colIdx.get("SAID");
        Integer idxGSE     = colIdx.get("GSE");
        Integer idxGSM     = colIdx.get("GSM");
        Integer idxSpecies = colIdx.get("species");
        Integer idxDisease = colIdx.get("disease information");
        Integer idxTissue  = colIdx.get("tissue information");
        Integer idxChar    = colIdx.get("Characteristics");
        Integer idxTitle   = colIdx.get("Title");
        Integer idxSummary = colIdx.get("Summary");
        Integer idxDesign  = colIdx.get("Overall Design");
        if (idxSAID == null) {
            out.println("<h2 style='color:red;'>Error: SAID column not found.</h2>");
            return;
        }

        // 6) 查找对应行
        for (int i = 1; i <= sheet.getLastRowNum(); i++) {
            Row row = sheet.getRow(i);
            if (row == null) continue;
            String v = fmt.formatCellValue(row.getCell(idxSAID));
            if (saidParam.equals(v)) {
                saidVal    = v;
                if (idxGSE     != null) gseVal     = fmt.formatCellValue(row.getCell(idxGSE));
                if (idxGSM     != null) gsmVal     = fmt.formatCellValue(row.getCell(idxGSM));
                if (idxSpecies != null) speciesVal = fmt.formatCellValue(row.getCell(idxSpecies));
                if (idxDisease != null) diseaseVal = fmt.formatCellValue(row.getCell(idxDisease));
                if (idxTissue  != null) tissueVal  = fmt.formatCellValue(row.getCell(idxTissue));
                if (idxChar    != null) charVal    = fmt.formatCellValue(row.getCell(idxChar));
                if (idxTitle   != null) titleVal   = fmt.formatCellValue(row.getCell(idxTitle));
                if (idxSummary != null) summaryVal = fmt.formatCellValue(row.getCell(idxSummary));
                if (idxDesign  != null) designVal  = fmt.formatCellValue(row.getCell(idxDesign));
                break;
            }
        }
    } catch (Exception e) {
        out.println("<h2 style='color:red;'>Error loading Excel: " + e.getMessage() + "</h2>");
        return;
    } finally {
        if (workbook != null) try { workbook.close(); } catch (Exception ignore) {}
        if (fis      != null) try { fis.close(); } catch (Exception ignore) {}
    }
%>

<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Dataset Details - scSAID</title>

    <!-- Google Fonts -->
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Cormorant+Garamond:wght@400;500;600;700&family=Source+Sans+3:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">

    <!-- Stylesheets -->
    <link rel="stylesheet" href="CSS/design-system.css">
    <link rel="stylesheet" href="CSS/header.css">
    <link rel="stylesheet" href="CSS/details.css">
    <link rel="stylesheet" href="CSS/animations.css">
    <link rel="stylesheet" href="CSS/construction-modal.css">
    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">

    <!-- Scripts -->
    <script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
    <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/xlsx/0.18.5/xlsx.full.min.js"></script>
    <script src="JS/micro-interactions.js"></script>
</head>
<body style="background: #faf8f5;" class="content-fade-in">

<!-- Header -->
<header class="site-header">
    <div class="container">
        <a href="index.jsp" class="site-logo">scSAID</a>
        <nav class="main-nav">
            <a href="index.jsp" class="main-nav__link">Home</a>
            <a href="browse.jsp" class="main-nav__link">Browse</a>
            <a href="search.jsp" class="main-nav__link">Search</a>
            <a href="gene-search.jsp" class="main-nav__link">Gene Search</a>
            <a href="download.jsp" class="main-nav__link">Download</a>
        </nav>
        <a href="https://zje.zju.edu.cn/zje/main.htm" target="_blank">
            <img src="images/ZJE_Logo.png"
                 alt="ZJE - Zhejiang University" class="university-logo">
        </a>
    </div>
</header>
<div class="details-box">
    <!-- Sidebar Navigation -->
    <aside class="sidebar">
        <h1 class="sidebar__title">Dataset Navigation</h1>
        <nav class="sidebar__nav">
            <a href="#ExperimentInformation" class="nav-item active">
                <svg class="nav-item__icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                    <circle cx="12" cy="12" r="10"></circle>
                    <path d="M12 16v-4M12 8h.01"></path>
                </svg>
                General Information
            </a>
            <a href="#CellClustering" class="nav-item">
                <svg class="nav-item__icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                    <circle cx="12" cy="12" r="3"></circle>
                    <circle cx="19" cy="12" r="2"></circle>
                    <circle cx="5" cy="12" r="2"></circle>
                    <circle cx="12" cy="5" r="2"></circle>
                    <circle cx="12" cy="19" r="2"></circle>
                </svg>
                Cell Clustering
            </a>
            <a href="#DEGResults" class="nav-item">
                <svg class="nav-item__icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                    <path d="M3 3v18h18"></path>
                    <path d="M18.7 8l-5.1 5.2-2.8-2.7L7 14.3"></path>
                </svg>
                DEG Results
            </a>
            <a href="#CellPhoneDBAnalysis" class="nav-item">
                <svg class="nav-item__icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                    <path d="M17 21v-2a4 4 0 0 0-4-4H5a4 4 0 0 0-4 4v2"></path>
                    <circle cx="9" cy="7" r="4"></circle>
                    <path d="M23 21v-2a4 4 0 0 0-3-3.87"></path>
                    <path d="M16 3.13a4 4 0 0 1 0 7.75"></path>
                </svg>
                CellPhoneDB Analysis
            </a>
            <a href="#EnrichmentAnalysis" class="nav-item">
                <svg class="nav-item__icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                    <path d="M2 3h6a4 4 0 0 1 4 4v14a3 3 0 0 0-3-3H2z"></path>
                    <path d="M22 3h-6a4 4 0 0 0-4 4v14a3 3 0 0 1 3-3h7z"></path>
                </svg>
                Enrichment Analysis
            </a>
            <a href="#GeneSetScoring" class="nav-item">
                <svg class="nav-item__icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                    <path d="M12 2L2 7l10 5 10-5-10-5z"></path>
                    <path d="M2 17l10 5 10-5"></path>
                    <path d="M2 12l10 5 10-5"></path>
                </svg>
                Gene Set Scoring
            </a>
            <a href="visualization?dataset=<%= saidVal %>" class="nav-item nav-item--highlight">
                <svg class="nav-item__icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                    <rect x="3" y="3" width="18" height="18" rx="2" ry="2"></rect>
                    <circle cx="8.5" cy="8.5" r="1.5"></circle>
                    <polyline points="21 15 16 10 5 21"></polyline>
                </svg>
                Interactive Visualizations
                <span class="nav-badge">NEW</span>
            </a>
        </nav>
    </aside>
    <div class="basic" id="ExperimentInformation" data-panel-enter>
        <div class="general_info" style="height: 600px;">
            <div class="header">General Information</div>
            <div class="general_info_part">
                <div style="width: 40%">
                    <div class="title_1">Overview<div class="separator"></div></div>

                    <div class="detail_container_1"><div class="subtitle">Data ID: </div><div class="text_2"><%= saidVal %></div></div>
                    <div class="detail_container_1"><div class="subtitle">GSE: </div><div class="text_2"><% if (gseVal != null && !gseVal.trim().isEmpty()) { %><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<%= gseVal %>" target="_blank" rel="noopener noreferrer" class="geo-link"><%= gseVal %></a><% } else { %><%= gseVal %><% } %></div></div>
                    <div class="detail_container_1"><div class="subtitle">GSM: </div><div class="text_2"><% if (gsmVal != null && !gsmVal.trim().isEmpty()) { %><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<%= gsmVal %>" target="_blank" rel="noopener noreferrer" class="geo-link"><%= gsmVal %></a><% } else { %><%= gsmVal %><% } %></div></div>
                    <div class="detail_container_1"><div class="subtitle">Species: </div><div class="text_2"><%= speciesVal %></div></div>
                    <div class="detail_container_1"><div class="subtitle">Disease: </div><div class="text_2"><%= diseaseVal %></div></div>

                    <div class="detail_container_1"><div class="subtitle">Tissue: </div><div class="text_2"><%= tissueVal %></div></div>

                    <div class="title_1" style="margin-top: 20px;">Characteristics<div class="separator"></div></div>
                    <% String[] charItems = charVal != null ? charVal.split("\\s*;\\s*") : new String[0]; for (String item : charItems) { if (item == null || item.trim().isEmpty()) continue; %>
                    <div class="detail_container_1"><div class="subtitle"></div><div class="text_2"><%= item.trim() %></div></div>
                    <% } %>
                </div>
                <div style="width: 60%">
                    <div class="title_1">Experiment Information<div class="separator"></div></div>

                    <div style="max-height: 500px; overflow-y: auto; padding-right: 10px;">
                        <div class="detail_container_2"><div class="subtitle">Title: </div><div class="text_2"><%= titleVal %></div></div>
                        <div class="detail_container_2"><div class="subtitle">Summary: </div><div class="text_2"><%= summaryVal %></div></div>
                        <div class="detail_container_2"><div class="subtitle">Overall Design: </div><div class="text_2"><%= designVal %></div></div>

                    </div>
                </div>
            </div>
        </div>
        <div class="CellClustering" id="CellClustering" data-panel-enter>
            <div class="cluster">
                <div class="header">Cell Clustering</div>
                <div id="dash-container" style="width:1000px; height:800px;">
                    <iframe
                            src="/dash/?saids=<%= java.net.URLEncoder.encode(saidVal, "UTF-8") %>"
                            style="width:100%; height:100%; border:0;"
                    ></iframe>
                </div>
            </div>

            <div class="cluster" id="DEGResults" data-panel-enter>
                <div class="header">
                    <div class="header-content">
                        <div>
                            <div class="header-title">Differentially Expressed Genes</div>
                        </div>
                        <button id="exportExcelBtn" class="export-btn" data-btn-morph>
                            <svg class="export-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"></path>
                                <polyline points="7 10 12 15 17 10"></polyline>
                                <line x1="12" y1="15" x2="12" y2="3"></line>
                            </svg>
                            Export Excel
                        </button>
                    </div>
                </div>
                <div class="panel-body">
                    <div class="deg-controls">
                        <div class="filter-grid">
                            <div class="filter-card">
                                <div class="filter-label">
                                    <span class="filter-name">p-value threshold</span>
                                    <span class="filter-value" id="pvalLabel">0.05</span>
                                </div>
                                <input type="range" id="pvalSlider" class="elegant-slider" min="0" max="0.1" step="0.001" value="0.05">
                                <div class="filter-hint">Maximum adjusted p-value</div>
                            </div>
                            <div class="filter-card">
                                <div class="filter-label">
                                    <span class="filter-name">Log fold change</span>
                                    <span class="filter-value" id="fcLabel">1.0</span>
                                </div>
                                <input type="range" id="fcSlider" class="elegant-slider" min="0" max="10" step="0.1" value="1.0">
                                <div class="filter-hint">Minimum log₂ fold change</div>
                            </div>
                            <div class="filter-card">
                                <div class="filter-label">
                                    <span class="filter-name">Cell type group</span>
                                </div>
                                <select id="groupSelect" class="elegant-select">
                                    <option value="">All groups</option>
                                </select>
                                <div class="filter-hint">Filter by cell type</div>
                            </div>
                        </div>
                    </div>
                    <div class="table-wrapper">
                        <table id="degTable" class="elegant-table" style="width:100%">
                            <thead>
                                <tr>
                                    <th>Gene</th>
                                    <th>logFC</th>
                                    <th>p-value</th>
                                    <th>Score</th>
                                    <th>Group</th>
                                </tr>
                            </thead>
                            <tbody></tbody>
                        </table>
                    </div>
                </div>
            </div>

            <div class="cluster" id="CellPhoneDBAnalysis" data-panel-enter>
                <div class="header">
                    <div class="header-content">
                        <div>
                            <div class="header-title">CellPhoneDB Cell-Cell Communication Inference</div>
                        </div>
                    </div>
                </div>
                <div class="panel-body">
                    <!-- Summary Plot Section -->
                    <div class="cpdb-section">
                        <div class="cpdb-section-header">
                            <div class="cpdb-section-info">
                                <h3 class="cpdb-section-title">Interaction Overview</h3>
                                <p class="cpdb-section-desc">Heatmap of significant ligand-receptor interactions across cell types</p>
                            </div>
                            <button id="cpdbSummaryBtn" class="generate-btn" data-btn-morph>
                                <svg class="btn-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <polyline points="23 4 23 10 17 10"></polyline>
                                    <path d="M20.49 15a9 9 0 1 1-2.12-9.36L23 10"></path>
                                </svg>
                                Generate Plot
                            </button>
                        </div>
                        <div id="cpdbSummaryPlotContainer" class="plot-container">
                            <div class="plot-placeholder">
                                <svg class="placeholder-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5">
                                    <rect x="3" y="3" width="7" height="7"></rect>
                                    <rect x="14" y="3" width="7" height="7"></rect>
                                    <rect x="14" y="14" width="7" height="7"></rect>
                                    <rect x="3" y="14" width="7" height="7"></rect>
                                </svg>
                                <p class="placeholder-text">Click "Generate Plot" to visualize interaction summary</p>
                            </div>
                        </div>
                    </div>

                    <div class="section-divider"></div>

                    <!-- Receiver Plot Section -->
                    <div class="cpdb-section">
                        <div class="cpdb-section-header">
                            <div class="cpdb-section-info">
                                <h3 class="cpdb-section-title">Top Receptor Interactions</h3>
                                <p class="cpdb-section-desc">Most significant interactions for a specific receiving cell type</p>
                            </div>
                        </div>
                        <div class="cpdb-receiver-controls">
                            <div class="receiver-select-wrapper">
                                <label for="cpdbCellTypeSelect" class="receiver-label">Receiving cell type</label>
                                <select id="cpdbCellTypeSelect" class="elegant-select receiver-select">
                                    <option value="">Loading cell types...</option>
                                </select>
                            </div>
                            <button id="cpdbReceiverBtn" class="generate-btn" data-btn-morph>
                                <svg class="btn-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <polyline points="23 4 23 10 17 10"></polyline>
                                    <path d="M20.49 15a9 9 0 1 1-2.12-9.36L23 10"></path>
                                </svg>
                                Generate Plot
                            </button>
                        </div>
                        <div id="cpdbReceiverPlotContainer" class="plot-container">
                            <div class="plot-placeholder">
                                <svg class="placeholder-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5">
                                    <circle cx="12" cy="12" r="3"></circle>
                                    <path d="M12 1v6m0 6v6M1 12h6m6 0h6"></path>
                                    <circle cx="12" cy="12" r="10"></circle>
                                </svg>
                                <p class="placeholder-text">Select a cell type and generate to view top interactions</p>
                            </div>
                        </div>
                    </div>
                </div>
            </div>

            <div class="cluster" id="EnrichmentAnalysis" data-panel-enter>
                <div class="header">
                    <div class="header-content">
                        <div>
                            <div class="header-title">Cell Type Enrichment Analysis</div>
                            <div class="header-subtitle">Gene set enrichment using MSigDB databases</div>
                        </div>
                    </div>
                </div>
                <div class="panel-body">
                    <div class="enrichment-iframe-container">
                        <iframe
                            id="enrichmentFrame"
                            src="http://localhost:8051/enrichment/"
                            style="width: 100%; height: 1200px; border: none; border-radius: 8px;"
                            title="Enrichment Analysis">
                        </iframe>
                    </div>
                </div>
            </div>

            <div class="cluster" id="GeneSetScoring" data-panel-enter>
                <div class="header">
                    <div class="header-content">
                        <div>
                            <div class="header-title">Gene Set Scoring Analysis</div>
                            <div class="header-subtitle">Score cell types using Reactome pathways with AUCell, GSVA, or Mean Expression</div>
                        </div>
                    </div>
                </div>
                <div class="panel-body">
                    <div class="scoring-iframe-container">
                        <iframe
                            id="scoringFrame"
                            src="http://localhost:8052/gene-scoring/?dataset=<%= saidVal %>"
                            style="width: 100%; height: 1400px; border: none; border-radius: 8px;"
                            title="Gene Set Scoring">
                        </iframe>
                    </div>
                </div>
            </div>
        </div>

        <script>
            $(document).ready(function() {
                document.querySelectorAll('a[href^="#"]').forEach(a => {
                    a.addEventListener('click', function (e) {
                        e.preventDefault();
                        const target = document.querySelector(this.getAttribute('href'));
                        const top = target.getBoundingClientRect().top + window.pageYOffset - 60 - 15;
                        window.scrollTo({ top: top, behavior: 'smooth' });
                    });
                });
                const offset = 120;
                const links  = document.querySelectorAll('.nav-item[href^="#"]');
                const sections = Array.from(links, a => document.querySelector(a.getAttribute('href')))
                    .filter(el => el);

                function highlight() {
                    let curr = '';
                    sections.forEach(sec => {
                        const rect = sec.getBoundingClientRect();
                        if (rect.top <= offset && rect.bottom > offset) curr = sec.id;
                    });
                    links.forEach(link => link.classList.toggle('active', link.getAttribute('href') === '#' + curr));
                }

                window.addEventListener('scroll', highlight);
                highlight();

            })
            $(function(){

                // =========================================================================
                // Original DEG Script
                // =========================================================================
                console.log("🌟 Original script start");
                const table = $('#degTable').DataTable({ paging:true, searching:false, info:true });

                const said = '<%= saidParam %>';
                const gse = '<%= gseVal %>';
                const gsm = '<%= gsmVal %>';

                function initGroupOptions() {
                    console.log("🔍 Initializing DEG group dropdown");

                    $.getJSON('/scrna_website_test_war/deg', { said: said, pval: 1.0, fc: 0.0 })
                        .done(function(data){
                            console.log("✅ DEG group data fetched");
                            const groups = Array.from(new Set(data.map(r => typeof r.group === "string" ? r.group.trim() : null).filter(g => g)));
                            const select = $('#groupSelect');
                            select.empty().append('<option value="">All</option>');
                            groups.forEach(g => select.append('<option value="' + g + '">' + g + '</option>'));
                            console.log("✅ DEG group dropdown populated");
                        })
                        .fail(function(xhr){ console.error("❌ DEG group data failed:", xhr.status, xhr.statusText); });
                }

                function loadDEG(){
                    const pval = $('#pvalSlider').val();
                    const fc = $('#fcSlider').val();
                    const group = $('#groupSelect').val();
                    $('#pvalLabel').text(pval);
                    $('#fcLabel').text(fc);
                    const params = { said: said, pval: pval, fc: fc };
                    if (group) params.group = group;
                    console.log("📡 Requesting DEG data:", params);
                    $.getJSON('/scrna_website_test_war/deg', params)
                        .done(function(data){
                            console.log("✅ DEG data received");
                            table.clear();

                            data.forEach(r => table.row.add([r.gene, r.logfoldchanges, r.pvals_adj, r.scores, r.group]));
                            table.draw();

                            // Animate table rows after draw
                            setTimeout(() => {
                                MicroInteractions.animateTableRows(document.getElementById('degTable'));
                            }, 100);
                        })
                        .fail(function(xhr){ console.error("❌ DEG data loading failed:", xhr.status, xhr.statusText); });
                }

                function exportTableToExcel() {
                    const btn = document.getElementById('exportExcelBtn');
                    MicroInteractions.setButtonState(btn, 'processing');

                    setTimeout(() => {
                        try {
                            const exportData = [["Gene", "logFC", "p-value", "Score", "Group"]];
                            table.rows({ search: 'applied' }).every(function () { exportData.push(this.data()); });
                            const ws = XLSX.utils.aoa_to_sheet(exportData);
                            const wb = XLSX.utils.book_new();
                            XLSX.utils.book_append_sheet(wb, ws, "Filtered_DEG");
                            XLSX.writeFile(wb, "filtered_DEG_results.xlsx");

                            MicroInteractions.setButtonState(btn, 'success');
                            MicroInteractions.showFeedback('success', 'Excel file exported successfully!');
                        } catch (error) {
                            MicroInteractions.setButtonState(btn, 'error');
                            MicroInteractions.showFeedback('error', 'Export failed. Please try again.');
                        }
                    }, 300);
                }

                initGroupOptions();
                loadDEG();
                $('#pvalSlider, #fcSlider, #groupSelect').on('input change', loadDEG);
                $('#exportExcelBtn').on('click', exportTableToExcel);

                // =========================================================================
                // SECTION D: NEW JAVASCRIPT FOR CELLPHONEDB
                // =========================================================================
                console.log("🌟 CellPhoneDB script start");
                function initCpdb() {
                    console.log("🔍 Initializing CellPhoneDB cell type dropdown");
                    // Use the gse and gsm values fetched by the main JSP logic
                    $.getJSON('', { action: 'get_cell_types', said: said, gse: gse, gsm: gsm })
                        .done(function(data) {
                            if (data.error) {
                                console.error("❌ CPDB Error:", data.error);
                                $('#cpdbCellTypeSelect').empty().append('<option value="">Error loading types</option>');
                                alert("Could not load CellPhoneDB cell types: " + data.error);
                                return;
                            }
                            console.log("✅ CPDB cell types fetched:", data.cell_types);
                            const select = $('#cpdbCellTypeSelect');
                            select.empty().append('<option value="">-- Select a Cell Type --</option>');
                            data.cell_types.forEach(ct => {
                                const optionHtml = '<option value="' + ct + '">' + ct + '</option>';
                                select.append(optionHtml);
                            });
                            console.log("✅ CPDB cell type dropdown populated");
                        })
                        .fail(function(xhr) {
                            console.error("❌ CPDB cell type request failed:", xhr.status, xhr.statusText);
                            $('#cpdbCellTypeSelect').empty().append('<option value="">Request failed</option>');
                            alert('Server error while loading Cell Types: ' + xhr.status + ' ' + xhr.statusText);
                        });
                }

                function generateCpdbPlot(plotType, cellType, containerId) {
                    // 修改：确保 id 匹配
                    const container = $('#' + containerId);
                    container.html('<div class="spinner-precision"></div>'); // Show loader

                    const btnId = plotType === 'summary' ? 'cpdbSummaryBtn' : 'cpdbReceiverBtn';
                    const btn = document.getElementById(btnId);
                    MicroInteractions.setButtonState(btn, 'processing');

                    const params = {
                        action: 'generate_plot',
                        plot_type: plotType,
                        said: said,
                        gse: gse,
                        gsm: gsm
                    };
                    if (cellType) {
                        params.cell_type = cellType;
                    }

                    console.log(`📡 Requesting CPDB plot '${plotType}':`, params);
                    $.getJSON('', params)
                        .done(function(data) {
                            if (data.error) {
                                console.error(`❌ CPDB plot generation failed: ${data.error}`);
                                container.html('<p style="color:red;"><strong>Error:</strong> ' + data.error + '</p>');
                                MicroInteractions.setButtonState(btn, 'error');
                                MicroInteractions.showFeedback('error', 'Plot generation failed: ' + data.error);
                                return;
                            }

                            if (data.imageUrl) {
                                console.log(`✅ CPDB plot generated. URL: ${data.imageUrl}`);
                                // 修改：使用传统字符串拼接以提高兼容性，并添加日志
                                const imgHtml = '<img src="' + data.imageUrl + '" alt="Generated ' + plotType + ' plot" class="chart-enter">';
                                console.log('Generated image HTML:', imgHtml);
                                container.html(imgHtml);
                                MicroInteractions.setButtonState(btn, 'success');
                                MicroInteractions.showFeedback('success', 'Plot generated successfully!');
                            } else {
                                console.error("❌ No imageUrl found in server response:", data);
                                container.html('<p style="color:red;"><strong>Error:</strong> Server response did not contain an image URL.</p>');
                                MicroInteractions.setButtonState(btn, 'error');
                                MicroInteractions.showFeedback('error', 'Plot generation failed.');
                            }
                        })
                        .fail(function(xhr) {
                            console.error(`❌ CPDB plot request failed: ${xhr.status} ${xhr.statusText}`);
                            container.html('<p style="color:red;"><strong>Request Failed:</strong> Server returned an error.</p>');
                            MicroInteractions.setButtonState(btn, 'error');
                            MicroInteractions.showFeedback('error', 'Server error. Please try again.');
                        });
                }

                // Bind events
                $('#cpdbSummaryBtn').on('click', function() {
                    generateCpdbPlot('summary', null, 'cpdbSummaryPlotContainer');
                });
                $('#cpdbReceiverBtn').on('click', function() {
                    const selectedCellType = $('#cpdbCellTypeSelect').val();
                    if (!selectedCellType) {
                        alert("Please select a cell type first.");
                        return;
                    }
                    generateCpdbPlot('receiver', selectedCellType, 'cpdbReceiverPlotContainer');
                });
                // Initialize the CellPhoneDB section
                initCpdb();
            });
        </script>
    </div>
</div>

<!-- Under Construction Modal Script -->
<script src="JS/construction-modal.js"></script>

</body>
</html>
