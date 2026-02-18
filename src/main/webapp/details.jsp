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
    <link href="https://fonts.googleapis.com/css2?family=Cormorant+Garamond:wght@400;500;600;700&family=Montserrat:wght@200;300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">

    <!-- Stylesheets -->
    <link rel="stylesheet" href="CSS/design-system.css">
    <link rel="stylesheet" href="CSS/header.css">
    <link rel="stylesheet" href="CSS/details.css">
    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">

    <!-- Scripts -->
    <script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
    <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/xlsx/0.18.5/xlsx.full.min.js"></script>
    <script src="https://cdn.plot.ly/plotly-2.20.0.min.js"></script>
</head>
<body style="background: #faf8f5;">

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
    <div class="basic"  id="ExperimentInformation">
        <div class="general_info" style="height: 600px;">
            <div class="header">General Information</div>
            <div class="general_info_part">
                <div style="width: 40%">
                    <div class="title_1">Overview<div class="separator"></div></div>

                    <div class="detail_container_1"><div class="subtitle">Data ID: </div><div class="text_2"><%= saidVal %></div></div>
                    <div class="detail_container_1"><div class="subtitle">GSE: </div><div class="text_2"><%= gseVal %></div></div>
                    <div class="detail_container_1"><div class="subtitle">GSM: </div><div class="text_2"><%= gsmVal %></div></div>
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
        <div class="CellClustering" id="CellClustering">
            <div class="cluster">
                <div class="header">Cell Clustering</div>
                <div id="dash-container" style="width:1000px; height:800px;">
                    <iframe src="http://localhost:8050/dash/?sample_id=${sid}" style="width:100%; height:100%; border:0;" scrolling="no"></iframe>
                </div>
            </div>

            <div class="cluster" id="DEGResults">
                <div class="header">
                    <div class="header-content">
                        <div>
                            <div class="header-title">Differentially Expressed Genes</div>
                        </div>
                        <button id="exportExcelBtn" class="export-btn">
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

            <div class="cluster" id="CellPhoneDBAnalysis">
                <div class="header">
                    <div class="header-content">
                        <div>
                            <div class="header-title">CellPhoneDB Cell-Cell Communication Analysis</div>
                        </div>
                        <span class="cpdb-badge">Dynamic Analysis</span>
                    </div>
                </div>
                <div class="panel-body">
                    <!-- Cell Type Selection -->
                    <div class="cpdb-config-section">
                        <h3 class="cpdb-section-title">Select Cell Types</h3>
                        <p class="cpdb-section-desc">Choose 2 or more cell types to analyze ligand-receptor interactions</p>

                        <div class="cpdb-cell-selector">
                            <select id="cpdbCellTypeMultiSelect" multiple class="cpdb-multiselect">
                                <option value="">Loading cell types...</option>
                            </select>
                        </div>

                        <div class="cpdb-mode-toggle">
                            <label class="cpdb-radio">
                                <input type="radio" name="cpdbMode" value="all" checked>
                                <span>All Combinations</span>
                            </label>
                            <label class="cpdb-radio">
                                <input type="radio" name="cpdbMode" value="directed">
                                <span>Sender → Receiver</span>
                            </label>
                        </div>

                        <!-- Sender/Receiver (shown when directed mode selected) -->
                        <div id="cpdbDirectedControls" style="display:none;">
                            <div class="cpdb-directed-row">
                                <div class="cpdb-directed-col">
                                    <label>Sender Cell Types</label>
                                    <select id="cpdbSenderTypes" multiple></select>
                                </div>
                                <div class="cpdb-directed-arrow">→</div>
                                <div class="cpdb-directed-col">
                                    <label>Receiver Cell Types</label>
                                    <select id="cpdbReceiverTypes" multiple></select>
                                </div>
                            </div>
                        </div>

                        <button id="runCpdbAnalysisBtn" class="generate-btn">
                            <svg class="btn-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                <polygon points="5 3 19 12 5 21 5 3"></polygon>
                            </svg>
                            Run Analysis
                        </button>
                    </div>

                    <!-- Progress Section -->
                    <div id="cpdbProgressSection" class="cpdb-progress" style="display:none;">
                        <div class="cpdb-progress-bar">
                            <div class="cpdb-progress-fill"></div>
                        </div>
                        <p class="cpdb-progress-text">Running CellPhoneDB analysis...</p>
                        <p class="cpdb-progress-hint">This may take several minutes for large datasets</p>
                    </div>

                    <!-- Results Section -->
                    <div id="cpdbResultsSection" style="display:none;">
                        <div class="cpdb-results-tabs">
                            <button class="cpdb-tab active" data-tab="heatmap">Interaction Heatmap</button>
                            <button class="cpdb-tab" data-tab="dotplot">Dot Plot</button>
                            <button class="cpdb-tab" data-tab="table">Results Table</button>
                        </div>

                        <div id="cpdbHeatmapTab" class="cpdb-tab-content active">
                            <div id="cpdbHeatmapPlot" class="cpdb-plot-container"></div>
                        </div>

                        <div id="cpdbDotplotTab" class="cpdb-tab-content">
                            <div id="cpdbDotplot" class="cpdb-plot-container"></div>
                        </div>

                        <div id="cpdbTableTab" class="cpdb-tab-content">
                            <div class="table-wrapper">
                                <table id="cpdbResultsTable" class="elegant-table">
                                    <thead>
                                        <tr>
                                            <th>Interaction</th>
                                            <th>Ligand</th>
                                            <th>Receptor</th>
                                            <th>Sender</th>
                                            <th>Receiver</th>
                                            <th>Score</th>
                                        </tr>
                                    </thead>
                                    <tbody></tbody>
                                </table>
                            </div>
                        </div>
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
                        })
                        .fail(function(xhr){ console.error("❌ DEG data loading failed:", xhr.status, xhr.statusText); });
                }

                function exportTableToExcel() {
                    const exportData = [["Gene", "logFC", "p-value", "Score", "Group"]];
                    table.rows({ search: 'applied' }).every(function () { exportData.push(this.data()); });
                    const ws = XLSX.utils.aoa_to_sheet(exportData);
                    const wb = XLSX.utils.book_new();
                    XLSX.utils.book_append_sheet(wb, ws, "Filtered_DEG");
                    XLSX.writeFile(wb, "filtered_DEG_results.xlsx");
                }

                initGroupOptions();
                loadDEG();
                $('#pvalSlider, #fcSlider, #groupSelect').on('input change', loadDEG);
                $('#exportExcelBtn').on('click', exportTableToExcel);

                // =========================================================================
                // SECTION D: CELLPHONEDB DYNAMIC ANALYSIS
                // =========================================================================
                console.log("🌟 CellPhoneDB Dynamic Analysis script start");

                // Initialize cell type multi-select
                function initCpdbCellTypes() {
                    console.log("🔍 Loading cell types for CellPhoneDB");
                    $.getJSON('/scrna_website_test_war/cpdb-api', {
                        action: 'cell-types',
                        said: said
                    }).done(function(data) {
                        if (data.error) {
                            console.error('CPDB cell types error:', data.error);
                            $('#cpdbCellTypeMultiSelect').html('<option value="">Error: ' + data.error + '</option>');
                            return;
                        }
                        console.log("✅ CPDB cell types loaded:", data.cell_types);
                        const select = $('#cpdbCellTypeMultiSelect');
                        select.empty();
                        data.cell_types.forEach(function(ct) {
                            const count = data.cell_counts ? data.cell_counts[ct] || '' : '';
                            const label = count ? ct + ' (' + count + ' cells)' : ct;
                            select.append('<option value="' + ct + '">' + label + '</option>');
                        });
                    }).fail(function(xhr) {
                        console.error("❌ CPDB cell type request failed:", xhr.status, xhr.statusText);
                        $('#cpdbCellTypeMultiSelect').html('<option value="">Failed to load cell types</option>');
                    });
                }

                // Toggle directed mode controls
                $('input[name="cpdbMode"]').change(function() {
                    const isDirected = $(this).val() === 'directed';
                    $('#cpdbDirectedControls').toggle(isDirected);

                    if (isDirected) {
                        // Copy selected types to sender/receiver
                        const selected = $('#cpdbCellTypeMultiSelect').val() || [];
                        ['#cpdbSenderTypes', '#cpdbReceiverTypes'].forEach(function(sel) {
                            $(sel).empty();
                            selected.forEach(function(ct) {
                                $(sel).append('<option value="' + ct + '">' + ct + '</option>');
                            });
                        });
                    }
                });

                // Update sender/receiver when main selection changes
                $('#cpdbCellTypeMultiSelect').change(function() {
                    const selected = $(this).val() || [];
                    if ($('input[name="cpdbMode"]:checked').val() === 'directed') {
                        ['#cpdbSenderTypes', '#cpdbReceiverTypes'].forEach(function(sel) {
                            const current = $(sel).val() || [];
                            $(sel).empty();
                            selected.forEach(function(ct) {
                                const isSelected = current.indexOf(ct) !== -1;
                                $(sel).append('<option value="' + ct + '"' + (isSelected ? ' selected' : '') + '>' + ct + '</option>');
                            });
                        });
                    }
                });

                // Run analysis
                $('#runCpdbAnalysisBtn').click(function() {
                    const selectedTypes = $('#cpdbCellTypeMultiSelect').val();

                    if (!selectedTypes || selectedTypes.length < 2) {
                        alert('Please select at least 2 cell types');
                        return;
                    }

                    const mode = $('input[name="cpdbMode"]:checked').val();
                    const params = {
                        action: 'run-analysis',
                        said: said,
                        cell_types: JSON.stringify(selectedTypes)
                    };

                    if (mode === 'directed') {
                        const senders = $('#cpdbSenderTypes').val();
                        const receivers = $('#cpdbReceiverTypes').val();
                        if (senders && senders.length > 0) {
                            params.senders = JSON.stringify(senders);
                        }
                        if (receivers && receivers.length > 0) {
                            params.receivers = JSON.stringify(receivers);
                        }
                    }

                    // Show progress
                    $('#cpdbProgressSection').show();
                    $('#cpdbResultsSection').hide();
                    $('.cpdb-progress-fill').css('width', '10%');

                    console.log("📡 Starting CPDB analysis:", params);

                    $.post('/scrna_website_test_war/cpdb-api', params)
                        .done(function(response) {
                            console.log("✅ CPDB analysis started:", response);
                            if (response.job_id) {
                                pollCpdbStatus(response.job_id);
                            } else if (response.error) {
                                showCpdbError(response.error);
                            }
                        })
                        .fail(function(xhr) {
                            console.error("❌ CPDB analysis request failed:", xhr.status, xhr.statusText);
                            showCpdbError('Request failed: ' + xhr.statusText);
                        });
                });

                // Poll job status
                function pollCpdbStatus(jobId) {
                    const poll = setInterval(function() {
                        $.getJSON('/scrna_website_test_war/cpdb-api', {
                            action: 'status',
                            job_id: jobId
                        }).done(function(data) {
                            console.log("📊 CPDB job status:", data);
                            if (data.status === 'completed') {
                                clearInterval(poll);
                                loadCpdbResults(jobId);
                            } else if (data.status === 'failed') {
                                clearInterval(poll);
                                showCpdbError(data.error || 'Analysis failed');
                            } else {
                                // Update progress bar if available
                                if (data.progress) {
                                    $('.cpdb-progress-fill').css('width', data.progress + '%');
                                }
                            }
                        }).fail(function() {
                            clearInterval(poll);
                            showCpdbError('Failed to check job status');
                        });
                    }, 3000);
                }

                // Load and display results
                function loadCpdbResults(jobId) {
                    $.getJSON('/scrna_website_test_war/cpdb-api', {
                        action: 'results',
                        job_id: jobId
                    }).done(function(data) {
                        console.log("✅ CPDB results loaded:", data);
                        $('#cpdbProgressSection').hide();
                        $('#cpdbResultsSection').show();

                        // Render heatmap
                        renderCpdbHeatmap(data.heatmap_data);

                        // Render dot plot
                        renderCpdbDotplot(data.dotplot_data);

                        // Populate table
                        populateCpdbTable(data.interactions);
                    }).fail(function(xhr) {
                        console.error("❌ Failed to load CPDB results:", xhr.status);
                        showCpdbError('Failed to load results');
                    });
                }

                // Render heatmap using Plotly
                function renderCpdbHeatmap(heatmapData) {
                    if (!heatmapData || !heatmapData.z || heatmapData.z.length === 0) {
                        $('#cpdbHeatmapPlot').html('<div class="cpdb-error"><p class="cpdb-error-text">No interaction data to display</p></div>');
                        return;
                    }

                    const trace = {
                        z: heatmapData.z,
                        x: heatmapData.x,
                        y: heatmapData.y,
                        type: 'heatmap',
                        colorscale: [
                            [0, '#faf8f5'],
                            [0.5, '#e8927c'],
                            [1, '#8B0000']
                        ],
                        hoverongaps: false
                    };

                    const layout = {
                        title: 'Ligand-Receptor Interaction Scores',
                        font: { family: 'Montserrat, sans-serif' },
                        xaxis: {
                            title: 'Cell Type Pairs (Sender|Receiver)',
                            tickangle: -45,
                            tickfont: { size: 10 }
                        },
                        yaxis: {
                            title: 'Interactions',
                            tickfont: { size: 10 }
                        },
                        margin: { l: 150, r: 50, t: 80, b: 150 }
                    };

                    Plotly.newPlot('cpdbHeatmapPlot', [trace], layout, { responsive: true });
                }

                // Render dot plot using Plotly
                function renderCpdbDotplot(dotplotData) {
                    if (!dotplotData || !dotplotData.interactions || dotplotData.interactions.length === 0) {
                        $('#cpdbDotplot').html('<div class="cpdb-error"><p class="cpdb-error-text">No interaction data to display</p></div>');
                        return;
                    }

                    const trace = {
                        x: dotplotData.cell_pairs,
                        y: dotplotData.interactions,
                        mode: 'markers',
                        marker: {
                            size: dotplotData.sizes.map(function(s) { return Math.min(30, Math.max(5, s * 10)); }),
                            color: dotplotData.scores,
                            colorscale: 'RdBu',
                            reversescale: true,
                            showscale: true,
                            colorbar: { title: 'Score' }
                        },
                        type: 'scatter',
                        text: dotplotData.scores.map(function(s) { return 'Score: ' + s.toFixed(3); }),
                        hoverinfo: 'text+x+y'
                    };

                    const layout = {
                        title: 'Top Ligand-Receptor Interactions',
                        font: { family: 'Montserrat, sans-serif' },
                        xaxis: {
                            title: 'Cell Type Pairs',
                            tickangle: -45,
                            tickfont: { size: 10 }
                        },
                        yaxis: {
                            title: 'Interactions',
                            tickfont: { size: 10 }
                        },
                        margin: { l: 150, r: 100, t: 80, b: 150 }
                    };

                    Plotly.newPlot('cpdbDotplot', [trace], layout, { responsive: true });
                }

                // Populate results table
                function populateCpdbTable(interactions) {
                    const tbody = $('#cpdbResultsTable tbody');
                    tbody.empty();

                    if (!interactions || interactions.length === 0) {
                        tbody.append('<tr><td colspan="6" style="text-align:center;">No interactions found</td></tr>');
                        return;
                    }

                    interactions.forEach(function(int) {
                        tbody.append(
                            '<tr>' +
                            '<td>' + int.interaction + '</td>' +
                            '<td>' + int.ligand + '</td>' +
                            '<td>' + int.receptor + '</td>' +
                            '<td>' + int.cell_type_sender + '</td>' +
                            '<td>' + int.cell_type_receiver + '</td>' +
                            '<td>' + int.score.toFixed(4) + '</td>' +
                            '</tr>'
                        );
                    });
                }

                // Show error
                function showCpdbError(message) {
                    $('#cpdbProgressSection').hide();
                    $('#cpdbResultsSection').show();
                    $('#cpdbHeatmapPlot').html(
                        '<div class="cpdb-error">' +
                        '<svg class="cpdb-error-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">' +
                        '<circle cx="12" cy="12" r="10"></circle>' +
                        '<line x1="12" y1="8" x2="12" y2="12"></line>' +
                        '<line x1="12" y1="16" x2="12.01" y2="16"></line>' +
                        '</svg>' +
                        '<p class="cpdb-error-text">' + message + '</p>' +
                        '</div>'
                    );
                }

                // Tab switching
                $('.cpdb-tab').click(function() {
                    const tab = $(this).data('tab');
                    $('.cpdb-tab').removeClass('active');
                    $(this).addClass('active');
                    $('.cpdb-tab-content').removeClass('active');
                    $('#cpdb' + tab.charAt(0).toUpperCase() + tab.slice(1) + 'Tab').addClass('active');

                    // Resize Plotly charts when tab becomes visible
                    if (tab === 'heatmap') {
                        Plotly.Plots.resize('cpdbHeatmapPlot');
                    } else if (tab === 'dotplot') {
                        Plotly.Plots.resize('cpdbDotplot');
                    }
                });

                // Initialize
                initCpdbCellTypes();
            });
        </script>
    </div>
</div>

</body>
</html>