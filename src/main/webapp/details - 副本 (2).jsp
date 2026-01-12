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
            jsonOut.print("{\"error\": \"CPDB data path not found for the given sample. Please check server configuration and sample identifiers.\"}");
            jsonOut.flush();
            return;
        }

        try {
            if ("get_cell_types".equals(action)) {
                File pvaluesFile = new File(cpdbPath, "statistical_analysis_pvalues.txt");
                if (!pvaluesFile.exists()) {
                    jsonOut.print("{\"error\": \"statistical_analysis_pvalues.txt not found.\"}");
                    return;
                }
                List<String> cellTypes = new ArrayList<>();
                try (BufferedReader reader = new BufferedReader(new FileReader(pvaluesFile))) {
                    String header = reader.readLine();
                    if (header != null) {
                        String[] columns = header.split("\\t");
                        // Cell types start from the 12th column in a standard CellPhoneDB pvalues.txt file
                        for (int i = 11; i < columns.length; i++) {
                            cellTypes.add(columns[i].trim());
                        }
                    }
                }
                // Build JSON response manually to avoid new library dependencies
                StringBuilder jsonResponse = new StringBuilder();
                jsonResponse.append("{\"cell_types\": [");
                for (int i = 0; i < cellTypes.size(); i++) {
                    jsonResponse.append("\"").append(cellTypes.get(i)).append("\"");
                    if (i < cellTypes.size() - 1) {
                        jsonResponse.append(",");
                    }
                }
                jsonResponse.append("]}");
                jsonOut.print(jsonResponse.toString());
            } else if ("generate_plot".equals(action)) {
                String plotType = request.getParameter("plot_type");
                String scriptName = "";
                List<String> command = new ArrayList<>();
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
                try (InputStream is = process.getInputStream(); InputStream es = process.getErrorStream()) {
                    while (is.read() != -1);
                    while (es.read() != -1);
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
        Map<String,Integer> colIdx = new HashMap<>();
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
    <title>Details-SR001</title>
    <link rel="stylesheet" href="CSS/header.css">
    <link rel="stylesheet" href="CSS/details.css">
    <link rel="stylesheet" href="CSS/degtest.css">
    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
    <script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
    <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/xlsx/0.18.5/xlsx.full.min.js"></script>

    <style>
        select#groupSelect, select#cpdbCellTypeSelect { min-width: 180px; font-size: 14px; padding: 6px; }
        .deg-controls, .cpdb-controls { display: flex; justify-content: space-between; align-items: stretch; margin-bottom: 10px; gap: 20px; }
        .filter-controls, .cpdb-filter-controls { display: flex; gap: 20px; align-items: center; flex: 1; }
        .export-button-wrapper { display: flex; align-items: center; justify-content: center; padding: 0 20px; min-width: 180px; }
        #exportExcelBtn, .cpdb-generate-btn { padding: 12px 24px; font-size: 16px; border: none; background-color: #007BFF; color: white; border-radius: 8px; cursor: pointer; box-shadow: 0 4px 6px rgba(0,0,0,0.1); transition: background-color 0.3s ease; width: 100%; }
        #exportExcelBtn:hover, .cpdb-generate-btn:hover { background-color: #0056b3; }
        .cpdb-plot-container { text-align: center; margin-top: 20px; padding: 10px; border: 1px solid #ddd; min-height: 400px; }
        .cpdb-plot-container img { max-width: 100%; height: auto; }
        .loader { border: 8px solid #f3f3f3; border-top: 8px solid #3498db; border-radius: 50%; width: 60px; height: 60px; animation: spin 2s linear infinite; margin: 50px auto; }
        @keyframes spin { 0% { transform: rotate(0deg); } 100% { transform: rotate(360deg); } }
    </style>
</head>
<body>
<header>
    <nav>
        <ul>
            <db_logo><a href="#">sDSSA</a></db_logo>
            <li><a href="index.jsp">Home</a></li>
            <li><a href="browse.jsp">Browse</a></li>
            <li><a href="search.jsp">Search</a><div class="top_list"><lib><a href="#">Gene</a></lib><lib><a href="#">Cell</a></lib></div></li>
            <li><a href="#">Help</a><div class="top_list"><lib><a href="#">Method</a></lib><lib><a href="#">Tutorial</a></lib></div></li>

            <li><a href="#">Download</a></li>
            <un_logo><img src="https://www.zju.edu.cn/_upload/tpl/0b/bf/3007/template3007/static/js/../../static/media/mlogo.80e02913954185729616ab6a1c6ae12d.svg" alt="" width="196" height="54"></un_logo>
        </ul>
    </nav>
</header>
<div class="basic">
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
                <div class="title_1">Characteristics<div class="separator"></div></div>
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
    <div class="basic" style="width: 1000px">
        <div class="cluster">
            <div class="header">Cell Clustering</div>
            <div id="dash-container" style="width:1000px; height:800px;">
                <iframe src="http://localhost:8050/dash/?sample_id=${sid}" style="width:100%; height:100%; border:0;" scrolling="no"></iframe>
            </div>
        </div>

        <div class="cluster">
            <div class="header">DEG Results</div>
            <div class="deg-controls">
                <div class="filter-controls">
                    <div>
                        <label for="pvalSlider">p-value ≤ <span id="pvalLabel">0.05</span></label><br>
                        <input type="range" id="pvalSlider" min="0" max="0.1" step="0.001" value="0.05">
                    </div>
                    <div>
                        <label for="fcSlider">logFC ≥ <span id="fcLabel">1.0</span></label><br>
                        <input type="range" id="fcSlider" min="0" max="10" step="0.1" value="1.0">
                    </div>
                    <div>
                        <label for="groupSelect">Group</label><br>
                        <select id="groupSelect"><option value="">All</option></select>
                    </div>
                </div>
                <div class="export-button-wrapper">
                    <button id="exportExcelBtn">Export as Excel</button>
                </div>
            </div>
            <table id="degTable" class="display" style="width:100%">
                <thead><tr><th>Gene</th><th>logFC</th><th>p-value</th><th>Score</th><th>Group</th></tr></thead>
                <tbody></tbody>

            </table>
        </div>

        <div class="cluster">
            <div class="header">CellPhoneDB Analysis</div>

            <div class="cpdb-controls">
                <div class="cpdb-filter-controls" style="justify-content: space-between;">
                    <span>Generate a summary plot of significant interactions.</span>
                    <button id="cpdbSummaryBtn" class="cpdb-generate-btn" style="width: 250px;">Generate Summary Plot</button>
                </div>
            </div>
            <div id="cpdbSummaryPlotContainer" class="cpdb-plot-container">
                Click the button above to generate the summary plot.
            </div>

            <hr style="margin: 40px 0;">

            <div class="cpdb-controls">
                <div class="cpdb-filter-controls">
                    <div>
                        <label for="cpdbCellTypeSelect">Select Cell Type (Receiver)</label><br>
                        <select id="cpdbCellTypeSelect">
                            <option value="">Loading cell types...</option>
                        </select>
                    </div>
                </div>
                <div class="export-button-wrapper">
                    <button id="cpdbReceiverBtn" class="cpdb-generate-btn">Generate Receiver Plot</button>
                </div>
            </div>
            <div id="cpdbReceiverPlotContainer" class="cpdb-plot-container">
                Select a cell type and click the button to generate the plot.
            </div>
        </div>
    </div>

    <script>
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
                            select.append(`<option value="${ct}">${ct}</option>`);
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
                container.html('<div class="loader"></div>'); // Show loader

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
                            return;
                        }

                        if (data.imageUrl) {
                            console.log(`✅ CPDB plot generated. URL: ${data.imageUrl}`);
                            // 修改：使用传统字符串拼接以提高兼容性，并添加日志
                            const imgHtml = '<img src="' + data.imageUrl + '" alt="Generated ' + plotType + ' plot">';
                            console.log('Generated image HTML:', imgHtml);
                            container.html(imgHtml);
                        } else {
                            console.error("❌ No imageUrl found in server response:", data);
                            container.html('<p style="color:red;"><strong>Error:</strong> Server response did not contain an image URL.</p>');
                        }
                    })
                    .fail(function(xhr) {
                        console.error(`❌ CPDB plot request failed: ${xhr.status} ${xhr.statusText}`);
                        container.html('<p style="color:red;"><strong>Request Failed:</strong> Server returned an error.</p>');
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
</body>
</html>