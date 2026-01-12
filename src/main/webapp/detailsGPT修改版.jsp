<%@ page contentType="text/html;charset=UTF-8" language="java" %>
<%@ page import="java.io.*" %>
<%@ page import="java.util.*" %>
<%@ page import="org.json.JSONObject" %>

<%!
    // ✅ Python 路径（全局常量）
    private static final String PYTHON_COMMAND = "C:\\\\Anaconda3\\\\envs\\\\dash_env\\\\python.exe";

    // ✅ 工具函数：执行 Python
    public String runPython(String script, String arg, String basePath) throws Exception {
        ProcessBuilder pb;
        if (arg != null && !arg.isEmpty()) {
            pb = new ProcessBuilder(PYTHON_COMMAND, script, arg);
        } else {
            pb = new ProcessBuilder(PYTHON_COMMAND, script);
        }
        pb.directory(new File(basePath));
        pb.redirectErrorStream(true);

        Process process = pb.start();
        BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
        StringBuilder output = new StringBuilder();
        String line;
        while ((line = reader.readLine()) != null) {
            output.append(line).append("\n");
        }
        int exitCode = process.waitFor();
        if (exitCode != 0) {
            throw new RuntimeException("❌ Python 脚本执行失败 (exitCode=" + exitCode + ")\n\n" + output);
        }
        return output.toString();
    }
%>

<%
    String gse = request.getParameter("gse");
    String gsm = request.getParameter("gsm");

    if (gse == null || gsm == null || gse.isEmpty() || gsm.isEmpty()) {
        response.setContentType("application/json;charset=UTF-8");
        System.out.println("{\"error\": \"Missing gse or gsm parameter.\"}");
        return;
    }

    // cpdb_out 基础目录
    String basePath = application.getRealPath("/") + "SkinDB_New/10X/human/" + gse + "/" + gsm + "/cpdb_out";

    // Debug JSON 输出
    JSONObject debugInfo = new JSONObject();
    debugInfo.put("gse", gse);
    debugInfo.put("gsm", gsm);
    debugInfo.put("cpdbPath", basePath);
%>

<html>
<head>
    <title>Details Page</title>
</head>
<body>
<h2>DEG 结果展示</h2>
<div id="deg-section">
    <!-- 这里放你原本的 DEG 展示逻辑 -->
</div>

<hr>

<!-- ✅ Summary Plot 区域 -->
<h2>Generate a summary plot of significant interactions</h2>
<%
    try {
        String scriptPath = basePath + File.separator + "plot_cpdb_sum_sig.py";
        runPython(scriptPath, null, basePath);

        String imgPath = "SkinDB_New/10X/human/" + gse + "/" + gsm + "/cpdb_out/summary_plot.png";
%>
<div>
    <img src="<%= imgPath %>" alt="Summary Plot" style="max-width:800px; border:1px solid #ccc;">
</div>
<%
    } catch (Exception e) {
        System.out.println("<pre style='color:red'>" + e.getMessage() + "</pre>");
    }
%>

<hr>

<!-- ✅ Receiver Plot 区域 -->
<h2>Generate Receiver Plot</h2>
<form method="get" action="details.jsp">
    <input type="hidden" name="gse" value="<%= gse %>">
    <input type="hidden" name="gsm" value="<%= gsm %>">

    <label for="celltype">选择细胞类型：</label>
    <select name="celltype" id="celltype">
        <%
            try {
                File cellFile = new File(basePath, "celltypes.txt");
                if (cellFile.exists()) {
                    BufferedReader br = new BufferedReader(new FileReader(cellFile));
                    String line;
                    while ((line = br.readLine()) != null) {
        %>
        <option value="<%= line.trim() %>"><%= line.trim() %></option>
        <%
                    }
                    br.close();
                } else {
                    System.out.println("<option disabled>celltypes.txt 不存在</option>");
                }
            } catch (Exception e) {
                System.out.println("<option disabled>读取 celltypes.txt 出错</option>");
            }
        %>
    </select>
    <button type="submit">生成图</button>
</form>

<%
    String selectedCell = request.getParameter("celltype");
    if (selectedCell != null && !selectedCell.isEmpty()) {
        try {
            String scriptPath = basePath + File.separator + "plot_cpdb_receiver_top15.py";
            runPython(scriptPath, selectedCell, basePath);

            String imgPath = "SkinDB_New/10X/human/" + gse + "/" + gsm + "/cpdb_out/receiver_plot.png";
%>
<div>
    <img src="<%= imgPath %>" alt="Receiver Plot" style="max-width:800px; border:1px solid #ccc;">
</div>
<%
        } catch (Exception e) {
            System.out.println("<pre style='color:red'>" + e.getMessage() + "</pre>");
        }
    }
%>

<hr>
<h3>调试信息</h3>
<pre><%= debugInfo.toString(2) %></pre>
</body>
</html>
