<%@ page language="java"
         contentType="text/html; charset=UTF-8"
         pageEncoding="UTF-8" %>
<%@ page import="java.io.FileInputStream, java.io.InputStream, org.apache.poi.ss.usermodel.*" %>
<!DOCTYPE html>
<html lang="en">
<head>
    <link rel="stylesheet"
          href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css"/>
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>

    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Details-SR001</title>
    <link rel="stylesheet" href="CSS/header.css">
    <link rel="stylesheet" href="CSS/details.css">
    <link rel="stylesheet" href="CSS/degtest.css">

    <meta charset="UTF-8"/>
    <title>Data Browser</title>
    <style>

        .pagination { margin:10px 0; }
        .pagination a, .pagination span { margin-right:8px; text-decoration:none; }
        .pagination span.disabled { color:#999; }
        .pagination input { width:40px; text-align:center; }
    </style>
</head>

<header>
    <nav>
        <ul>
            <db_logo>
                <a href="#">sDSSA</a>
            </db_logo>
            <li><a href="index.jsp">Home</a>
            </li>
            <li><a href="browse.jsp">Browse</a>
            </li>
            <li><a href="search.jsp">Search</a>
                <div class="top_list">
                    <lib><a href="#">Gene</a></lib>
                    <lib><a href="#">Cell</a></lib>
                </div>
            </li>
            <li><a href="#">Help</a>
                <div class="top_list">
                    <lib><a href="#">Method</a></lib>
                    <lib><a href="#">Tutorial</a></lib>
                </div>
            </li>
            <li><a href="#">Download</a>
            </li>
            <un_logo>
                <img src="https://www.zju.edu.cn/_upload/tpl/0b/bf/3007/template3007/static/js/../../static/media/mlogo.80e02913954185729616ab6a1c6ae12d.svg" alt="" width="196" height="54">
            </un_logo>

        </ul>
    </nav>
</header>

<body>
<div class="main">
    <div class="title-div">
        <h2 class="page-title">Dataset Preview</h2>
    </div>

<%
    // 读取 Excel
    String excelPath = application.getRealPath("/WEB-INF/BrowseShow.xlsx");
    InputStream input = null;
    Workbook workbook = null;
    try {
        input = new FileInputStream(excelPath);
        workbook = WorkbookFactory.create(input);
        Sheet sheet = workbook.getSheetAt(0);

        // 分页参数
        int rowsPerPage = 10;
        int totalRows   = sheet.getLastRowNum();
        int totalPages  = (int) Math.ceil((double) totalRows / rowsPerPage);

        // 当前页
        String pageParam = request.getParameter("page");
        int pageNum = 1;
        try { pageNum = Integer.parseInt(pageParam); } catch(Exception ignore){}
        if (pageNum < 1)          pageNum = 1;
        if (pageNum > totalPages) pageNum = totalPages;

        int startRow = (pageNum - 1) * rowsPerPage + 1;
        int endRow   = Math.min(startRow + rowsPerPage - 1, totalRows);
%>
<table class="table-style-1">
    <thead>
    <tr>
        <th>SAID</th>
        <th>GSE</th>
        <th>GSM</th>
        <th>species</th>
        <th>disease information</th>
        <th>tissue information</th>
        <th>Details</th> <!-- 新增这一列 -->
    </tr>
    </thead>
    <tbody>
    <%
        for (int r = startRow; r <= endRow; r++) {
            Row row = sheet.getRow(r);
            if (row == null) continue;
            String said    = row.getCell(0, Row.MissingCellPolicy.CREATE_NULL_AS_BLANK).toString();
            String gse     = row.getCell(1, Row.MissingCellPolicy.CREATE_NULL_AS_BLANK).toString();
            String gsm     = row.getCell(2, Row.MissingCellPolicy.CREATE_NULL_AS_BLANK).toString();
            String species = row.getCell(3, Row.MissingCellPolicy.CREATE_NULL_AS_BLANK).toString();
            String disease = row.getCell(4, Row.MissingCellPolicy.CREATE_NULL_AS_BLANK).toString();
            String tissue  = row.getCell(5, Row.MissingCellPolicy.CREATE_NULL_AS_BLANK).toString();
    %>
    <tr>
        <td><%= said %></td>
        <td><%= gse %></td>
        <td><%= gsm %></td>
        <td><%= species %></td>
        <td><%= disease %></td>
        <td><%= tissue %></td>
        <td>
            <!-- 详情链接，点击跳转到 details.jsp -->
            <a href="details.jsp?said=<%= java.net.URLEncoder.encode(said, "UTF-8") %>">
                Details
            </a>
        </td>
    </tr>
    <%
        }
    %>
    </tbody>
</table>




<!-- 分页控件：Previous / Page 输入 / Next -->
<div class="pagination">
    <% if (pageNum > 1) { %>
    <a href="?page=<%= pageNum - 1 %>">Previous</a>
    <% } else { %>
    <span class="disabled">Previous</span>
    <% } %>

    <span>
        Page
        <form method="get" style="display:inline;">
          <input type="number" name="page"
                 min="1" max="<%= totalPages %>"
                 value="<%= pageNum %>"/>
          <button type="submit">Go</button>
        </form>
      </span>

    <% if (pageNum < totalPages) { %>
    <a href="?page=<%= pageNum + 1 %>">Next</a>
    <% } else { %>
    <span class="disabled">Next</span>
    <% } %>
</div>
<%
} catch (Exception e) {
%>
<p style="color:red;">Error loading data: <%= e.getMessage() %></p>
<%
    } finally {
        if (workbook != null) try { workbook.close(); } catch(Exception ignore){}
        if (input    != null) try { input.close();    } catch(Exception ignore){}
    }
%>
</div>
</body>
</html>
