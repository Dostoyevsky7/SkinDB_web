<%@ page contentType="text/html;charset=UTF-8" %>
<%@ taglib uri="http://java.sun.com/jsp/jstl/core" prefix="c" %>

<html>
<head>
    <title>Dataset Integrate</title>
    <style>
        body { font-family: sans-serif; margin: 30px; }
        iframe { border: 1px solid #ccc; margin-top: 10px; }
        .dataset { margin-bottom: 50px; padding-bottom: 20px; border-bottom: 1px dashed #ccc; }
    </style>
</head>
<body>

<h1>All Datasets</h1>

<c:forEach var="item" items="${datasetList}">
    <div class="dataset">
        <!-- 数据编号，点击后跳转到 Dash 详情页 -->
        <h2>
            <a href="http://localhost:5000/details/${item.said}" target="_blank">
                    ${item.said}
            </a>
        </h2>

        <!-- GSE / GSM 信息 -->
        <p>
            <strong>GSE:</strong> ${item.gse}
            &nbsp;|&nbsp;
            <strong>GSM:</strong> ${item.gsm}
        </p>

        <!-- 标题与摘要 -->
        <p><strong>Title:</strong> ${item.title}</p>
        <p><strong>Summary:</strong> ${item.summary}</p>

        <!-- 交互式 UMAP 图 -->
        <h3>Cell Clustering</h3>
        <iframe
                src="http://localhost:5000/dash/?sample_id=${item.said}"
                width="100%"
                height="600px"
                frameborder="0">
        </iframe>
    </div>
</c:forEach>

</body>
</html>
