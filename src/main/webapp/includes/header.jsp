<%--
  Shared Header Component for scSAID
  Include this file in all JSP pages: <%@ include file="includes/header.jsp" %>
--%>
<!-- Header -->
<header class="site-header">
    <div class="container">
        <a href="index.jsp" class="site-logo">scSAID</a>
        <nav class="main-nav">
            <a href="index.jsp" class="main-nav__link">Home</a>

            <%-- Browse Dropdown --%>
            <div class="main-nav__item">
                <a href="browse.jsp" class="main-nav__link">Browse</a>
                <div class="main-nav__dropdown">
                    <a href="browse.jsp" class="main-nav__dropdown-link">All Data</a>
                    <a href="browse.jsp?filter=scrna" class="main-nav__dropdown-link">scRNA-seq</a>
                    <a href="browse.jsp?filter=spatial" class="main-nav__dropdown-link">Spatial</a>
                </div>
            </div>

            <%-- Tools Dropdown --%>
            <div class="main-nav__item">
                <a href="gene-search.jsp" class="main-nav__link">Tools</a>
                <div class="main-nav__dropdown">
                    <a href="gene-search.jsp" class="main-nav__dropdown-link">Gene Search</a>
                    <a href="visualization.jsp" class="main-nav__dropdown-link">Cell Clustering</a>
                </div>
            </div>

            <%-- Help Dropdown --%>
            <div class="main-nav__item">
                <a href="help.jsp?topic=faq" class="main-nav__link">Help</a>
                <div class="main-nav__dropdown">
                    <a href="help.jsp?topic=faq" class="main-nav__dropdown-link">FAQ</a>
                    <a href="help.jsp?topic=pipeline" class="main-nav__dropdown-link">Pipeline</a>
                    <a href="contact.jsp" class="main-nav__dropdown-link">Contact</a>
                    <a href="feedback.jsp" class="main-nav__dropdown-link">Feedback</a>
                </div>
            </div>

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
