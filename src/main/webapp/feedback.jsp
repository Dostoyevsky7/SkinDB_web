<%@ page language="java"
         contentType="text/html; charset=UTF-8"
         pageEncoding="UTF-8" %>
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Feedback - scSAID</title>

    <!-- Google Fonts -->
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Cormorant+Garamond:wght@400;500;600;700&family=Urbanist:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">

    <!-- Design System -->
    <link rel="stylesheet" href="CSS/design-system.css">
    <link rel="stylesheet" href="CSS/header.css">
    <link rel="stylesheet" href="CSS/construction-modal-simple.css">

    <style>
        body {
            background-color: #faf8f5;
        }

        .feedback-page {
            min-height: 100vh;
            padding-top: 72px;
        }

        .page-header {
            background: #1a2332;
            padding: 4.5rem 0;
            text-align: center;
        }

        .page-header__content {
            max-width: 800px;
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
            margin-bottom: 1.5rem;
        }

        .page-header__title {
            font-family: 'Cormorant Garamond', Georgia, serif;
            font-size: clamp(2.5rem, 5vw, 3.5rem);
            font-weight: 500;
            color: #ffffff;
            margin: 0 0 1rem;
        }

        .page-header__description {
            font-size: 1.15rem;
            color: rgba(255, 255, 255, 0.75);
            margin: 0;
            line-height: 1.7;
        }

        .feedback-content {
            max-width: 1200px;
            margin: 0 auto;
            padding: 3rem 2rem 4rem;
        }

        .feedback-card {
            background: #ffffff;
            border-radius: 16px;
            box-shadow: 0 4px 12px rgba(26, 35, 50, 0.08);
            padding: 2rem;
        }

        .feedback-iframe {
            width: 100%;
            height: 1100px;
            border: 0;
            border-radius: 12px;
        }

        .feedback-fallback {
            margin-top: 1.5rem;
            text-align: center;
        }

        @media (max-width: 768px) {
            .page-header {
                padding: 3.5rem 0;
            }

            .feedback-content {
                padding: 2.5rem 1.25rem 3rem;
            }

            .feedback-card {
                padding: 1.25rem;
            }

            .feedback-iframe {
                height: 1350px;
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
            <a href="gene-search.jsp" class="main-nav__link">Search</a>
            <a href="download.jsp" class="main-nav__link">Download</a>
            <div class="main-nav__item">
                <a href="help?topic=faq" class="main-nav__link">Help</a>
                <div class="main-nav__dropdown">
                    <a href="help?topic=faq" class="main-nav__dropdown-link">FAQ</a>
                    <a href="help?topic=methods" class="main-nav__dropdown-link">Methods</a>
                    <a href="help?topic=markers" class="main-nav__dropdown-link">Markers</a>
                    <a href="help?topic=pipeline" class="main-nav__dropdown-link">Pipeline</a>
                    <a href="help?topic=usage" class="main-nav__dropdown-link">Usage</a>
                </div>
            </div>
            <a href="feedback" class="main-nav__link main-nav__link--active">Feedback</a>
            <a href="contact" class="main-nav__link">Contact</a>
        </nav>
        <a href="https://zje.zju.edu.cn/zje/main.htm" target="_blank">
            <img src="images/ZJE_Logo.png" alt="ZJE - Zhejiang University" class="university-logo">
        </a>
    </div>
</header>

<main class="feedback-page">
    <div class="page-header">
        <div class="page-header__content">
            <span class="page-header__eyebrow">Community</span>
            <h1 class="page-header__title">Feedback</h1>
            <p class="page-header__description">
                Tell us what you think. Your feedback helps us improve.
            </p>
        </div>
    </div>

    <div class="feedback-content">
        <div class="feedback-card">
            <iframe
                    class="feedback-iframe"
                    title="Feedback form"
                    src="https://docs.google.com/forms/d/e/1FAIpQLSflIMVjxnvApZ7I0uUwPdvt9_C7self4p-a3K2NoC6T8YLgLg/viewform?embedded=true"
                    loading="lazy"
            ></iframe>

            <div class="feedback-fallback">
                <a
                        class="btn btn--outline"
                        href="https://docs.google.com/forms/d/e/1FAIpQLSflIMVjxnvApZ7I0uUwPdvt9_C7self4p-a3K2NoC6T8YLgLg/viewform?usp=dialog"
                        target="_blank"
                        rel="noopener noreferrer"
                >
                    Open the feedback form in a new tab
                </a>
            </div>
        </div>
    </div>
</main>

<script src="JS/construction-modal-simple.js"></script>
</body>
</html>
