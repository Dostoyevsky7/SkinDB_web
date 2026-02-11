<%@ page language="java"
         contentType="text/html; charset=UTF-8"
         pageEncoding="UTF-8" %>
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Contact - scSAID</title>

    <!-- Google Fonts -->
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Cormorant+Garamond:wght@400;500;600;700&family=Source+Sans+3:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">

    <!-- Design System -->
    <link rel="stylesheet" href="CSS/design-system.css">
    <link rel="stylesheet" href="CSS/header.css">
    <link rel="stylesheet" href="CSS/construction-modal-simple.css">

    <style>
        body {
            background-color: #faf8f5;
        }

        .contact-page {
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

        .contact-content {
            max-width: 1000px;
            margin: 0 auto;
            padding: 3rem 2rem 4rem;
        }

        .contact-card {
            background: #ffffff;
            border-radius: 16px;
            box-shadow: 0 4px 12px rgba(26, 35, 50, 0.08);
            padding: 2rem;
        }

        .contact-card p {
            margin: 0 0 1rem;
            line-height: 1.7;
            color: #4b5563;
        }

        .contact-email {
            font-family: 'JetBrains Mono', monospace;
            color: #1a2332;
        }

        @media (max-width: 768px) {
            .page-header {
                padding: 3.5rem 0;
            }

            .contact-content {
                padding: 2.5rem 1.25rem 3rem;
            }

            .contact-card {
                padding: 1.5rem;
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
            <a href="feedback" class="main-nav__link">Feedback</a>
            <a href="contact" class="main-nav__link main-nav__link--active">Contact</a>
        </nav>
        <a href="https://zje.zju.edu.cn/zje/main.htm" target="_blank">
            <img src="images/ZJE_Logo.png" alt="ZJE - Zhejiang University" class="university-logo">
        </a>
    </div>
</header>

<main class="contact-page">
    <div class="page-header">
        <div class="page-header__content">
            <span class="page-header__eyebrow">Get in Touch</span>
            <h1 class="page-header__title">Contact</h1>
            <p class="page-header__description">
                Reach out for dataset questions, collaborations, or publication guidance.
            </p>
        </div>
    </div>

    <div class="contact-content">
        <div class="contact-card">
            <p>Contact our active developer (contact one and cc another):</p>
            <p class="contact-email">ethanshen111 [at] gmail [dot] com</p>
            <p class="contact-email">yixiangren99 [at] gmail [dot] com</p>
            <p>For publication-related matters contact our correspondence author.</p>
        </div>
    </div>
</main>

<script src="JS/construction-modal-simple.js"></script>
</body>
</html>
