<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>scSAID - Single-Cell Skin & Appendages Integrated Database</title>
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
            font-family: 'Arial', sans-serif;
        }

        body {
            background-color: #f8f9fa;
            color: #333;
            line-height: 1.6;
        }

        /* 导航栏样式 */
        header {
            background: linear-gradient(135deg, #2c3e50 0%, #4a6491 100%);
            color: white;
            padding: 0;
            position: sticky;
            top: 0;
            z-index: 100;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
        }

        nav ul {
            display: flex;
            justify-content: space-between;
            align-items: center;
            list-style: none;
            padding: 0 5%;
        }

        nav ul li {
            position: relative;
            padding: 15px 0;
        }

        nav a {
            color: white;
            text-decoration: none;
            padding: 10px 15px;
            border-radius: 4px;
            transition: background-color 0.3s;
        }

        nav a:hover {
            background-color: rgba(255, 255, 255, 0.1);
        }

        db_logo a {
            font-size: 24px;
            font-weight: bold;
            color: #ffd700 !important;
        }

        un_logo img {
            filter: brightness(0) invert(1);
        }

        /* 欢迎区域样式 */
        .index_welcome {
            position: relative;
            text-align: center;
            margin-bottom: 40px;
        }

        .index_welcome img {
            width: 100%;
            height: 490px;
            object-fit: cover;
            display: block;
        }

        .welcome {
            position: absolute;
            top: 50%;
            left: 50%;
            transform: translate(-50%, -50%);
            font-size: 48px;
            font-weight: bold;
            color: white;
            text-shadow: 2px 2px 4px rgba(0, 0, 0, 0.7);
        }

        .description {
            position: absolute;
            top: 65%;
            left: 50%;
            transform: translate(-50%, -50%);
            font-size: 20px;
            color: white;
            text-shadow: 1px 1px 3px rgba(0, 0, 0, 0.7);
        }

        /* 主要内容区域样式 */
        .basic {
            background-color: white;
            padding: 40px 5%;
            margin: 30px 0;
            border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.05);
        }

        .header {
            text-align: center;
            font-size: 32px;
            margin-bottom: 30px;
            color: #2c3e50;
            position: relative;
            padding-bottom: 15px;
        }

        .header:after {
            content: '';
            position: absolute;
            bottom: 0;
            left: 50%;
            transform: translateX(-50%);
            width: 80px;
            height: 4px;
            background: linear-gradient(to right, #3498db, #2c3e50);
            border-radius: 2px;
        }

        .content-container {
            display: flex;
            flex-wrap: wrap;
            gap: 30px;
        }

        .image-column {
            flex: 1;
            min-width: 300px;
        }

        .text-column {
            flex: 1;
            min-width: 300px;
        }

        .image-card {
            margin-bottom: 20px;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
            transition: transform 0.3s ease;
        }

        .image-card:hover {
            transform: translateY(-5px);
        }

        .image-card img {
            width: 100%;
            height: 296px;
            object-fit: cover;
            display: block;
        }

        .image-caption {
            padding: 15px;
            background-color: white;
        }

        .image-caption h3 {
            color: #2c3e50;
            margin-bottom: 8px;
        }

        .image-caption p {
            color: #7f8c8d;
            font-size: 14px;
        }

        .text_1 {
            font-size: 16px;
            line-height: 1.8;
            color: #34495e;
        }

        .text_1 b {
            color: #2c3e50;
        }

        /* 响应式设计 */
        @media (max-width: 768px) {
            .content-container {
                flex-direction: column;
            }

            .welcome {
                font-size: 36px;
            }

            .description {
                font-size: 16px;
            }

            nav ul {
                flex-wrap: wrap;
                justify-content: center;
            }

            nav ul li {
                padding: 10px 5px;
            }
        }

        /* 页脚样式 */
        footer {
            background: linear-gradient(135deg, #2c3e50 0%, #4a6491 100%);
            color: white;
            text-align: center;
            padding: 30px 5%;
            margin-top: 50px;
        }

        .footer-content {
            max-width: 1200px;
            margin: 0 auto;
        }

        .counter {
            background-color: white;
            color: #2c3e50;
            padding: 15px;
            border-radius: 8px;
            margin-top: 20px;
            display: inline-block;
        }
        .circles-container {
            display: flex;
            justify-content: center;
            gap: 150px;
            flex-wrap: wrap;
            margin: 60px 0;
            --blue-1: #e6f2ff;
            /* 最浅 */
            --blue-2: #cce5ff;
            --blue-3: #66b3ff;
            --blue-4: #0066cc;
            /* 主蓝 */
            --blue-5: #004d99;
            /* 最深 */

            .circles-container {
                display: flex;
                justify-content: center;
                gap: 40px;
                flex-wrap: wrap;
            }

            .circle {
                width: 130px;
                height: 130px;
                border-radius: 50%;
                background: rgba(255, 255, 255, .35);
                backdrop-filter: blur(10px);
                -webkit-backdrop-filter: blur(10px);
                border: 1.5px solid rgba(255, 255, 255, .6);
                display: flex;
                align-items: center;
                justify-content: center;
                color: var(--blue-5);
                font-size: 18px;
                font-weight: 600;
                cursor: pointer;
                transition: all .3s ease;
                box-shadow: 0 8px 24px rgba(0, 102, 204, .12);
                position: relative;
                overflow: hidden;
            }

            /* 蓝→白渐变外圈 */
            .circle::before {
                content: "";
                position: absolute;
                inset: -3px;
                /* 比原来大 3px */
                border-radius: 50%;
                /*background: linear-gradient(135deg, var(--blue-4) 0%, #ffffff 100%);*/
                z-index: -1;
            }

            .circle:hover {
                transform: translateY(-6px) scale(1.05);
                box-shadow: 0 12px 32px rgba(0, 102, 204, .2);
                background: rgba(255, 255, 255, .5);
            }
        }
        .block-container {
            --blue-1: #e6f2ff;
            /* 背景浅蓝 */
            --blue-4: #0066cc;
            /* 主蓝 */
            --blue-5: #004d99;
            /* 文字深蓝 */
            /*background: linear-gradient(135deg, var(--blue-1) 0%, #fff 100%);*/
            padding: 30px 0;
            display: flex;
            justify-content: center;
            gap: 30px;
            flex-wrap: wrap;

            .block {
                width: 160px;
                height: 90px;
                background: #fff;
                border: 1.5px solid var(--blue-4);
                border-radius: 12px;
                display: flex;
                align-items: center;
                justify-content: center;
                color: var(--blue-4);
                font-size: 18px;
                font-weight: 600;
                cursor: pointer;
                transition: all .3s ease;
                box-shadow: 0 4px 16px rgba(0, 102, 204, .1);
            }

            .block:hover {
                transform: translateY(-4px);
                box-shadow: 0 8px 24px rgba(0, 102, 204, .2);
                background: var(--blue-1);
            }
        }
    </style>
</head>
<body>

<header>
    <nav>
        <ul>
            <db_logo>
                <a href="#">scSAID</a>
            </db_logo>
            <li><a href="index.jsp">Home</a></li>
            <li><a href="browse.jsp">Browse</a></li>
            <li><a href="search.jsp">Search</a></li>
            <li><a href="#">Help</a></li>
            <li><a href="#">Download</a></li>
            <un_logo>
                <img src="https://www.zju.edu.cn/_upload/tpl/0b/bf/3007/template3007/static/js/../../static/media/mlogo.80e02913954185729616ab6a1c6ae12d.svg" alt="University Logo" width="196" height="54">
            </un_logo>
        </ul>
    </nav>
</header>

<div class="index_welcome">
    <img src="images/campus.png" alt="Banner">
    <div class="welcome">
        scSAID
    </div>
    <div class="description">
        Single-Cell <B>S</B>kin & <B>A</B>ppendages <B>I</B>ntegrated <B>D</B>atabase
    </div>
</div>
<div class="circles-container">
    <div class="circle">Browse</div>
    <div class="circle">Search</div>
    <div class="circle">Help</div>
    <div class="circle">Download</div>
</div>

<%--<div class="block-container">--%>
<%--    <div class="block">Browse</div>--%>
<%--    <div class="block">Search</div>--%>
<%--    <div class="block">Help</div>--%>
<%--    <div class="block">Download</div>--%>
<%--</div>--%>
<div class="basic">
    <div class="header">Data Overview</div>
    <div class="content-container">
        <!-- 左侧图片列 -->
        <div class="image-column">
            <div class="image-card">
                <img src="images/proportion.png" alt="Description">
                <div class="image-caption">
                    <h3>Sample Distribution</h3>
                    <p>Overview of samples by species and tissue type</p>
                </div>
            </div>

            <div class="image-card">
                <img src="images/date.png" alt="Date">
                <div class="image-caption">
                    <h3>Data Date</h3>
                    <p>Collected data are mainly from 2020 to 2025</p>
                </div>
            </div>

            <div class="image-card">
                <img src="images/tissue.png" alt="Tissue Sources">
                <div class="image-caption">
                    <h3>Tissue Sources</h3>
                    <p>Multiple tissues of skin and appendages are included</p>
                </div>
            </div>

            <div class="image-card">
                <img src="images/disease.png" alt="Disease type">
                <div class="image-caption">
                    <h3>Disease type</h3>
                    <p>Multiple types of disease are included</p>
                </div>
            </div>
        </div>

        <!-- 右侧文本列 -->
        <div class="text-column">
            <div class="text_1">
                Welcome to our comprehensive scRNA-seq database dedicated to skin and its appendages. This database contains data from over <b>1,000,000</b> cells derived from more than <b>600</b> samples across <b>100</b> independent experiments, including both human and mouse datasets. As one of the most extensive collections available to date, our database provides a platform for exploring the complex cellular landscapes and molecular mechanisms underlying skin biology and its associated structures.
                <br><br>
                The database includes detailed annotations for each sample, including species, gender, age, anatomical region, and experimental conditions. All data has been processed through a standardized pipeline to ensure consistency and comparability across studies.
                <br><br>
            </div>
        </div>
    </div>
</div>

<footer>
    <div class="footer-content">
        <p>© 2023 scSAID - Single-Cell Skin & Appendages Integrated Database</p>
        <p>Contact: info@scsaid.org | Privacy Policy | Terms of Use</p>
        <div class="counter">
            Total Visits: 12,487 | Today: 143
        </div>
    </div>
</footer>

</body>
</html>