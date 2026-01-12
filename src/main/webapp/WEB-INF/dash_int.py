import os
import pandas as pd
import scanpy as sc
from flask import Flask
from dash import Dash, html, dcc
from dash.dependencies import Input, Output
from urllib.parse import parse_qs
import plotly.express as px
import plotly.colors as pcolors

# Flask + Dash 初始化
server = Flask(__name__)
app = Dash(__name__, server=server, url_base_pathname="/dash/")

# Path to the master h5ad file (this should match the path in IntegrateServlet)
# 你需要根据你的实际部署环境来设置这个路径
H5AD_PATH = "SkinDB_combined_human_integrated.h5ad"

def load_integrated_adata(saids: list):
    """
    Loads the master .h5ad file, selects specified SAIDs, and performs integration and UMAP.
    """
    if not os.path.exists(H5AD_PATH):
        raise FileNotFoundError(f"Master H5AD file not found: {H5AD_PATH}")

    # Load the entire master AnnData object
    adata = sc.read_h5ad(H5AD_PATH)

    # 确保用于筛选的列存在于 adata.obs
    ACTUAL_SAMPLE_IDENTIFIER_COLUMN = 'GSM'
    if ACTUAL_SAMPLE_IDENTIFIER_COLUMN not in adata.obs.columns:
        raise ValueError(f"The '{ACTUAL_SAMPLE_IDENTIFIER_COLUMN}' column is missing in the master H5AD's .obs. "
                         f"Cannot filter by SAIDs. Available columns: {adata.obs.columns.tolist()}")

    # 过滤选定的 SAIDs (现在是 GSM 编号)
    adata_filtered = adata[adata.obs[ACTUAL_SAMPLE_IDENTIFIER_COLUMN].isin(saids)].copy()

    if adata_filtered.n_obs == 0:
        raise ValueError(f"No data found for the selected SAIDs ({saids}). "
                         f"Please check SAID values or master H5AD content and the '{ACTUAL_SAMPLE_IDENTIFIER_COLUMN}' column.")

    # 执行集成和 UMAP 步骤
    sc.pp.normalize_total(adata_filtered, target_sum=1e4)
    sc.pp.log1p(adata_filtered)
    sc.pp.highly_variable_genes(adata_filtered, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata_filtered = adata_filtered[:, adata_filtered.var['highly_variable']]
    sc.pp.scale(adata_filtered, max_value=10)
    sc.tl.pca(adata_filtered, svd_solver='arpack')
    sc.pp.neighbors(adata_filtered)
    sc.tl.umap(adata_filtered)
    sc.tl.leiden(adata_filtered, key_added="integrated_clusters")

    return adata_filtered

# Dash 应用布局
app.layout = html.Div(
    style={
        'width': '100%',
        'height': '100vh',
        'display': 'flex',
        'flexDirection': 'column'
    },
    children=[
        dcc.Location(id="url", refresh=False),
        html.Div(
            style={
                'padding': '10px',
                'textAlign': 'right', # 让按钮靠右
                'background': '#f8f8f8',
                'borderBottom': '1px solid #eee'
            },
            children=[
                # 返回整合页面按钮
                html.A(
                    html.Button('Back', style={'marginRight': '10px'}),
                    # 确保这里的 href 路径是正确的，替换为你的实际应用上下文路径
                    href="http://localhost:8080/scrna_website_test_war/search.jsp", # 根据你的实际上下文路径调整
                    target="_self"
                ),
                # 导出为PNG功能通过 dcc.Graph 的 config 属性实现，无需额外按钮
                # 用户可以通过图表右上角的相机图标下载
            ]
        ),
        html.Div(
            style={
                'flex': 1,
                'display': 'flex'
            },
            children=[
                dcc.Graph(
                    id='umap-plot',
                    style={'width': '100%', 'height': '100%'},
                    # 启用 modebar 使得用户可以通过图表右上角的相机图标下载PNG
                    config={'responsive': True, 'displayModeBar': True, 'toImageButtonOptions': {'format': 'png', 'filename': 'integrated_umap_plot'}}
                )
            ]
        )
    ]
)

@app.callback(
    Output("umap-plot", "figure"),
    [Input("url", "search")]
)
def update_umap(search):
    # 解析 URL 参数
    params = parse_qs(search.lstrip("?"))
    saids_str = params.get("saids", [None])[0]

    if not saids_str:
        return px.scatter(title="Error: No SAIDs (GSMs) provided in URL.")

    selected_saids = saids_str.split(',')

    if len(selected_saids) < 2:
        return px.scatter(title="Error: Please select at least two datasets for integration.")

    # 加载并集成 AnnData
    try:
        adata = load_integrated_adata(selected_saids)
    except Exception as e:
        print(f"Error loading or processing data: {e}")
        return px.scatter(title=f"Error loading or processing data: {e}. Please check server logs.")

    # ----------------------------------------------------------------------
    # 新的着色和图例逻辑
    # ----------------------------------------------------------------------
    REQUIRED_COLS_FOR_COLORING = ['GSM', 'scanvi_pred']
    for col in REQUIRED_COLS_FOR_COLORING:
        if col not in adata.obs.columns:
            raise ValueError(
                f"Required column '{col}' for coloring is missing in adata.obs. "
                f"Available columns: {adata.obs.columns.tolist()}"
            )

    # 构造 UMAP DataFrame
    df_umap = pd.DataFrame(
        adata.obsm["X_umap"],
        index=adata.obs.index,
        columns=["UMAP1", "UMAP2"]
    )

    # 将需要的元数据列添加到 df_umap 中
    df_umap['GSM'] = adata.obs['GSM']
    df_umap['scanvi_pred'] = adata.obs['scanvi_pred']

    # 如果存在 integrated_clusters，也将其添加到 df_umap 以便悬停显示
    if 'integrated_clusters' in adata.obs.columns:
        df_umap['integrated_clusters'] = adata.obs['integrated_clusters'].astype(str)

    # ----------------------------------------------------------------------
    # 生成一个包含足够多独特颜色的颜色序列
    # 拼接多个 Plotly 内置颜色序列以确保有足够多的颜色
    color_sequence = pcolors.qualitative.Alphabet + pcolors.qualitative.Light24 + pcolors.qualitative.Plotly + pcolors.qualitative.Bold
    # ----------------------------------------------------------------------

    # 准备悬停数据显示的列
    hover_data_cols = {
        "GSM": True,          # 鼠标悬停时显示原始 GSM
        "scanvi_pred": True   # 鼠标悬停时显示 scanvi_pred 细胞类型
    }
    # 如果有集成聚类结果，也添加到悬停信息
    if 'integrated_clusters' in df_umap.columns:
        hover_data_cols['integrated_clusters'] = True

    # 创建散点图
    fig = px.scatter(
        df_umap,
        x="UMAP1",
        y="UMAP2",
        color="scanvi_pred",  # 仅根据 'scanvi_pred' 着色
        title=f"Integrated UMAP Plot for GSMs: {', '.join(selected_saids)}",
        hover_data=hover_data_cols,  # 设置悬停信息
        labels={"color": "Cell Type"}, # 设置图例标题
        color_discrete_sequence=color_sequence # 使用自定义的颜色序列
    )

    fig.update_layout(
        autosize=True,
        margin=dict(l=10, r=10, t=40, b=10),
        legend=dict(
            orientation="v", # 垂直排列图例项
            yanchor="top",
            y=1,
            xanchor="left",
            x=1.02 # 将图例放置在图表右侧，稍微超出一点以避免覆盖数据
        )
    )

    return fig

if __name__ == "__main__":
    server.run(debug=True, host="0.0.0.0", port=8050)