import os
import pandas as pd
import scanpy as sc
from flask import Flask
from dash import Dash, html, dcc
from dash.dependencies import Input, Output
from urllib.parse import parse_qs
import plotly.express as px

# Flask + Dash 初始化
server = Flask(__name__)
app = Dash(__name__, server=server, url_base_pathname="/dash/")

# Path to the master h5ad file (this should match the path in IntegrateServlet)
# 你需要根据你的实际部署环境来设置这个路径
# 例如：如果你的 dash_app.py 和 all_data.h5ad 在同一个WEB-INF/data目录下，可以这样设置
# H5AD_PATH = os.path.join(os.path.dirname(__file__), 'data', 'SkinDB_combined_human_integrated.h5ad')
# 或者如果你手动放置在项目根目录的某个位置
H5AD_PATH = "SkinDB_combined_human_integrated.h5ad" # 保持你现有的相对路径，但请确保此路径在dash app运行环境中可访问

def load_integrated_adata(saids: list):
    """
    Loads the master .h5ad file, selects specified SAIDs, and performs integration and UMAP.
    """
    if not os.path.exists(H5AD_PATH):
        raise FileNotFoundError(f"Master H5AD file not found: {H5AD_PATH}")

    # Load the entire master AnnData object
    adata = sc.read_h5ad(H5AD_PATH)

    # ----------------------------------------------------------------------
    # 核心修改点：将用于筛选的列名从 'sample_id' 更改为 'GSM'
    # 根据你提供的 .h5ad 结构，'GSM' 列包含你想要筛选的 GSM 编号
    ACTUAL_SAMPLE_IDENTIFIER_COLUMN = 'GSM'
    # ----------------------------------------------------------------------

    # 确保用于筛选的列存在于 adata.obs
    if ACTUAL_SAMPLE_IDENTIFIER_COLUMN not in adata.obs.columns:
        raise ValueError(f"The '{ACTUAL_SAMPLE_IDENTIFIER_COLUMN}' column is missing in the master H5AD's .obs. "
                         f"Cannot filter by SAIDs. Available columns: {adata.obs.columns.tolist()}")

    # 过滤选定的 SAIDs (现在是 GSM 编号)
    # adata.obs[ACTUAL_SAMPLE_IDENTIFIER_COLUMN] 将使用 'GSM' 列的值进行筛选
    adata_filtered = adata[adata.obs[ACTUAL_SAMPLE_IDENTIFIER_COLUMN].isin(saids)].copy()

    if adata_filtered.n_obs == 0:
        raise ValueError(f"No data found for the selected SAIDs ({saids}). "
                         f"Please check SAID values or master H5AD content and the '{ACTUAL_SAMPLE_IDENTIFIER_COLUMN}' column.")

    # 执行集成和 UMAP 步骤
    # 这些步骤对于正确的集成至关重要。你可能需要根据你使用的具体集成方法进行调整
    # （例如 harmony, scVI 等）。这里提供一个基本的 Scanpy 工作流。

    # 1. 归一化和对数转换（如果尚未完成）
    sc.pp.normalize_total(adata_filtered, target_sum=1e4)
    sc.pp.log1p(adata_filtered)

    # 2. 识别高变基因（可选，但对于大型数据集是良好实践）
    sc.pp.highly_variable_genes(adata_filtered, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata_filtered = adata_filtered[:, adata_filtered.var['highly_variable']]

    # 3. 缩放数据
    sc.pp.scale(adata_filtered, max_value=10)

    # 4. 运行 PCA
    sc.tl.pca(adata_filtered, svd_solver='arpack')

    # 5. 计算邻居（UMAP 的重要步骤）
    # 如果要集成多个数据集，你可以在这里考虑使用集成方法
    # 为了简化，我们暂时使用标准的邻居和 UMAP，如果数据预处理得当，它将隐式集成
    sc.pp.neighbors(adata_filtered)

    # 6. 运行 UMAP
    sc.tl.umap(adata_filtered)

    # 7. 聚类数据（可选，但对于按细胞类型着色 UMAP 非常有用）
    sc.tl.leiden(adata_filtered, key_added="integrated_clusters") # 或 'louvain'

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
                'flex': 1,
                'display': 'flex'
            },
            children=[
                dcc.Graph(
                    id='umap-plot',
                    style={'width': '100%', 'height': '100%'},
                    config={'responsive': True}
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

    # 从 URL 获取 SAID (GSM 编号) 列表。IntegrateServlet 将以 'saids=GSM1,GSM2,GSM3' 格式传递
    saids_str = params.get("saids", [None])[0]

    if not saids_str:
        return px.scatter(title="Error: No SAIDs (GSMs) provided in URL.")

    # 将逗号分隔的 SAID (GSM 编号) 字符串转换为列表
    selected_saids = saids_str.split(',')

    if len(selected_saids) < 2:
        return px.scatter(title="Error: Please select at least two datasets for integration.")

    # 加载并集成 AnnData
    try:
        adata = load_integrated_adata(selected_saids)
    except Exception as e:
        # 记录完整的错误以便调试
        print(f"Error loading or processing data: {e}")
        return px.scatter(title=f"Error loading or processing data: {e}. Please check server logs.")

    # 构造 UMAP DataFrame
    df_umap = pd.DataFrame(
        adata.obsm["X_umap"],
        index=adata.obs.index,
        columns=["UMAP1", "UMAP2"]
    )

    # ----------------------------------------------------------------------
    # 优化着色和悬停信息：
    # 优先使用 'integrated_clusters' 进行着色（显示聚类结果）
    # 如果没有 'integrated_clusters'，则回退到 'GSM'（显示原始数据集来源）
    color_by_key = "integrated_clusters"
    if color_by_key not in adata.obs.columns:
        color_by_key = "GSM" # <--- 回退到 'GSM'

    df_umap["color_group"] = adata.obs[color_by_key].astype(str)
    df_umap['GSM'] = adata.obs['GSM']

    # 确保悬停数据包含 'GSM' 以显示原始样本 ID
    hover_data_cols = ["color_group", "GSM"] # <--- 确保包含 'GSM'
    # 如果你的 adata.obs 中还有其他你想在悬停时看到的重要信息，可以添加到这里
    # 例如：hover_data_cols.append("cell_type")

    # ----------------------------------------------------------------------

    color_by_key = "integrated_clusters"
    legend_title = "Integrated Clusters" # 默认图例标题
    if color_by_key not in adata.obs.columns:
        color_by_key = "GSM"
        legend_title = "Original GSM" # 如果回退到GSM，则图例标题变为Original GSM

    df_umap["color_group"] = adata.obs[color_by_key].astype(str)


    # 创建散点图
    fig = px.scatter(
        df_umap,
        x="UMAP1",
        y="UMAP2",
        color="color_group", # 使用 color_group 进行着色
        title=f"Integrated UMAP Plot for GSMs: {', '.join(selected_saids)}",
        hover_data={"color_group": True, "GSM": True}, # 明确显示 color_group 和 GSM，True表示显示值
        # ----------------------------------------------------------------------
        # 添加 labels 参数来设置图例标题
        labels={"color_group": legend_title} # 将 color_group 的图例标题设置为 legend_title
        # ----------------------------------------------------------------------
    )

    fig.update_layout(
        autosize=True,
        margin=dict(l=10, r=10, t=40, b=10),
        # ----------------------------------------------------------------------
        # 调整图例位置，使其在右侧
        legend=dict(
            orientation="v", # 垂直排列
            yanchor="top",
            y=1,
            xanchor="left",
            x=1.02 # 将图例放置在图表右侧（稍微超出一点以避免覆盖数据）
        )
        # ----------------------------------------------------------------------
    )

    return fig

if __name__ == "__main__":
    server.run(debug=True, host="0.0.0.0", port=8050)