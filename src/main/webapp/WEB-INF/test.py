import scanpy as sc
import pandas as pd

# 替换为你的 all_data.h5ad 文件的实际路径
# 这个路径应该和你 dash_app.py 中 H5AD_PATH 设置的路径一致
h5ad_file_path = "SkinDB_combined_human_integrated.h5ad" # 示例路径

try:
    # 加载 .h5ad 文件
    adata = sc.read_h5ad(h5ad_file_path)

    print("AnnData object loaded successfully.")
    print(f"Number of cells (observations): {adata.n_obs}")
    print(f"Number of genes (variables): {adata.n_vars}")

    print("\nColumns in adata.obs (observation metadata):")
    # 打印 adata.obs 中的所有列名
    print(adata.obs.columns.tolist())

    # 打印 adata.obs 的前几行，查看数据示例
    print("\nFirst 5 rows of adata.obs:")
    print(adata.obs.head())

    # 检查是否存在类似 'sample_id' 的列
    if 'sample_id' in adata.obs.columns:
        print("\n'sample_id' column found!")
        print(adata.obs['sample_id'].value_counts()) # 打印每个样本的细胞数量
    else:
        print("\n'sample_id' column NOT found. Looking for alternative sample identifiers...")
        # 尝试打印所有包含 'sample', 'batch', 'dataset' 等关键词的列
        for col in adata.obs.columns:
            if 'sample' in col.lower() or 'batch' in col.lower() or 'dataset' in col.lower() or 'orig' in col.lower():
                print(f"Potential sample identifier column: '{col}'")
                print(adata.obs[col].value_counts().head()) # 打印前几个值

except FileNotFoundError:
    print(f"Error: The file '{h5ad_file_path}' was not found. Please check the path.")
except Exception as e:
    print(f"An error occurred while loading or inspecting the .h5ad file: {e}")