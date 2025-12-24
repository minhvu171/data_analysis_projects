"""
Generates a heatmap of differential expression for orphan GPCR genes in GBM (glioblastoma) tumor samples
relative to normal samples using TCGA_GBM expression data.

Data sources:
- TCGA_GBM expression data (GlioVis)
- Orphan GPCR gene list (Class A)
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# === Data files ===
# Expression data is obtained by filtering the genes of interest from the TCGA_GBM dataset from GlioVis,
# https://gliovis.bioinfo.cnio.es/
expression_file = "adult_genes_and_expression_data.xlsx" # Each row corresponds to a sample, each column corresponds to a gene
orphan_gpcr_file = "Orphan_GPCR_Class_A.xlsx" # Only consist one column of genes

# Load gene expression data to a dataframe (sample IDs as rows, genes as columns)
expr_df = pd.read_excel(expression_file, index_col=0)

# Load orphan GPCR gene list
orphan_gpcr_genes = pd.read_excel(orphan_gpcr_file, usecols=[0]).iloc[:, 0].dropna().unique()

# Filter expression data to include only orphan GPCR genes
orphanGPCR_expr_df = expr_df.loc[:, expr_df.columns.intersection(orphan_gpcr_genes)]
print(f"Number of orphan GPCR genes in the TCGA_GBM datase is {orphanGPCR_expr_df.shape[1]}")

# === Organize samples ===
normal_samples = [
    "TCGA.06.0673", "TCGA.06.0675", "TCGA.06.0676",
    "TCGA.06.0678", "TCGA.06.0680", "TCGA.06.0681",
    "TCGA.08.0623", "TCGA.08.0625", "TCGA.08.0626", "TCGA.08.0627"
]

# Ensure sample IDs are strings and index is consistent
orphanGPCR_expr_df.index = orphanGPCR_expr_df.index.astype(str)

# Split samples into normal and tumor groups
normal_df = orphanGPCR_expr_df.loc[orphanGPCR_expr_df.index.isin(normal_samples)]
tumor_df = orphanGPCR_expr_df.loc[~orphanGPCR_expr_df.index.isin(normal_samples)]

# === For each gene, compute differential expression between each tumor sample and the mean of normal group ===
# Mean expression across normal samples for each gene
normal_mean = normal_df.mean()

# Subtract the mean of normal samples from each tumor sample -> This is the differential expression value for each gene
diff_expr = tumor_df.subtract(normal_mean)

# === Plotting heatmap ===
sns.set_theme(font_scale=0.8)
g = sns.clustermap(
    diff_expr.T,  # Transpose: genes as rows
    cmap="viridis", 
    figsize=(12, 10), 
    row_cluster=True, 
    col_cluster=True,
     xticklabels=False
)

g.ax_col_dendrogram.set_visible(False)
g.ax_heatmap.set_xlabel("GBM Samples")
colorbar = g.ax_heatmap.collections[0].colorbar
colorbar.set_ticks([-4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7])
g.savefig("GBM_samples_orphan_gpcr_heatmap.png", dpi=600, bbox_inches='tight')

plt.show()
