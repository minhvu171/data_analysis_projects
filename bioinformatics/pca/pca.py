import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# LOAD DATA 
df = pd.read_csv("2025-07-23_TCGA_GBM_expression.txt", sep="\t")
# Each sample (each row) is a data point, each gene (each column) is a dimension

# First column is sample IDs
sample_ids = df.iloc[:, 0]
expression_matrix = df.iloc[:, 1:]

# RUN PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(expression_matrix)

# Variance explained
explained_var = pca.explained_variance_ratio_ * 100

# DEFINE SAMPLE GROUPS
normal_samples = {
    "TCGA.06.0673", "TCGA.06.0675", "TCGA.06.0676",
    "TCGA.06.0678", "TCGA.06.0680", "TCGA.06.0681",
    "TCGA.08.0623", "TCGA.08.0625", "TCGA.08.0626",
    "TCGA.08.0627"
}

is_normal = sample_ids.isin(normal_samples)

# PLOT 
plt.figure(figsize=(8, 6))

# Normal samples
plt.scatter(
    pca_result[is_normal, 0],
    pca_result[is_normal, 1],
    c="teal", label="Normal", alpha=0.7
)

# GBM samples
plt.scatter(
    pca_result[~is_normal, 0],
    pca_result[~is_normal, 1],
    c="orange", label="GBM", alpha=0.7
)

plt.xlabel(f"PC1 ({explained_var[0]:.2f}%)", fontsize=18)
plt.ylabel(f"PC2 ({explained_var[1]:.2f}%)", fontsize=18)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend()
plt.tight_layout()
plt.savefig("TCGA_GBM dataset - pca plot normal vs GBM.png", dpi=600)
plt.show()
