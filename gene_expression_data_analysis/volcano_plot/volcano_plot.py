import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# LOAD DATA
df = pd.read_excel("DE_results_all_genes.xlsx")

# Extract log2FC and adjusted P-value columns
df["log2FC"] = df.iloc[:, 1]
df["adj_pval"] = df.iloc[:, 5]

# CALCULATE -LOG10(P-VALUE) 
df["neg_log10_pval"] = -np.log10(df["adj_pval"])

#  APPLY THRESHOLDS 
sig_fc = abs(df["log2FC"]) > 1
sig_p = df["adj_pval"] < 0.05

# Coloring rules
df["color"] = "black"  # default
df.loc[sig_p & ~sig_fc, "color"] = "blue"   # significant p-value only
df.loc[sig_p & sig_fc, "color"] = "red"     # both significant

# PLOT 
plt.figure(figsize=(10,7))
plt.scatter(df["log2FC"], df["neg_log10_pval"], c=df["color"], alpha=0.7, s=10)


# Threshold lines
plt.axhline(-np.log10(0.05), color="gray", linestyle="--", linewidth=1)
plt.axvline(-1, color="gray", linestyle="--", linewidth=1)
plt.axvline(1, color="gray", linestyle="--", linewidth=1)


# Label the top 20 most significant red genes
top_red = df[(df["color"] == "red")].nsmallest(20, "adj_pval")
for _, row in top_red.iterrows():
    plt.text(row["log2FC"], row["neg_log10_pval"], str(row.iloc[0]), 
fontsize=8, ha='right', va='bottom')


# Labels
plt.xlabel("log2(Fold Change)", fontsize=18)
plt.ylabel("-log10(Adjusted P-value)", fontsize=18)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig("volcano_plot_all_genes_TCGA_GBM.png", dpi=600, bbox_inches="tight")
