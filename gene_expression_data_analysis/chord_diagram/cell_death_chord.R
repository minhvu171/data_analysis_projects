# === R SCRIPT FOR CHORD DIAGRAM ===
# There are 3 main parts: READ DATA, CREAT CHORD DATA, PLOT CHORD DIAGRAM.

# Required libraries
library(readxl)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
library(grDevices)

# === CONFIG ===
sheet_map <- list(
  "4" = "Apoptosis",
  "5" = "Autophagy",
  "6" = "Necroptosis",
  "7" = "Ferroptosis",
  "8" = "Pyroptosis"
)

# Fixed colors for cell death types
cell_death_colors <- c(
  "Apoptosis" = "#e41a1c",
  "Autophagy" = "#377eb8",
  "Necroptosis" = "#4daf4a",
  "Ferroptosis" = "#984ea3",
  "Pyroptosis" = "#ff7f00"
)

# === READ DATA ===
file <- "Common_genes_from_three_different_treatment_conditions.xlsx"
# The task here is to read from each sheet (each represents a cell death program that
# the genes take part in) in this excel file a gene list. There will be
# duplicated genes between these gene lists, and that's what we want to visualize.
edges <- data.frame(CellDeath = character(), Gene = character(), stringsAsFactors = FALSE)
gene_info <- data.frame(Gene = character(), log2FC = numeric(), stringsAsFactors = FALSE)

for (sheet_idx in names(sheet_map)) {
  death_type <- sheet_map[[sheet_idx]]
  df <- read_excel(file, sheet = as.numeric(sheet_idx))
  colnames(df)[1:2] <- c("Gene", "log2FC")
  df <- df[!is.na(df$Gene), ]
  
  edges <- rbind(edges, data.frame(CellDeath = death_type, Gene = df$Gene))
  gene_info <- rbind(gene_info, df[, 1:2])
}

# Remove duplicates in gene_info
gene_info <- gene_info[!duplicated(gene_info$Gene), ]


# === CREATE CHORD DATA ===
links <- edges

# Enforce a stable sector order: cell deaths (fixed order) first, then genes (alphabetical)
cell_deaths <- unname(unlist(sheet_map))
all_genes <- sort(gene_info$Gene)
all_nodes <- c(cell_deaths, all_genes)

# Assign colors: cell deaths fixed, genes gradient by log2FC
min_fc <- min(gene_info$log2FC, na.rm = TRUE)
max_fc <- max(gene_info$log2FC, na.rm = TRUE)

# Gradient color function for log2FC (blue -> white -> red)
gene_colors <- colorRamp2(c(min_fc, 0, max_fc), c("blue", "white", "red"))

node_colors <- sapply(all_nodes, function(x) {
  if (x %in% names(cell_death_colors)) cell_death_colors[x] else gene_colors(gene_info$log2FC[match(x, gene_info$Gene)])
})
names(node_colors) <- all_nodes

# Link colors inherit from cell death type (first column of links)
link_cols <- cell_death_colors[links$CellDeath]


# === PLOT CHORD DIAGRAM ===
# Save as high-quality PNG 
png("cell_death_chord.png", width = 7500, height = 6500, res = 450)

circos.clear()
par(mar = c(1, 1, 1, 1))
# Slight gaps and a bit more space after the last sector
circos.par(gap.after = c(rep(2, length(all_nodes) - 1), 8), track.margin = c(0.005, 0.005))

chordDiagram(
  x = links,
  order = all_nodes, # stable order
  grid.col = node_colors, # named mapping ensures consistent colors
  annotationTrack = c("grid"), 
  preAllocateTracks = list(track.height = 0.06),
  col = link_cols, # chord color from cell-death section
  directional = 0
)

circos.track(track.index = 1, panel.fun = function(x, y) {
  sector_name <- CELL_META$sector.index
  # is_type <- sector_name %in% names(cell_death_colors)
  cex_val <- 1.3 
  font_val <- 2
  adj_val <- c(0, 0.5)
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], sector_name,
              facing = "clockwise", niceFacing = TRUE, adj = adj_val,
              cex = cex_val, font = font_val)
}, bg.border = NA)


# === LEGENDS ===
lgd_cell <- Legend(
  at = names(cell_death_colors), type = "grid",
  legend_gp = gpar(fill = unname(cell_death_colors)),
  title = "Cell Death Types",
  labels_gp = gpar(fontsize = 14, fontface = "bold"),  # legend labels
  title_gp  = gpar(fontsize = 16, fontface = "bold")   # legend title
)

lgd_fc <- Legend(
  col_fun = gene_colors, title = "log2FC",
  labels_gp = gpar(fontsize = 14, fontface = "bold"),
  title_gp  = gpar(fontsize = 16, fontface = "bold")
)

# Pack and draw legends near the right edge
leg <- packLegend(lgd_cell, lgd_fc, direction = "vertical")
draw(leg, x = unit(1, "npc") - unit(0.01, "cm"), y = unit(0.45, "npc"), just = "right")

dev.off()
