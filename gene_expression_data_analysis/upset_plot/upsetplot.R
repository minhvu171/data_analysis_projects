# Load required libraries
library(readxl)
library(UpSetR)

# Set file paths. Each of these files only has only column containing the genes.
# We want to see the overlapping content among the 5 files, so we create an upset plot.
apoptosis_file <- "apoptosis_genes.xlsx"
autophagy_file <- "autophagy_genes.xlsx"
ferroptosis_file <- "ferroptosis_genes.xlsx"
necroptosis_file <- "necroptosis_genes.xlsx"
pyroptosis_file <- "pyroptosis_genes.xlsx"

# Read gene lists 
apoptosis_genes <- read_excel(apoptosis_file)[[1]]
autophagy_genes <- read_excel(autophagy_file)[[1]]
ferroptosis_genes <- read_excel(ferroptosis_file)[[1]]
necroptosis_genes <- read_excel(necroptosis_file)[[1]]
pyroptosis_genes <- read_excel(pyroptosis_file)[[1]]

# Create a named list of gene sets
gene_lists <- list(
  Apoptosis = unique(apoptosis_genes),
  Autophagy = unique(autophagy_genes),
  Ferroptosis = unique(ferroptosis_genes),
  Necroptosis = unique(necroptosis_genes),
  Pyroptosis = unique(pyroptosis_genes)
)

# Convert the list into a binary membership matrix
# using fromList() from the UpSetR package
input_matrix <- fromList(gene_lists)

# Plot the UpSet plot
upset(input_matrix, 
      order.by = "freq", 
      sets = c("Apoptosis", "Autophagy", "Ferroptosis", "Necroptosis", "Pyroptosis"),
      sets.bar.color = "#1f77b4", 
      main.bar.color = "#ff7f0e", 
      text.scale = c(4, 4, 1, 4, 4, 4),
      matrix.color = "#2ca02c")

