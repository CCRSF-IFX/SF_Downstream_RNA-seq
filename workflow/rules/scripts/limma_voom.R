# Load necessary libraries
# These libraries are essential for performing differential expression analysis (limma and edgeR) 
# and visualizing results (pheatmap).
library(limma)
library(edgeR)
library(pheatmap)

# Parse Command-Line Arguments:
# These arguments provide the file paths for raw counts, sample metadata, contrasts, and output directory.
args <- commandArgs(trailingOnly = TRUE)
raw_counts_path <- args[1]
coldata_path <- args[2]
contrasts_path <- args[3]
output_dir <- args[4]

# Read input data
# The raw counts data, sample metadata, and contrasts are loaded into R for processing.
counts <- read.table(raw_counts_path, header = TRUE, row.names = 1)
coldata <- read.table(coldata_path, header = TRUE)
contrasts <- read.table(contrasts_path, header = FALSE, stringsAsFactors = FALSE)

# Print a preview of the input files for debugging purposes
# This helps ensure the data is correctly loaded before proceeding with the analysis.
cat("Counts data:\n")
print(head(counts))
cat("Coldata:\n")
print(head(coldata))
cat("Contrasts:\n")
print(head(contrasts))

# Setup group factors based on experimental conditions
# 'group' represents experimental groups (e.g., treated vs control) and is converted to a factor.
# The factor levels are printed to ensure the groups are correct.
group <- factor(coldata$condition)  # Ensure the correct column name is used
cat("Factor levels in 'group':\n")
print(levels(group))

# Ensure the group factor has at least two levels for valid contrast comparisons
if (length(levels(group)) < 2) {
  stop("The group factor must have at least two levels.")
}

# Create DGEList object
# The raw count data is encapsulated in a DGEList object, and TMM normalization is applied to adjust for 
# differences in library sizes across samples.
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)

# Create the design matrix
# The design matrix specifies how the experimental groups are structured (e.g., control vs treated).
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
cat("Design matrix:\n")
print(design)

# Estimate dispersions
# This step estimates the variability of gene expression between replicates.
# Dispersion estimates are used in the model to detect true biological differences.
dge <- estimateDisp(dge, design)

# Transform raw counts using the voom method
# Voom transforms the raw counts into log-scaled values (log2 counts per million) and calculates weights for genes
# based on their variability. This helps stabilize the variance across different expression levels.
v <- voom(dge, design, plot = FALSE)

# Fit a linear model
# A linear model is fitted to the voom-transformed data to assess differential expression across conditions.
fit <- lmFit(v, design)

# Iterate through each contrast to perform differential expression analysis
# For each contrast (comparison between two groups), differential expression analysis is performed.
for (i in 1:nrow(contrasts)) {
  contrast <- contrasts[i, ]
  contrast_name <- paste(contrast[1], "vs", contrast[2], sep = "_")
  cat("Processing contrast:", contrast_name, "\n")
  
  # Create the contrast matrix
  # This matrix defines the specific comparison (e.g., treated vs control) for the analysis.
  cont.matrix <- makeContrasts(contrasts = paste(contrast[1], "-", contrast[2], sep = ""), levels = design)

  # Fit the contrast using the linear model
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)  # Empirical Bayes moderation to improve statistical accuracy

  # Extract results
  # The topTable function retrieves the top differentially expressed genes based on statistical significance.
  results <- topTable(fit2, number = Inf)
  results <- cbind(GeneID = rownames(results), results)
  
  # Filter results based on adjusted p-value and log fold change
  # Genes with an adjusted p-value < 0.05 and a log fold change > 1 or < -1 (FC > 2 or < 0.5) are retained.
  filtered_results <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]
  filtered_results <- filtered_results[order(-abs(filtered_results$logFC)), ]  # Sort by absolute logFC

  # Save the results to files
  output_file <- file.path(output_dir, paste0("limma_voom_results_", contrast_name, ".txt"))
  write.table(results, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

  filtered_output_file <- file.path(output_dir, paste0("limma_voom_filtered_results_", contrast_name, ".txt"))
  write.table(filtered_results, file = filtered_output_file, sep = "\t", quote = FALSE, row.names = FALSE)

  # Generate MA plot
  # An MA plot visualizes the relationship between mean expression and log fold change, highlighting DE genes.
  ma_plot_file <- file.path(output_dir, paste0("MAplot_limma_voom_", contrast_name, "_mqc.png"))
  png(ma_plot_file)
  plotMA(fit2, main = paste("MA Plot -", contrast_name))
  dev.off()

  # Generate heatmap for top 50 differentially expressed genes
  # A heatmap shows the expression of the top 50 DE genes across samples.
  top_genes <- rownames(filtered_results)[1:min(50, nrow(filtered_results))]
  heatmap_data <- v$E[top_genes, ]
  heatmap_file <- file.path(output_dir, paste0("heatmap_limma_voom_", contrast_name, "_mqc.png"))
  png(heatmap_file, width = 1000, height = 800)
  pheatmap(heatmap_data, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE, show_colnames = TRUE,
           fontsize_row = 10, fontsize_col = 10)
  dev.off()
}

