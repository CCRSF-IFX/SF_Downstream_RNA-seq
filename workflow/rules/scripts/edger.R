# Load necessary libraries for differential expression analysis and visualization
# edgeR: for differential gene expression analysis
# statmod: statistical modeling
# RColorBrewer: color palettes for plots
# gplots: additional plotting tools (including heatmaps)
# reshape2: data reshaping (melt)
# ggplot2: data visualization
# ggfortify: fortify ggplot with PCA plotting
library(edgeR)
library(statmod)
library(RColorBrewer)
library(gplots)
library(reshape2)
library(ggplot2)
library(ggfortify)

# Parse command-line arguments
# args[1]: path to raw counts file
# args[2]: path to coldata file (sample metadata)
# args[3]: path to contrasts file (comparisons to make)
# args[4]: output directory path
args <- commandArgs(trailingOnly = TRUE)
raw_counts_path <- args[1]
coldata_path <- args[2]
contrasts_path <- args[3]
output_dir <- args[4]

# Print paths to check correct input
cat("Raw counts file path: ", raw_counts_path, "\n")
cat("Coldata file path: ", coldata_path, "\n")
cat("Contrasts file path: ", contrasts_path, "\n")
cat("Output directory: ", output_dir, "\n")

# Create the output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read the sample metadata (coldata) and count matrix (raw counts)
sampleinfo <- read.delim(coldata_path, header = TRUE, stringsAsFactors = FALSE)
counts <- read.delim(raw_counts_path, row.names = 1)
colnames(counts) <- as.character(sampleinfo$samplename)

# Set up the condition factor based on experimental groups
# condition: factors defining the experimental groups
condition <- as.factor(sampleinfo$condition)

# Create DGEList object to hold the counts and experimental group information
y <- DGEList(counts = counts, group = condition)

# Normalize the data using TMM (Trimmed Mean of M-values)
# This method adjusts for differences in library sizes and compositional biases
y <- calcNormFactors(y, method = "TMM")

# Generate plots to visualize library size distribution and MDS plots (BCV and logFC)
png(file.path(output_dir, "libdistrib_mqc.png"), width = 800, height = 600)
barplot(y$samples$lib.size * 1e-6, main = "Library size distribution", names = colnames(y$counts), ylab = "Library size (millions)", las = 2, cex.names = 0.6)
dev.off()

png(file.path(output_dir, "MDS_bcv_mqc.png"), width = 800, height = 600)
plotMDS(y, method = "bcv", main = "MDS plot BCV")
dev.off()

png(file.path(output_dir, "MDS_logFC_mqc.png"), width = 800, height = 600)
plotMDS(y, method = "logFC", main = "MDS plot logFC", cex.lab = 0.7, cex.axis = 0.7)
dev.off()

# Estimate common and tagwise dispersions
# Dispersions help model variability in the data and are essential for accurate differential expression testing
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)

# Generate BCV plot (Biological Coefficient of Variation)
png(file.path(output_dir, "BCVplot_mqc.png"), width = 800, height = 600)
plotBCV(y, main = "BCV plot")
dev.off()

# Differential expression analysis using exactTest
# Loop through each contrast in the contrasts file
contrasts <- read.delim(contrasts_path, header = FALSE, stringsAsFactors = FALSE)
for (i in 1:nrow(contrasts)) {
  contrast <- unlist(strsplit(contrasts[i, 1], split = " "))
  
  # Perform differential expression test between the two specified conditions
  y_test <- exactTest(y, pair = c(contrast[2], contrast[1]))
  
  # Get all DEGs (differentially expressed genes)
  n <- dim(y$counts)[1]
  tt <- topTags(y_test, n = n)
  res1 <- as.data.frame(tt)
  
  # Filter for significant genes based on FDR < 0.05
  significant_genes <- res1[res1$FDR < 0.05, ]
  
  # Sort by absolute logFC (log fold change) to get top genes
  final <- significant_genes[order(abs(significant_genes$logFC), decreasing = TRUE), ]
  
  # Convert logFC to fold change (positive for upregulated, negative for downregulated)
  final$FC <- ifelse(final$logFC < 0, -1 / (2^final$logFC), 2^final$logFC)
  
  # Add GeneID column for clarity
  final <- cbind(GeneID = rownames(final), final)
  
  # Output all DEGs to a file
  output_filename <- file.path(output_dir, paste0("DEG_EdgeR_", contrast[1], "_vs_", contrast[2], ".txt"))
  write.table(final, file = output_filename, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Generated file: ", output_filename, "\n")
  
  # Filter DEGs for significant fold changes (FC >= 2 or FC <= 0.5)
  fc_threshold <- 2
  filtered_final <- final[final$FDR < 0.05 & (final$FC >= fc_threshold | final$FC <= 1 / fc_threshold), ]
  filtered_output_filename <- file.path(output_dir, paste0("Filtered_DEG_EdgeR_", contrast[1], "_vs_", contrast[2], ".txt"))
  write.table(filtered_final, file = filtered_output_filename, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Generated filtered file: ", filtered_output_filename, "\n")
  
  # Generate MA plot (mean-difference plot) to visualize differential expression
  deg1sel <- decideTestsDGE(y_test, p = 0.05, adjust = "BH")
  detags <- rownames(y)[as.logical(deg1sel)]
  if (length(detags) == 0) {
    cat("No DE tags found for contrast: ", contrast, "\n")
  } else {
    ma_plot_filename <- file.path(output_dir, paste0("Smearplot_EdgeR_", contrast[1], "_vs_", contrast[2], "_mqc.png"))
    png(ma_plot_filename, width = 800, height = 600)
    plotSmear(y_test, de.tags = detags, main = paste("Smearplot FDR<0.05 ", contrast[1], "_vs_", contrast[2]))
    abline(h = c(-2, 2), col = "blue")
    dev.off()
    cat("Generated file: ", ma_plot_filename, "\n")
  }
}

# Optional: GO and KEGG pathway analysis if Entrez gene IDs are available
if (file.exists(file.path(output_dir, "entrez_gene_ids.txt"))) {
  gene_ids <- read.delim(file.path(output_dir, "entrez_gene_ids.txt"), header = FALSE, stringsAsFactors = FALSE)
  rownames(filtered_final) <- gene_ids$V1
  
  go <- goana(y_test, species = "Mm")  # "Mm" for mouse; change to appropriate species if needed
  go_output_filename <- file.path(output_dir, paste0("GO_EdgeR_", contrast[1], "_vs_", contrast[2], ".txt"))
  write.table(topGO(go), file = go_output_filename, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Generated GO file: ", go_output_filename, "\n")
  
  kegg <- kegga(y_test, species = "Mm")
  kegg_output_filename <- file.path(output_dir, paste0("KEGG_EdgeR_", contrast[1], "_vs_", contrast[2], ".txt"))
  write.table(topKEGG(kegg), file = kegg_output_filename, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Generated KEGG file: ", kegg_output_filename, "\n")
} else {
  cat("Entrez Gene IDs file not found. Skipping GO and KEGG pathway analysis.\n")
}

# Log transformation of normalized counts
# Normalized counts are calculated in both log and raw CPM formats
ylog2 <- cpm(y, log = TRUE, normalized.lib.sizes = TRUE, prior.count = 2)
ndata <- cpm(y, log = FALSE, normalized.lib.sizes = TRUE) * 1e6

# Save normalized counts to files
write.table(ylog2, file = file.path(output_dir, "edgeR_normalized_counts_log.txt"), sep = "\t", quote = FALSE, row.names = TRUE)
write.table(ndata, file = file.path(output_dir, "edgeR_normalized_counts.txt"), sep = "\t", quote = FALSE, row.names = TRUE)
cat("Generated normalized counts files.\n")

# Generate density plot of normalized counts
png(file.path(output_dir, "HistEdgeRnormFilter_mqc.png"), width = 800, height = 600)
df.m <- melt(as.data.frame(ndata))
print(ggplot(df.m) + geom_density(aes(x = value, colour = variable)) + labs(x = NULL) + theme(legend.position = 'top', legend.key.width = unit(2, "line")) + scale_x_log10())
dev.off()

# Generate heatmap of sample-to-sample distances based on log2-normalized counts
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
distylog2 <- dist(t(ylog2))
mat <- as.matrix(distylog2)
heatmap_filename <- file.path(output_dir, paste0("edgeR_heatmaps_samplebysample_mqc_", contrast[1], "_vs_", contrast[2], ".png"))
png(heatmap_filename, width = 800, height = 600)
heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(16, 16), cexRow = 0.7, cexCol = 0.7)
dev.off()
cat("Generated heatmap: ", heatmap_filename, "\n")

# Perform PCA and visualize sample clustering based on log-transformed data
pr2 <- prcomp(t(ylog2))
dd <- cbind(t(ylog2), condition = as.character(condition))

# Standard PCA plot
pca_plot_filename <- file.path(output_dir, paste0("edgeR_prcomp_mqc_", contrast[1], "_vs_", contrast[2], ".png"))
png(pca_plot_filename, width = 800, height = 600)
plot(pr2$x[,1], pr2$x[,2], col = "red", main = "PCA plot using prcomp and Logcpm data", xlab = "PC1", ylab = "PC2")
text(pr2$x[,1], pr2$x[,2], labels = colnames(ylog2), cex = 0.7, pos = 4)
dev.off()
cat("Generated PCA plot: ", pca_plot_filename, "\n")

# PCA plot with ggplot2
autoplot_pca_filename <- file.path(output_dir, paste0("edgeR_pca_mqc_", contrast[1], "_vs_", contrast[2], ".png"))
png(autoplot_pca_filename, width = 800, height = 600)
autoplot(pr2, data = dd, colour = 'condition') + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
cat("Generated autoplot PCA: ", autoplot_pca_filename, "\n")

