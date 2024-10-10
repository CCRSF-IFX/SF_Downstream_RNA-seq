# Load necessary libraries
# These libraries are used for differential gene expression analysis (DESeq2), plotting (ggplot2, pheatmap), 
# reading/writing Excel files (openxlsx), reshaping data (reshape2), heatmaps (gplots, pheatmap), and data manipulation (SummarizedExperiment).
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(openxlsx)
library(reshape2)
library(gplots)
library(SummarizedExperiment)

# Parse command-line arguments
# These variables will hold the file paths for the raw counts, sample metadata (coldata), contrasts, and output directory.
args <- commandArgs(trailingOnly = TRUE)
raw_counts_path <- args[1]
coldata_path <- args[2]
contrasts_path <- args[3]
output_dir <- args[4]

# Display the file paths for verification
cat("Raw counts file path: ", raw_counts_path, "\n")
cat("Coldata file path: ", coldata_path, "\n")
cat("Contrasts file path: ", contrasts_path, "\n")
cat("Output directory: ", output_dir, "\n")

# Create the output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read the sample metadata
# This file contains information about the samples, such as condition and batch, which will be used for the analysis design.
sampleinfo <- read.delim(coldata_path, header = TRUE, stringsAsFactors = FALSE)
sampleinfo$condition <- factor(sampleinfo$condition)
sampleinfo$batch <- factor(sampleinfo$batch)
rownames(sampleinfo) <- sampleinfo$samplename
sampleinfo$samplename <- NULL

# Print the condition levels to check that the metadata is correct
print(levels(sampleinfo$condition))

# Read the raw counts data
# The raw gene expression data (count matrix) is loaded, and the counts are rounded for consistency.
x <- read.delim(raw_counts_path, row.names = 1, check.names = FALSE)
x <- round(x)

# Create a DESeq2 dataset
# This prepares the dataset for differential expression analysis by specifying the counts and sample info and setting the design formula.
ddsHTSeq <- DESeqDataSetFromMatrix(countData = x, colData = sampleinfo, design = ~ condition)
dds <- DESeq(ddsHTSeq)

# Get normalized counts
# The raw counts are normalized by library size to ensure comparability across samples, and the results are saved to a file.
ndata <- as.data.frame(counts(dds, normalized = TRUE))
write.table(ndata, file = file.path(output_dir, "Deseq2_normalized_counts.txt"), sep = "\t", col.names = NA)
cat("Generated file: ", file.path(output_dir, "Deseq2_normalized_counts.txt"), "\n")

# Plot histogram of normalized counts
# This step generates a density plot showing the distribution of normalized expression values across samples.
png(file.path(output_dir, "HistDesq2normFilter_mqc.png"))
df.m <- melt(as.data.frame(ndata))
print(ggplot(df.m) + geom_density(aes(x = value, colour = variable)) + labs(x = NULL) + theme(legend.position = 'top', legend.key.width = unit(2, "line")) + scale_x_log10())
dev.off()
cat("Generated file: ", file.path(output_dir, "HistDesq2normFilter_mqc.png"), "\n")

# Read and process the contrasts file
# Contrasts specify the pairs of conditions to compare, which are read from the contrasts file.
contrasts <- readLines(contrasts_path)
contrast_pairs <- strsplit(contrasts, " ")

# First loop for DEG files and MA plots
# Loop through each contrast pair for DEG analysis and MA plot generation.
# For each contrast, differential expression is calculated, results are saved, and an MA plot is generated.
for (pair in contrast_pairs) {
  if (length(pair) == 2) {
    contrast_name <- paste(pair[1], "vs", pair[2], sep = "_")
    cat("Comparing:", pair[1], "vs", pair[2], "\n")
    
    # Obtain differential expression results for the given contrast
    res <- results(dds, contrast = c("condition", pair[1], pair[2]))

    # Remove rows with NA values in pvalue or log2FoldChange
    res <- res[!is.na(res$pvalue) & !is.na(res$log2FoldChange), ]

    # Filter results based on significance threshold (p-value < 0.05)
    res <- res[res$pvalue < 0.05, ]

    # Order results by absolute log2 fold change (most significant fold changes first)
    res <- res[order(abs(res$log2FoldChange), decreasing = TRUE), ]

    # Add FoldChange and Regulation columns (classify genes as up or down-regulated)
    res$FoldChange <- 2^res$log2FoldChange
    res$Regulation <- ifelse(res$log2FoldChange > 0, "Up", "Down")

    # Save results to a file
    res_df <- as.data.frame(res)
    res_df <- cbind(GeneID = rownames(res_df), res_df)
    
    output_filename <- file.path(output_dir, paste0("DEG_", contrast_name, ".txt"))
    write.table(res_df, file = output_filename, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Generated file: ", output_filename, "\n")

    # Generate MA plot for the contrast
    x <- res$log2FoldChange[which(!is.na(res$log2FoldChange))]
    png(file.path(output_dir, paste("MAplot_", contrast_name, "_mqc.png", sep = "")))
    plotMA(res, ylim = range(x), main = paste("MAplot_", contrast_name, sep = ""))
    dev.off()
    cat("Generated file: ", file.path(output_dir, paste("MAplot_", contrast_name, "_mqc.png", sep = "")), "\n")
  } else {
    warning("Invalid contrast pair:", pair)
  }
}

# Loop through contrasts again for filtered DEG analysis and heatmap generation
for (pair in contrast_pairs) {
  if (length(pair) == 2) {
    contrast_name <- paste(pair[1], "vs", pair[2], sep = "_")
    cat("Filtering results for:", pair[1], "vs", pair[2], "\n")

    # Obtain differential expression results for the given contrast
    res <- results(dds, contrast = c("condition", pair[1], pair[2]))

    # Remove rows with NA values in pvalue or log2FoldChange
    res <- res[!is.na(res$pvalue) & !is.na(res$log2FoldChange), ]

    # Filter based on p-value < 0.05 and apply fold change threshold (Fold Change >= 2 or Fold Change <= 0.5)
    filtered_res <- res[res$pvalue < 0.05 & (abs(res$log2FoldChange) >= log2(2)), ]

    # Save filtered results to a file
    filtered_res_df <- as.data.frame(filtered_res)
    filtered_res_df <- cbind(GeneID = rownames(filtered_res_df), filtered_res_df)

    output_filename <- file.path(output_dir, paste0("DEG_", contrast_name, "_filtered.txt"))
    write.table(filtered_res_df, file = output_filename, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Generated file: ", output_filename, "\n")

    # Generate heatmap for top differentially expressed genes
    top30_genes <- rownames(head(filtered_res_df, 30))
    vsd <- vst(dds, blind = FALSE)

    ordered_counts <- assay(vsd)[top30_genes, , drop = FALSE]
    relevant_samples <- colnames(ordered_counts)[sampleinfo$condition %in% pair]
    ordered_counts <- ordered_counts[, relevant_samples, drop = FALSE]

    # Adjust height and width dynamically based on the number of genes and samples
    num_genes <- nrow(ordered_counts)
    num_samples <- ncol(ordered_counts)
    plot_height <- max(8, num_genes * 0.5)
    plot_width <- max(10, num_samples * 2)

    my_palette <- colorRampPalette(c("green", "white", "red"))(n = 100)
    png(file.path(output_dir, paste("heatmap_", contrast_name, "_mqc.png", sep = "")), width = plot_width, height = plot_height, units = "in", res = 300)
    pheatmap(ordered_counts, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE, color = my_palette, show_rownames = TRUE, fontsize_row = 12, fontsize_col = 12)
    dev.off()
    cat("Generated file: ", file.path(output_dir, paste("heatmap_", contrast_name, "_mqc.png", sep = "")), "\n")
  } else {
    warning("Invalid contrast pair:", pair)
  }
}

