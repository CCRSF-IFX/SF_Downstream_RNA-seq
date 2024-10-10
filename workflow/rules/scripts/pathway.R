# Set library path to ensure the correct R libraries are used
.libPaths("enter_lib_path_here")

# Load required libraries for pathway enrichment and visualization
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(pathview)
library(knitr)
library(png)
library(grid)
library(AnnotationDbi)
library(DOSE)

# Parse command-line arguments passed from Snakemake
args <- commandArgs(trailingOnly = TRUE)
deseq2_file <- args[1]           # Input DESeq2 file
contrast <- args[2]              # Contrast name (comparison)
output_prefix <- args[3]         # Output directory path
genome <- args[4]                # Genome type (human, mouse, monkey)
organism_db <- args[5]           # Organism database (e.g., org.Hs.eg.db)
kegg_organism <- args[6]         # KEGG organism code (e.g., hsa, mmu)
pathway_id <- args[7]            # Specific pathway ID (e.g., hsa04130)

# Ensure the output directory exists, create it if necessary
if (!dir.exists(output_prefix)) {
    dir.create(output_prefix, recursive = TRUE)
}

# Load the appropriate organism database based on the genome input
if (genome == "human") {
    library(org.Hs.eg.db)
    organism_db <- org.Hs.eg.db
} else if (genome == "mouse") {
    library(org.Mm.eg.db)
    organism_db <- org.Mm.eg.db
} else if (genome == "monkey") {
    library(org.Mmu.eg.db)
    organism_db <- org.Mmu.eg.db
} else {
    stop("Unsupported genome")
}

# Read in the DESeq2 results file
deseq2 <- read.table(deseq2_file, header = TRUE, row.names = 1, sep = "\t")

# Clean gene IDs by removing any suffix (e.g., versions after the dot)
clean_gene_ids <- function(gene_names) {
    sapply(gene_names, function(x) {
        sub("\\..*", "", x)
    })
}
cleaned_ids <- clean_gene_ids(rownames(deseq2))
unique_cleaned_ids <- make.unique(cleaned_ids)

# Filter the DESeq2 results to include only valid ENSEMBL IDs
valid_ensembl_ids <- keys(organism_db, keytype = "ENSEMBL")
valid_ids <- unique_cleaned_ids[unique_cleaned_ids %in% valid_ensembl_ids]
deseq2_filtered <- deseq2[unique_cleaned_ids %in% valid_ids, ]
rownames(deseq2_filtered) <- make.unique(valid_ids)  # Ensure row names are unique

# Prepare the gene list for enrichment analysis
deseq2_gene_list <- deseq2_filtered$log2FoldChange
names(deseq2_gene_list) <- rownames(deseq2_filtered)
deseq2_gene_list <- sort(na.omit(deseq2_gene_list), decreasing = TRUE)

# Run GSEA (Gene Set Enrichment Analysis) on the DESeq2 gene list
gse <- gseGO(geneList = deseq2_gene_list, 
             ont = "ALL", 
             keyType = "ENSEMBL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism_db, 
             pAdjustMethod = "none")

# Determine the number of significant categories based on p.adjust < 0.05
n_categories <- length(gse@result$Description[gse@result$p.adjust < 0.05])

# Adjust the height of the plot based on the number of categories
base_height <- 0.5
plot_height <- n_categories * base_height
plot_height <- min(max(plot_height, 6), 15)  # Set minimum and maximum height

# Generate dot plot for pathway enrichment
dotplot_gg <- dotplot(gse, showCategory = min(15, n_categories), split = ".sign") +  
              facet_grid(. ~ .sign) +
              theme(axis.text.y = element_text(size = 8, hjust = 1, vjust = 1)) +
              theme(axis.text.x = element_text(size = 10)) +
              theme(axis.title = element_text(size = 12)) +
              theme(strip.text = element_text(size = 12)) +
              theme(plot.margin = margin(10, 10, 10, 10))

# Save the dot plot with dynamic height adjustment
ggsave(file.path(output_prefix, paste0("pathway_analysis_deseq2_", contrast, "_dotplot_mqc.png")),
       plot = dotplot_gg, width = 10, height = plot_height, limitsize = FALSE)

# Save additional GSEA visualizations (emapplot, cnetplot, ridgeplot)
ggsave(file.path(output_prefix, paste0("pathway_analysis_deseq2_", contrast, "_emapplot_mqc.png")),
       plot = emapplot(pairwise_termsim(gse), showCategory = 10))
ggsave(file.path(output_prefix, paste0("pathway_analysis_deseq2_", contrast, "_cnetplot_mqc.png")),
       plot = cnetplot(gse, categorySize = "pvalue", color.params = list(foldChange = deseq2_gene_list), showCategory = 3))

# Adjust ridge plot settings and save
ridgeplot_gg <- ridgeplot(gse) + 
                labs(x = "enrichment distribution") +
                theme(axis.text.y = element_text(size = 8, hjust = 1, vjust = 0.5)) +
                theme(axis.text.x = element_text(size = 10)) +
                theme(axis.title = element_text(size = 12)) +
                theme(strip.text = element_text(size = 12)) +
                theme(plot.margin = margin(10, 10, 10, 10))

plot_height <- max(min(n_categories * 0.4, 15), 6)  # Adjust height for ridge plot
ggsave(file.path(output_prefix, paste0("pathway_analysis_deseq2_", contrast, "_ridgeplot_mqc.png")),
       plot = ridgeplot_gg, width = 10, height = plot_height, limitsize = FALSE)

# GSEA plot for the first category
ggsave(file.path(output_prefix, paste0("pathway_analysis_deseq2_", contrast, "_gseaplot_mqc.png")),
       plot = gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1))

# Convert ENSEMBL IDs to ENTREZ IDs for KEGG enrichment analysis
ids <- bitr(names(deseq2_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = organism_db)
df2 <- merge(deseq2_filtered, ids, by.x = "row.names", by.y = "ENSEMBL", all.x = TRUE)
rownames(df2) <- make.unique(df2$Row.names)
df2 <- df2[, -1]

# Prepare gene list for KEGG enrichment
kegg_gene_list <- na.omit(setNames(df2$log2FoldChange, df2$ENTREZID))
kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)

# Run KEGG enrichment
kk2 <- gseKEGG(geneList = kegg_gene_list,
               organism = kegg_organism,
               minGSSize = 3,
               maxGSSize = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType = "ncbi-geneid")

# Save KEGG visualizations
ggsave(file.path(output_prefix, paste0("pathway_analysis_kegg_", contrast, "_dotplot_mqc.png")),
       plot = dotplot(kk2, showCategory = 10, title = "Enriched Pathways", split = ".sign"))
ggsave(file.path(output_prefix, paste0("pathway_analysis_kegg_", contrast, "_emapplot_mqc.png")),
       plot = emapplot(pairwise_termsim(kk2), showCategory = 10))
ggsave(file.path(output_prefix, paste0("pathway_analysis_kegg_", contrast, "_cnetplot_mqc.png")),
       plot = cnetplot(kk2, categorySize = "pvalue", color.params = list(foldChange = kegg_gene_list), showCategory = 3))

# Ridge plot for KEGG
kegg_ridgeplot <- ridgeplot(kk2) + 
                  labs(x = "enrichment distribution") +
                  theme(axis.text.y = element_text(size = 8))

# Save KEGG ridge plot
ggsave(file.path(output_prefix, paste0("pathway_analysis_kegg_", contrast, "_ridgeplot_mqc.png")),
       plot = kegg_ridgeplot, width = 10, height = plot_height, limitsize = FALSE)

# GSEA plot for KEGG analysis
ggsave(file.path(output_prefix, paste0("pathway_analysis_kegg_", contrast, "_gseaplot_mqc.png")),
       plot = gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1))
