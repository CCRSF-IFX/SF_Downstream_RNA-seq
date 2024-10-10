# Load required libraries
# VennDiagram is used for generating Venn diagrams
# futile.logger is used for logging messages during script execution
library(VennDiagram)
library(futile.logger)

# Parse command-line arguments
# The first arguments are input files, followed by the output file and the log file
args <- commandArgs(trailingOnly = TRUE)
input_files <- args[1:(length(args)-2)]   # Input files containing gene lists
output_file <- args[length(args)-1]       # The output Venn diagram file path
log_file <- args[length(args)]            # The log file for logging script events

# Set logger to output to the specified log file and define logging threshold
flog.appender(appender.file(log_file))
flog.threshold(INFO)

# Log the input files and output file paths
flog.info("Input files: %s", paste(input_files, collapse = ", "))
flog.info("Output file: %s", output_file)

# Load and process input files to extract gene lists
# This step reads each input file and retrieves the "GeneID" column
gene_lists <- lapply(input_files, function(file) {
  data <- read.delim(file, header = TRUE)
  if ("GeneID" %in% colnames(data)) {
    return(data$GeneID)  # Return gene list if "GeneID" column exists
  } else {
    stop("GeneID column not found in file: ", file)  # Stop if column is missing
  }
})

# Check if at least two non-empty gene lists are provided
# A Venn diagram requires at least two sets of data
non_empty_lists <- sum(sapply(gene_lists, length) > 0)
if (non_empty_lists < 2) {
  stop("At least two non-empty gene lists are required.")
}

# Create Venn diagram
# Customize the diagram size, label positions, text size, and resolution
venn.plot <- venn.diagram(
  x = gene_lists,                      # List of gene sets for the diagram
  category.names = c("DESeq2", "edgeR", "limma-voom"),  # Category names
  filename = NULL,                     # Don't automatically save the file yet
  output = TRUE,                       # Render the diagram
  height = 400,                        # Set height of the image
  width = 400,                         # Set width of the image
  resolution = 150,                    # Set image resolution
  margin = 0.1,                        # Add margin to prevent cutting off labels
  cat.cex = 0.6,                       # Set size of category labels
  cat.pos = c(-15, 15, -180),          # Adjust positions of category labels
  cat.dist = c(0.055, 0.055, 0.085),   # Adjust distance of labels from circles
  cex = 0.8,                           # Set size of the count labels
  fontface = "bold",                   # Set font face to bold
  fontfamily = "sans"                  # Set font family to sans-serif
)

# Save the Venn diagram to the specified output file in PNG format
png(output_file, width = 400, height = 400, res = 150)  # Define PNG output size and resolution
grid.draw(venn.plot)  # Draw the diagram
dev.off()  # Close the PNG device to save the file

