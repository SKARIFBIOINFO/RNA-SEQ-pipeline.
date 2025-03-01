# ------------------------------
# Install and Load DESeq2 Package
# ------------------------------

# Install DESeq2 package if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")  # Install DESeq2 from Bioconductor
library(DESeq2)  # Load the DESeq2 package

# ------------------------------
# Load Count Data
# ------------------------------

# Read the count data from a CSV file (gene expression counts from RNA-Seq)
counts_data <- read.csv("combined_gene_counts.csv")  
View(counts_data)  # Display the data for verification

# Convert decimal values to integers (DESeq2 requires integer count values)
counts_data_1 <- apply(counts_data, 2, function(x) as.integer(as.numeric(x)))

# Set gene_id column as row names (assuming the first column contains gene IDs)
rownames(counts_data_1) <- counts_data$gene_id  

# Remove the original gene_id column since it's now set as row names
counts_data_1 <- counts_data_1[, -1]
View(counts_data_1)  # Check the processed count matrix

# ------------------------------
# Define Sample Metadata
# ------------------------------

# Extract sample names from the column names of the count data
sample_data <- colnames(counts_data_1)
View(sample_data)  # Display the extracted sample names

# Define the experimental conditions (e.g., control vs. experimental groups)
# Assuming first 5 samples are controls and last 5 are experimental samples
condition <- c(rep("control", 5), rep("experiment", 5))

# Create a metadata dataframe for DESeq2 with sample names and condition labels
condition_table <- data.frame(row.names = sample_data, condition = factor(condition))
View(condition_table)  # Display metadata table

# ------------------------------
# Create DESeq2 Dataset
# ------------------------------

# Convert the count matrix into a DESeq2 dataset
dge_data <- DESeqDataSetFromMatrix(
  countData = counts_data_1,  # The gene expression count matrix
  colData = condition_table,  # Metadata describing the samples
  design = ~condition         # Define the experimental design (condition effect)
)

# ------------------------------
# Run Differential Expression Analysis
# ------------------------------

# Perform differential gene expression analysis
deg_run <- DESeq(dge_data)

# Extract results (log2 fold changes, p-values, adjusted p-values, etc.)
deg_result <- results(deg_run)
deg_result  # Print the results

# Order results by adjusted p-value (padj) for significance ranking
final_deg_result <- deg_result[order(deg_result$padj), ]

# Save the results as a CSV file
write.csv(as.data.frame(final_deg_result), "deg_glucose_vs_glycerol_1.csv", row.names = TRUE)

# View the saved results
View(read.csv("deg_glucose_vs_glycerol_1.csv"))
View("deg_glucose_vs_glycerol_1.csv")

# ------------------------------
# Plot the MA Plot for Results
# ------------------------------

# The MA plot visualizes changes in expression (log2 fold change vs. mean expression)
plotMA(deg_result)

# Generate a filtered results table using an adjusted p-value threshold of 0.05
deg_res_0.05 <- results(deg_run, alpha = 0.05)
deg_res_0.05  # Print filtered results

# Plot the MA plot for the filtered results
plotMA(deg_res_0.05)

# ------------------------------
# Create an MA Plot Using ggplot2
# ------------------------------

library(ggplot2)  # Load ggplot2 for visualization

ggplot(deg_res_0.05, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = (padj < 0.05 & abs(log2FoldChange) > 2)), alpha = 0.6, size = 1) + 
  scale_x_log10() +  # Log scale for mean expression values
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray"), 
                     labels = c("Not Significant", "Significant")) + 
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "blue") + 
  labs(title = "MA Plot of Differential Expression",
       x = "Mean Expression (baseMean)",
       y = "Log2 Fold Change",
       color = "Significance") +
  theme_minimal()  # Apply a clean theme

# Save the MA plot as an image file
ggsave("MA_plot.png", width = 6, height = 4, dpi = 300)

# ------------------------------
# Identify Upregulated and Downregulated Genes
# ------------------------------

# Extract upregulated genes (padj < 0.05 and log2FoldChange > 1)
up_regulated_genes <- deg_result[which(deg_result$padj < 0.05 & deg_result$log2FoldChange > 1), ]
up_regulated_genes  # Print upregulated genes

# Save upregulated genes to a CSV file
write.csv(as.data.frame(up_regulated_genes), "upregulated_genes.csv", row.names = TRUE)

# Extract downregulated genes (padj < 0.05 and log2FoldChange < -1)
down_regulated_genes <- deg_result[which(deg_result$padj < 0.05 & deg_result$log2FoldChange < -1), ]
down_regulated_genes  # Print downregulated genes

# Save downregulated genes to a CSV file
write.csv(as.data.frame(down_regulated_genes), "downregulated_genes.csv", row.names = TRUE)

# ------------------------------
# End of Analysis
# ------------------------------
