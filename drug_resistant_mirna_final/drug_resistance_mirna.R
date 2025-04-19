# Load required packages
library(limma)
library(ggplot2)
library(ggrepel)
library(pheatmap)

# Set working directory to where your Series Matrix file is located
setwd("D:/DEG analysis_oisharja/drug resistance/New folder/GSE161784_series_matrix.txt")

# Read the Series Matrix file.
# The file usually contains header lines starting with "!" that we can skip.
# The parameter comment.char="!" tells R to ignore lines starting with "!".
data <- read.delim("GSE161784_series_matrix.txt", header = TRUE, sep = "\t", comment.char = "!")
# Check the structure and column names
head(data)
dim(data)  # expression data dimensions

# Typically, the first column holds probe IDs. Set them as rownames:
rownames(data) <- data[,1]
data <- data[,-1]  # remove the first column now that it's in rownames

# At this point, 'data' is a matrix of normalized expression values.
# Inspect column names to see sample IDs:
print(colnames(data))


# grouping 
group <- factor(c(rep("Sensitive", 2), rep("Resistant", 4)))

# Check that length(group) equals the number of columns in 'data'
print(length(group))
print(ncol(data))

# If your dataset has more samples, adjust the group vector accordingly.
# You can also extract group information from the file's metadata if available.

# Create the design matrix (one column per group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
print(design)

# Fit the linear model using limma
fit <- lmFit(as.matrix(data), design)


print(contrast.matrix)

# defining contrast: Sensitive vs Resistant 
contrast.matrix <- makeContrasts(Sensitive - Resistant, levels = design)

# Fit contrasts and compute moderated statistics
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract differential expression results for one contrast at a time.
# For example, for Doxorubicin vs Control:
#results_doxo <- topTable(fit2, coef = "Doxorubicin_vs_Control", number = Inf, adjust.method = "BH")
# Add a column for gene/probe names (if not already present)
#results_doxo$Gene <- rownames(results_doxo)
#head(results_doxo)

# extract DE miRNA 
results <- topTable(fit2, adjust.method = "BH", number = Inf)
head(results)

# Filter for significant miRNAs (adjusted p-value < 0.05 and absolute logFC > 1)
#sig_doxo <- subset(results_doxo, adj.P.Val < 0.05 & abs(logFC) > 1)
#head(sig_doxo)

# significant genes
sig_results <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1,]

# export 
write.csv(sig_results, file = "significant_miRNA.csv", row.names = TRUE)


# Save the results if needed


# Annotation

Annotation <- read.table("GPL19117-74051.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)

# Convert rownames (the probe IDs) to a proper column called "ProbeName"
results$ID <- rownames(results)
sig_results$ID <- rownames(sig_results)
str(results)
colnames(results)


# merge using probename

annotated_degs <- merge(results, Annotation, by.x = "ID", by.y = "ID", all.x = TRUE)
head(annotated_degs)

annotated_sig_degs <- merge(sig_results, Annotation, by.x = "ID", by.y = "ID", all.x = TRUE)


deg_gene_names <- annotated_sig_degs$miRNA_ID

print(deg_gene_names)

#export
write.csv(annotated_degs, "Annotated_DEG_Genes.csv", row.names = FALSE)

write.csv(annotated_degs, "annotated_sig_DEG_genes.csv", row.names = FALSE)


# visualisation 

# #Add a column for significance 
annotated_sig_degs$Significance <- ifelse(annotated_sig_degs$adj.P.Val < 0.05 & abs(annotated_sig_degs$logFC) > 1, 
                                          "Significant", 
                                          "Not Significant")

# volcano plot

# Create the volcano plot and label only significant genes
ggplot(annotated_sig_degs, aes(x = logFC, y = -log10(adj.P.Val), color = Significance, label = miRNA_ID)) +
  geom_point(alpha = 0.6) +
  # Only label significant points to reduce clutter
  geom_text_repel(data = subset(annotated_sig_degs, Significance == "Significant"),
                  size = 3, max.overlaps = 10) +
  scale_color_manual(values = c("red", "gray")) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10 Adjusted P-value") +
  theme_minimal()



# Load the pheatmap package (if not already loaded)
library(pheatmap)

# 1. Extract the probe IDs for significant miRNAs
sig_ids <- rownames(sig_results)

# 2. Subset the expression matrix 'data' (which was created when reading the Series Matrix file)
sig_expr <- data[sig_ids, ]

# 3. Optionally, scale the expression values by row (z-score transformation) to enhance visualization.
# This will standardize each gene's expression across samples.
sig_expr_scaled <- t(scale(t(as.matrix(sig_expr))))

# 4. Create a column annotation using your group vector.
# Here 'group' is a factor (e.g., "Sensitive" and "Resistant").
# Ensure that the rownames of col_annotation match the column names of your data.
col_annotation <- data.frame(Group = group)
rownames(col_annotation) <- colnames(data)

# 5. Generate the heatmap.
# Adjust show_rownames to TRUE if you want to see gene names (might be cluttered if there are many genes).
pheatmap(sig_expr_scaled,
         annotation_col = col_annotation,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,   # Set to FALSE to reduce clutter
         fontsize_row = 6,       # Adjust font size of row names as needed
         fontsize_col = 8,       # Adjust font size of column names as needed
         main = "Heatmap of Significant miRNAs (Resistant vs. Sensitive)")
