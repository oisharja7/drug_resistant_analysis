library(clusterProfiler)
library(org.Hs.eg.db)   # Human gene annotation
library(enrichplot)     # Visualization for enrichment analysis
library(DOSE)           # For additional enrichment analysis functions

# working directory 
setwd("D:/DEG analysis_oisharja/drug resistance/New folder/GSE161784_series_matrix.txt/pathway_enrichement")

# load common genes
common_deg_genes <- read.csv("Common_Genes.csv", stringsAsFactors = FALSE)


# extract gene name vector
gene_vector <- common_deg_genes$Gene

# trim extra space 
gene_vector <- trimws(gene_vector)

# convert gene symbol to entrez ID
gene_ids <- bitr(gene_vector,
                 fromType = "SYMBOL",
                 toType = "ENTREZID", 
                 OrgDb = "org.Hs.eg.db")


# extract entrex ID vector for enrichment 
entrez_list <- gene_ids$ENTREZID

# GO enrichment analysis for biological process 
ego <- enrichGO(gene = entrez_list,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2)

head(ego)


# KEGG pathway enrichment analysis 
ekegg <- enrichKEGG(gene = entrez_list,
                      organism = "hsa",
                      pvalueCutoff = 0.05)

head(ekegg)


# visualise GO enrichment result 
dotplot(ego, showCategory = 10) + ggtitle("GO Biological Process Enrichment")

# visualise KEGG enrichment result 
barplot(ekegg, showCategory = 10) + ggtitle("KEGG Pathway Enrichment")
