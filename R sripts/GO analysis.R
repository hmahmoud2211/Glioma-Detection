# Assuming you already have the gene symbols and Entrez IDs in `entrez_ids`
# Extract Entrez IDs
gene_list_entrez <- entrez_ids$entrezgene_id

# Load necessary libraries
library(clusterProfiler)
library(org.Hs.eg.db)  # or another organism's database

# Perform GO enrichment analysis
go_enrich <- enrichGO(gene = gene_list_entrez,
                      OrgDb = org.Hs.eg.db,
                      ont = "ALL",  # You can use "BP", "MF", "CC" for specific ontologies
                      pAdjustMethod = "BH",  # Adjust p-values using Benjamini-Hochberg method
                      qvalueCutoff = 0.05)   # Adjust p-values for false discovery rate

# View results
summary(go_enrich)
# Visualizing the GO enrichment results
dotplot(go_enrich, showCategory = 20)  # Adjust the number of categories as needed

# Optional: You can also use a bar plot for visualization
barplot(go_enrich, showCategory = 20)

