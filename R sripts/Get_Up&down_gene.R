# Replace with the actual column names
gene_column <- "Symbol"               # Column with gene identifiers
log2fc_column <- "log2FoldChange"     # Column with log2 fold change
padj_column <- "padj"                 # Column with adjusted p-values

# Filter for upregulated genes (log2FC > 1 and padj < 0.05)
up_genes <- deseq2_results[
  deseq2_results[[log2fc_column]] > 1 & deseq2_results[[padj_column]] < 0.05,
  gene_column
]

# Filter for downregulated genes (log2FC < -1 and padj < 0.05)
down_genes <- deseq2_results[
  deseq2_results[[log2fc_column]] < -1 & deseq2_results[[padj_column]] < 0.05,
  gene_column
]
gene_list <- unique(c(up_genes, down_genes))
head(gene_list)
