# Load necessary library
library(dplyr)

# Step 1: Read the DESeq2 output CSV file
deseq2_results <- read.csv("C:\\Users\\user\\Downloads\\Bioinfo\\Desqoutput\\Output\\Output_Desq2.csv")  # Replace with your CSV file name

# Step 2: Filter significant DEGs
upregulated_genes <- deseq_output %>%
  filter(log2FoldChange > 1 & padj < 0.05)

downregulated_genes <- deseq_output %>%
  filter(log2FoldChange < -1 & padj < 0.05)

write.table(upregulated_genes$GeneID, "C:\\Users\\user\\Downloads\\Bioinfo\\Desqoutput\\up_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(downregulated_genes$GeneID, "C:\\Users\\user\\Downloads\\Bioinfo\\Desqoutput\\down_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Optional: Print summary
cat("Number of upregulated genes:", nrow(upregulated_genes), "\n")
cat("Number of downregulated genes:", nrow(downregulated_genes), "\n")

