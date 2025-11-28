library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
entrez_ids <- getBM(
  filters = "hgnc_symbol",                  # Use the correct filter for gene symbols
  values = gene_list,                      # List of gene symbols to convert
  attributes = c("hgnc_symbol", "entrezgene_id"),  # Retrieve gene symbols and Entrez IDs
  mart = ensembl                           # Database connection
)
write.csv(entrez_ids, "C:\\Users\\user\\Downloads\\Bioinfo\\The output of R\\gene_symbol_to_entrez.csv", row.names = FALSE)


