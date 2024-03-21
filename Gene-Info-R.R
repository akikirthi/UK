# Import gene data focusing on essential columns
gene_details <- fread("G:\\Homo_sapiens.gene_info\\Homo_sapiens.gene_info", select = c(2, 3, 5))

# Mapping gene symbols to Entrez IDs
id_map <- setNames(gene_details$V2, gene_details$V3)
gene_details$V5 <- as.character(gene_details$V5)

# Handling synonyms
split_synonyms <- strsplit(gene_details$V5, "\\|")
for (i in seq_along(split_synonyms)) {
  synonyms_list <- split_synonyms[[i]][split_synonyms[[i]] != ""]
  id_map[synonyms_list] <- gene_details$V2[i]
}

# Processing a GMT file
gmt_content <- readLines("G:\\h.all.v2023.1.Hs.symbols.gmt")

# Substituting gene symbols with IDs
updated_gmt <- sapply(gmt_content, function(line) {
  parts <- strsplit(line, "\t")[[1]]
  pathway_elements <- parts[1:2]
  symbols <- parts[-(1:2)]
  ids <- sapply(symbols, function(symbol) id_map[symbol])
  c(pathway_elements, ids)
})

# Saving the updated GMT file
writeLines(sapply(updated_gmt, function(item) paste(item, collapse = "\t")), "G:\\updated_gmt.gmt")


-------------->

  # Loading the protein sequences
  sequence_data <- readAAStringSet("G:\\NC_000913.faa\\NC_000913.faa")

# Computing the average sequence length
average_length <- mean(width(sequence_data))
print(average_length)


---------------->
  


# Load gene info
gene_info <- fread("G:\\Homo_sapiens.gene_info\\Homo_sapiens.gene_info", select = c(3, 7))



# Assuming gene_info has been correctly loaded and contains the necessary chromosome information
# Ensure the column names are known
colnames(gene_info) <- c("Symbol", "Chromosome") # Adjust these names based on the actual content

# Filter out ambiguous chromosome values
clean_gene_info <- gene_info[!grepl("\\|", gene_info$Chromosome),]

# Aggregate gene counts by chromosome
gene_counts <- clean_gene_info[, .N, by = .(Chromosome)]
gene_counts <- gene_counts[order(Chromosome),]

# Plotting the distribution of genes across chromosomes
p <- ggplot(gene_counts, aes(x = Chromosome, y = N)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Number of Genes Across Chromosomes", x = "Chromosome", y = "Gene Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Improve readability of chromosome labels

# Display the plot
print(p)

# Optionally, save the plot to a PDF file
ggsave("genes_per_chromosome.pdf", plot = p, path = "G:\\", width = 11, height = 8, units = "in")
