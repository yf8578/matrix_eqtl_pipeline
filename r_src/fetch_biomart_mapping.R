#!/usr/bin/env Rscript
# Fetch ID Mapping from BioMart
# Usage: Rscript fetch_biomart_mapping.R [output_file] [source_attr] [target_attr] [dataset]
# Example: Rscript fetch.R map.txt uniprot_gn_id ensembl_gene_id hsapiens_gene_ensembl

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop("Usage: fetch_biomart_mapping.R output_file source_attr target_attr [dataset]")
}

output_file <- args[1]
source_attr <- args[2]
target_attr <- args[3]
dataset_name <- if (length(args) >= 4) args[4] else "hsapiens_gene_ensembl"

cat(">>> Loading biomaRt...\n")
if (!require("biomaRt", quietly = TRUE)) {
    stop("biomaRt not installed.")
}

cat(paste(">>> Connecting to Ensembl (", dataset_name, ")...\n"))
ensembl <- tryCatch(
    {
        useEnsembl(biomart = "genes", dataset = dataset_name)
    },
    error = function(e) {
        cat("Error connecting to main Ensembl. Trying Asia mirror...\n")
        useEnsembl(biomart = "genes", dataset = dataset_name, mirror = "asia")
    }
)

attributes <- c(source_attr, target_attr)

cat(paste(">>> Fetching mapping:", source_attr, "->", target_attr, "...\n"))
# distinct=TRUE removes duplicates
mapping <- getBM(attributes = attributes, mart = ensembl, uniqueRows = TRUE)

# Filter empty strings (BioMart often returns empty strings for missing IDs)
mapping <- mapping[mapping[[source_attr]] != "" & mapping[[target_attr]] != "", ]
mapping <- na.omit(mapping)

cat(paste(">>> Writing", nrow(mapping), "mappings to", output_file, "...\n"))
write.table(mapping, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
