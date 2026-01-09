#!/usr/bin/env Rscript

# Script to fetch gene locations with Strand information from Ensembl BioMart
# Usage: Rscript fetch_biomart_positions.R [output_file] [dataset]
# Defaults: output=gene_positions_biomart.txt, dataset=hsapiens_gene_ensembl

args <- commandArgs(trailingOnly = TRUE)
output_file <- if (length(args) >= 1) args[1] else "gene_positions_biomart.txt"
dataset_name <- if (length(args) >= 2) args[2] else "hsapiens_gene_ensembl"

cat(">>> Loading biomaRt...\n")
if (!require("biomaRt")) {
    stop("biomaRt package is not installed. Please install it with: BiocManager::install('biomaRt')")
}

# Connect to Ensembl
cat(paste(">>> Connecting to Ensembl (dataset:", dataset_name, ")...\n"))
# Use tryCatch to handle potential connection issues (common with BioMart)
ensembl <- tryCatch(
    {
        useEnsembl(biomart = "genes", dataset = dataset_name)
    },
    error = function(e) {
        cat("Error connecting to main Ensembl site. Trying mirrors...\n")
        # Try Asia mirror if main fails (good for users in China/Asia)
        useEnsembl(biomart = "genes", dataset = dataset_name, mirror = "asia")
    }
)

# Define attributes to fetch
# ensembl_gene_id: Gene ID matches typical expression matrices
# chromosome_name: Chr
# start_position: Left
# end_position: Right
# strand: Strand (1/-1)
# external_gene_name: Gene Symbol (for reference)

attributes <- c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name")

cat(">>> Fetching data (this may take a minute)...\n")
genes <- getBM(attributes = attributes, mart = ensembl)

# Filter valid chromosomes (1-22, X, Y, MT) to avoid patches
valid_chrs <- c(as.character(1:22), "X", "Y", "MT")
genes <- genes[genes$chromosome_name %in% valid_chrs, ]

# Rename columns to match pipeline requirements
# Pipeline expects: geneid, chr, left, right, strand
colnames(genes) <- c("geneid", "chr", "left", "right", "strand", "gene_symbol")

# Sort by Chr and Pos
genes$chr <- factor(genes$chr, levels = valid_chrs)
genes <- genes[order(genes$chr, genes$left), ]

# Write to file
cat(paste(">>> Writing", nrow(genes), "genes to", output_file, "...\n"))
write.table(genes, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Done!\n")
