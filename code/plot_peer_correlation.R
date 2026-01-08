library(pheatmap)
library(ggplot2)
library(reshape2)

# Parse Args
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop("Usage: Rscript plot_peer_correlation.R <sample_info.csv> <peer_factors.tsv> <output_prefix>")
}

sample_file <- args[1]
peer_file <- args[2]
out_prefix <- args[3]

message(paste("Reading Sample Info:", sample_file))
# Check if CSV or TSV
if (grepl("\\.tsv$", sample_file) || grepl("\\.txt$", sample_file)) {
    cov_info <- read.table(sample_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
} else {
    cov_info <- read.csv(sample_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

# Clean Sample Info Columns (strip whitespace)
colnames(cov_info) <- trimws(colnames(cov_info))

# Ensure 'sample' column exists
if (!"sample" %in% colnames(cov_info)) {
    # Try to detect sample column: usually the first one or named 'ID'
    if ("ID" %in% colnames(cov_info)) {
        colnames(cov_info)[which(colnames(cov_info) == "ID")] <- "sample"
    } else {
        message("Warning: 'sample' column not found. Using the first column as sample ID.")
        colnames(cov_info)[1] <- "sample"
    }
}

# Fix: Robust cleanup
# 1. Trim whitespace from ALL character columns
cov_info[] <- lapply(cov_info, function(x) if (is.character(x)) trimws(x) else x)

# 2. Filter valid samples
cov_info <- cov_info[!is.na(cov_info$sample) & cov_info$sample != "", , drop = FALSE]

# 3. Handle Duplicates
if (any(duplicated(cov_info$sample))) {
    dups <- cov_info$sample[duplicated(cov_info$sample)]
    message(paste("Warning: Removing", length(dups), "duplicate sample entries (keeping first):", paste(head(dups), collapse = ", ")))
    cov_info <- cov_info[!duplicated(cov_info$sample), ]
}

rownames(cov_info) <- cov_info$sample

message(paste("Reading PEER Factors:", peer_file))
peer <- read.table(peer_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
# PEER file: Rows=Factors, Cols=Samples. Need to transpose.
peer_t <- t(peer)

# Intersect Samples
common_samples <- intersect(rownames(peer_t), rownames(cov_info))
message(paste("Matched samples:", length(common_samples)))

if (length(common_samples) == 0) {
    stop("Error: No matching samples found between Sample Info and PEER file.")
}

peer_dat <- peer_t[common_samples, , drop = FALSE]
cov_dat <- cov_info[common_samples, , drop = FALSE]

# Auto-detect Target Columns (All except Sample ID)
# Use all columns except 'sample'
target_cols <- setdiff(colnames(cov_info), "sample")

message(paste("Selected Covariates (All except ID):", paste(target_cols, collapse = ", ")))

if (length(target_cols) == 0) {
    stop("No covariates found in sample info file.")
}

bio_cov <- cov_dat[, target_cols, drop = FALSE]

# Ensure everything is numeric
# Strategy:
# 1. Try numeric conversion.
# 2. If it results in significant NAs, treat as factor and convert levels to numeric.
bio_cov[] <- lapply(bio_cov, function(x) {
    if (is.numeric(x)) {
        return(x)
    }

    num_x <- suppressWarnings(as.numeric(as.character(x)))
    na_orig <- sum(is.na(x))
    na_new <- sum(is.na(num_x))

    # If conversion drastically increases NAs (>50% of data lost), fallback to factor encoding
    if (na_new > na_orig && (na_new / length(x)) > 0.5) {
        message("Note: Converting categorical column to factor levels.")
        return(as.numeric(as.factor(x)))
    }
    return(num_x)
})
rownames(bio_cov) <- common_samples

# Remove constant columns (Variance = 0)
vars <- apply(bio_cov, 2, var, na.rm = TRUE)
valid_cols <- vars > 1e-6 & !is.na(vars)
bio_cov <- bio_cov[, valid_cols, drop = FALSE]

if (ncol(bio_cov) == 0) {
    stop("All covariates had near-zero variance. Cannot compute correlation.")
}

message("Calculating stats...")

# Calculate Correlation (PEER vs Bio)
# Result: Rows = PEER Factors, Cols = Bio Covariates
cor_mat <- cor(peer_dat, bio_cov, use = "pairwise.complete.obs", method = "pearson")
p_mat <- matrix(NA, nrow = nrow(cor_mat), ncol = ncol(cor_mat))
rownames(p_mat) <- rownames(cor_mat)
colnames(p_mat) <- colnames(cor_mat)

# Calculate P-values
for (i in rownames(cor_mat)) {
    for (j in colnames(cor_mat)) {
        res <- cor.test(peer_dat[, i], bio_cov[, j])
        p_mat[i, j] <- res$p.value
    }
}

# Save Raw Data
out_csv <- paste0(out_prefix, "_correlation.csv")
write.csv(cor_mat, out_csv)
message(paste("Saved correlation matrix:", out_csv))

# Plot Heatmap
out_pdf <- paste0(out_prefix, "_heatmap.pdf")
pdf(out_pdf, width = max(8, ncol(cor_mat) * 0.8), height = max(8, nrow(cor_mat) * 0.4))

# Mark significant correlations with *
display_txt <- matrix(ifelse(p_mat < 0.05, "*", ""), nrow(p_mat))

pheatmap(cor_mat,
    display_numbers = display_txt,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fontsize_number = 15,
    main = "Correlation: PEER Factors vs Biological Covariates (* p<0.05)",
    breaks = seq(-1, 1, length.out = 100),
    color = colorRampPalette(c("blue", "white", "red"))(100)
)

dev.off()
message(paste("Saved heatmap:", out_pdf))
