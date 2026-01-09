library(peer)

get_colname <- function(number) {
  return(paste0("PEER", number))
}

peer_function <- function(input_file, numfactors, cov_file = NULL) {
  # DEBUG: Print file head
  cat("DEBUG: Head of input file:\n")
  print(readLines(input_file, n = 5))

  # Try reading without row.names first to be safe
  expr <- read.table(input_file, header = T, check.names = F, sep = "\t", quote = "")
  if (ncol(expr) > 1 && class(expr[, 1]) != "numeric") {
    rownames(expr) <- expr[, 1]
    expr <- expr[, -1]
  }

  # Ensure numeric
  if (!all(sapply(expr, is.numeric))) {
    print("Warning: Expression matrix contains non-numeric data. Attempting conversion.")
    expr <- as.data.frame(sapply(expr, as.numeric))
  }

  # Handle NAs (Impute with row mean)
  if (any(is.na(expr))) {
    print(paste("Warning: Found", sum(is.na(expr)), "NAs in expression matrix. Imputing with row means."))
    # Row-wise imputation
    expr <- t(apply(expr, 1, function(x) {
      x[is.na(x)] <- mean(x, na.rm = TRUE)
      return(x)
    }))
    expr <- as.data.frame(expr)
  }

  print(paste("PEER Input Dims:", nrow(expr), "Genes x", ncol(expr), "Samples"))

  model <- PEER()
  PEER_setPhenoMean(model, t(as.matrix(expr)))
  PEER_setNk(model, numfactors)

  if (!is.null(cov_file) && cov_file != "") {
    cov <- read.table(cov_file, header = T, row.names = 1, check.names = F)
    # Align samples: keep only samples present in expression data, in the same order
    common_samples <- intersect(colnames(expr), rownames(cov))
    if (length(common_samples) < 2) {
      stop("Error: Too few common samples between expression and covariates.")
    }
    # Reorder covariates to match expression columns
    cov <- cov[colnames(expr), , drop = FALSE]

    # PEER expects numeric matrix for covariates
    PEER_setCovariates(model, as.matrix(cov))
    print(paste("Included", ncol(cov), "known covariates."))
  }

  PEER_update(model)

  factors <- PEER_getX(model)
  factors <- t(factors)
  colnames(factors) <- colnames(expr)


  rownames(factors) <- sapply(seq_along(1:nrow(factors)), get_colname)
  factors <- cbind(rownames(factors), factors)
  colnames(factors)[1] <- "ID"
  return(data.frame(factors))
}
