library(MatrixEQTL)

# Function to run Matrix_eQTL_engine without position files
# Function to run Matrix_eQTL_engine without position files
run_engine <- function(snp_file, expr_file, cov_file, output_file, p_value = 1e-5, no_fdr = FALSE, chunk_size = 2000) {
    # 1. Load Genotypes
    snps <- SlicedData$new()
    snps$fileDelimiter <- "\t"
    snps$fileOmitCharacters <- "NA"
    snps$fileSkipRows <- 1
    snps$fileSkipColumns <- 1
    snps$fileSliceSize <- chunk_size
    snps$LoadFile(snp_file)

    # 2. Load Expression/Metabolites
    gene <- SlicedData$new()
    gene$fileDelimiter <- "\t"
    gene$fileOmitCharacters <- "NA"
    gene$fileSkipRows <- 1
    gene$fileSkipColumns <- 1
    gene$fileSliceSize <- chunk_size
    gene$LoadFile(expr_file)

    # 3. Load Covariates (Optional)
    cvrt <- SlicedData$new()
    if (!is.null(cov_file) && cov_file != "") {
        cvrt$fileDelimiter <- "\t"
        cvrt$fileOmitCharacters <- "NA"
        cvrt$fileSkipRows <- 1
        cvrt$fileSkipColumns <- 1
        if (file.exists(cov_file)) {
            cvrt$LoadFile(cov_file)
        }
    }

    # 4. Run Matrix_eQTL_engine
    me <- Matrix_eQTL_engine(
        snps = snps,
        gene = gene,
        cvrt = cvrt,
        output_file_name = output_file,
        pvOutputThreshold = as.numeric(p_value),
        useModel = modelLINEAR,
        errorCovariance = numeric(),
        verbose = TRUE,
        pvalue.hist = TRUE,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = no_fdr
    )

    return(me)
}
