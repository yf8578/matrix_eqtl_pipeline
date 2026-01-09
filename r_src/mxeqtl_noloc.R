library(MatrixEQTL)

# Function to run Matrix_eQTL_engine without position files
# Function to run Matrix_eQTL_engine without position files
# Function to run Matrix_eQTL_engine without position files
run_engine <- function(snp_file, expr_file, cov_file, output_file, p_value = 1e-5,
                       no_fdr = FALSE, chunk_size = 2000,
                       sep = "\t", missing = "NA", header = TRUE, rownames = TRUE,
                       model = "linear") {
    # Model Selection
    useModel <- modelLINEAR
    model_lower <- tolower(model)
    if (model_lower == "anova") {
        useModel <- modelANOVA
    } else if (model_lower == "linear_cross") {
        useModel <- modelLINEAR_CROSS
    }

    # 1. Load Genotypes
    snps <- SlicedData$new()
    snps$fileDelimiter <- sep
    snps$fileOmitCharacters <- missing
    snps$fileSkipRows <- if (header) 1 else 0
    snps$fileSkipColumns <- if (rownames) 1 else 0
    snps$fileSliceSize <- chunk_size
    snps$LoadFile(snp_file)

    # 2. Load Expression/Metabolites
    gene <- SlicedData$new()
    gene$fileDelimiter <- sep
    gene$fileOmitCharacters <- missing
    gene$fileSkipRows <- if (header) 1 else 0
    gene$fileSkipColumns <- if (rownames) 1 else 0
    gene$fileSliceSize <- chunk_size
    gene$LoadFile(expr_file)

    # 3. Load Covariates (Optional)
    cvrt <- SlicedData$new()
    if (!is.null(cov_file) && cov_file != "") {
        cvrt$fileDelimiter <- sep
        cvrt$fileOmitCharacters <- missing
        cvrt$fileSkipRows <- if (header) 1 else 0
        cvrt$fileSkipColumns <- if (rownames) 1 else 0
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
        useModel = useModel,
        errorCovariance = numeric(),
        verbose = TRUE,
        pvalue.hist = TRUE,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = no_fdr
    )

    return(me)
}
