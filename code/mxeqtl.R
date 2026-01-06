library(MatrixEQTL)

# ... (previous code)
toSkip <- function(exists) {
    if (exists) {
        return(1)
    } else {
        return(0)
    }
}

getGenotypes <- function(sep, missing, header, rownames, snp_file, chunk_size) {
    snps <- SlicedData$new()
    snps$fileDelimiter <- sep
    snps$fileOmitCharacters <- missing
    snps$fileSkipRows <- if (header) 1 else 0
    snps$fileSkipColumns <- if (rownames) 1 else 0
    snps$fileSliceSize <- chunk_size
    snps$LoadFile(snp_file)
    return(snps)
}

getExpression <- function(sep, missing, header, rownames, expr_file, chunk_size) {
    gene <- SlicedData$new()
    gene$fileDelimiter <- sep
    gene$fileOmitCharacters <- missing
    gene$fileSkipRows <- if (header) 1 else 0
    gene$fileSkipColumns <- if (rownames) 1 else 0
    gene$fileSliceSize <- chunk_size
    gene$LoadFile(expr_file)
    return(gene)
}

getCovariates <- function(sep, missing, header, rownames, cov_file) {
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
    return(cvrt)
}

mafFilter <- function(snps, MAF) {
    maf <- snps$RowCheckFun(function(row) {
        m <- mean(row, na.rm = TRUE) / 2
        min(m, 1 - m)
    })
    snps$RowReorder(maf >= MAF)
    return(snps)
}

# ... (skipping unchanged code)

setModel <- function(model) {
    model_lower <- tolower(model)
    if (model_lower == "linear") {
        return(modelLINEAR)
    } else if (model_lower == "anova") {
        return(modelANOVA)
    } else if (model_lower == "linear_cross") {
        return(modelLINEAR_CROSS)
    } else {
        stop("Model must be linear, anova or linear_cross")
    }
}

getCovFileName <- function(covariates) {
    if (covariates != "") {
        return(covariates)
    } else {
        return(character())
    }
}

mxeqtl <-
    function(snp_file, snp_location, expr_file, expr_location, cis_output_file,
             cis_pval, covariates = "", trans_output_file = "", trans_pval = 0,
             model = "linear", MAF = 0, cis_dist = 1e6, qq = "", missing = "NA", sep = "\t",
             header = TRUE, rownames = TRUE, chunk_size = 2000) {
        # Matrix eQTL function based on the sample code by Andrey A. Shabalin
        # http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/

        useModel <- setModel(model) # nocov start

        # Covariates file name

        covariates_file_name <- getCovFileName(covariates)

        # Output file name
        output_file_name_cis <- cis_output_file
        if (trans_output_file != "") {
            output_file_name_tra <- trans_output_file
        } else {
            output_file_name_tra <- "onlyTRANSresults.txt"
        }

        # TODO: Create option for errorCovariance
        errorCovariance <- numeric()

        snps <- getGenotypes(sep, missing, header, rownames, snp_file, chunk_size)

        if (MAF > 0) {
            snps <- mafFilter(snps, MAF)
            snps$SaveFile("meQTL_filtered_input")
        }

        gene <- getExpression(sep, missing, header, rownames, expr_file, chunk_size)

        cvrt <- getCovariates(sep, missing, header, rownames, covariates_file_name)

        # Load the genotype and expression positions
        snpspos <- read.table(snp_location, header = TRUE, stringsAsFactors = FALSE)
        genepos <- read.table(expr_location, header = TRUE, stringsAsFactors = FALSE)

        # Call MatrixEQTL function
        me <- Matrix_eQTL_main(
            snps = snps,
            gene = gene,
            cvrt = cvrt,
            output_file_name = output_file_name_tra,
            pvOutputThreshold = trans_pval,
            useModel = useModel,
            errorCovariance = errorCovariance,
            verbose = TRUE,
            output_file_name.cis = output_file_name_cis,
            pvOutputThreshold.cis = cis_pval,
            snpspos = snpspos,
            genepos = genepos,
            cisDist = cis_dist,
            pvalue.hist = "qqplot",
            min.pv.by.genesnp = TRUE,
            noFDRsaveMemory = FALSE
        )

        cat("Analysis done in: ", me$time.in.sec, " seconds", "\n")
        cat("Detected ", me$cis$neqtls, " local eQTLs:", "\n")
        cat("Detected ", me$trans$neqtls, " distant eQTLs:", "\n")

        ## Plot the Q-Q plot of local and distant p-values
        if (qq != "") {
            pdf(qq)
            plot(me)
            dev.off()
        } else {
            plot(me)
        }

        return(me) # nocov end
    }
