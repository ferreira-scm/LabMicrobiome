#function from: https://joey711.github.io/phyloseq-extensions/edgeR.html

phyloseq_to_edgeR = function(physeq, method="RLE", ...){
    require("edgeR")
    require("phyloseq")
# Enforce orientation.
    if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
    x = as(otu_table(physeq), "matrix")
# Add one to protect against overflow, log(0) issues.
    x = x + 1
# Define gene annotations (`genes`) as tax_table
    taxonomy = tax_table(physeq, errorIfNULL=FALSE)
    if( !is.null(taxonomy) ){
        taxonomy = data.frame(as(taxonomy, "matrix"))
    }
# Now turn into a DGEList
    y = DGEList(counts=x, genes=taxonomy, remove.zeros = TRUE, ...)
# Calculate the normalization factors
    z = calcNormFactors(y, method=method)
# Check for division by zero inside `calcNormFactors`
    if( !all(is.finite(z$samples$norm.factors)) ){
            stop("Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument")
    }
# Estimate dispersions
    return(estimateTagwiseDisp(estimateCommonDisp(z)))
    }
