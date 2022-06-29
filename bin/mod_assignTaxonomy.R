

#this is a test


MA <- readRDS("/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/MATax18S.Rds")
seqs <- getSequencesFromTable(MA)
seqs <- seqs[[1]]

seqs

refFasta <- "/SAN/Susanas_den/AmpMarkers/SILVAdb/silva138_ENA_nosp.fa"



seqs
refFasta
minBoot=50
tryRC=FALSE
outputBootstraps=FALSE
taxLevels=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
multithread=90
verbose=FALSE


    MIN_REF_LEN <- 20 # Enforced minimum length of reference seqs. Must be bigger than the kmer-size used (8).
    MIN_TAX_LEN <- 50 # Minimum length of input sequences to get a taxonomic assignment
                                        # Get character vector of sequences

seqs <- getSequences(seqs)

if(min(nchar(seqs)) < MIN_TAX_LEN) {
        warning("Some sequences were shorter than ", MIN_TAX_LEN, " nts and will not receive a taxonomic classification.")
    }

# Read in the reference fasta
refsr <- readFasta(refFasta)

    lens <- width(sread(refsr))
    if(any(lens<MIN_REF_LEN)) {
        refsr <- refsr[lens>=MIN_REF_LEN]
        warning(paste0("Some reference sequences were too short (<", MIN_REF_LEN, "nts) and were excluded."))
    }

refs <- as.character(sread(refsr))

tax <- as.character(id(refsr))

tax <- sapply(tax, function(x) gsub("^\\s+|\\s+$", "", x)) # Remove leading/trailing whitespace



# Sniff and parse UNITE fasta format
    UNITE <- FALSE
    if(all(grepl("FU\\|re[pf]s", tax[1:10]))) {
        UNITE <- TRUE
        cat("UNITE fungal taxonomic reference detected.\n")
        tax <- sapply(strsplit(tax, "\\|"), `[`, 5)
        tax <- gsub("[pcofg]__unidentified;", "_DADA2_UNSPECIFIED;", tax)
        tax <- gsub(";s__(\\w+)_", ";s__", tax)
        tax <- gsub(";s__sp$", ";_DADA2_UNSPECIFIED", tax)
    }
# Crude format check
    if(!grepl(";", tax[[1]])) {
        if(length(unlist(strsplit(tax[[1]], "\\s")))==3) {
            stop("Incorrect reference file format for assignTaxonomy (this looks like a file formatted for assignSpecies).")
        } else {
            stop("Incorrect reference file format for assignTaxonomy.")
        }
    }

# Parse the taxonomies from the id string
    tax.depth <- sapply(strsplit(tax, ";"), length)


summary(as.factor(tax.depth))

length(tax)

td <- max(tax.depth)
    for(i in seq(length(tax))) {
        if(tax.depth[[i]] < td) {
            for(j in seq(td - tax.depth[[i]])) {
                tax[[i]] <- paste0(tax[[i]], "_DADA2_UNSPECIFIED;")
            }
        }
    }

# Create the integer maps from reference to type ("genus") and for each tax level
    genus.unq <- unique(tax)

genus.unq

ref.to.genus <- match(tax, genus.unq)

    tax.mat <- matrix(unlist(strsplit(genus.unq, ";")), ncol=td, byrow=TRUE)
    tax.df <- as.data.frame(tax.mat)
    for(i in seq(ncol(tax.df))) {
        tax.df[,i] <- factor(tax.df[,i])
        tax.df[,i] <- as.integer(tax.df[,i])
    }
    tax.mat.int <- as.matrix(tax.df)
### Assign
# Parse multithreading argument
    if(is.logical(multithread)) {
        if(multithread==TRUE) { RcppParallel::setThreadOptions(numThreads = "auto") }
        else { RcppParallel::setThreadOptions(numThreads = 1) }
    } else if(is.numeric(multithread)) {
        RcppParallel::setThreadOptions(numThreads = multithread)
    } else {
        warning("Invalid multithread parameter. Running as a single thread.")
        RcppParallel::setThreadOptions(numThreads = 1)
    }
# Run C assignemnt code
    assignment <- C_assign_taxonomy2(seqs, rc(seqs), refs, ref.to.genus, tax.mat.int, tryRC, verbose)
# Parse results and return tax consistent with minBoot
    bestHit <- genus.unq[assignment$tax]
    boots <- assignment$boot
    taxes <- strsplit(bestHit, ";")
    taxes <- lapply(seq_along(taxes), function(i) taxes[[i]][boots[i,]>=minBoot])
                                        # Convert to character matrix
    tax.out <- matrix(NA_character_, nrow=length(seqs), ncol=td)
    for(i in seq(length(seqs))) {
        if(length(taxes[[i]]) > 0) {
            tax.out[i,1:length(taxes[[i]])] <- taxes[[i]]
        }
    }
    rownames(tax.out) <- seqs
    colnames(tax.out) <- taxLevels[1:ncol(tax.out)]
    tax.out[tax.out=="_DADA2_UNSPECIFIED"] <- NA_character_
        if(outputBootstraps){
                                        # Convert boots to integer matrix
            boots.out <- matrix(boots, nrow=length(seqs), ncol=td)
            rownames(boots.out) <- seqs
            colnames(boots.out) <- taxLevels[1:ncol(boots.out)]
            list(tax=tax.out, boot=boots.out)
        } else {
            tax.out
        }
}

                                        # Helper function for assignSpecies
mapHits <- function(x, refs, keep, sep="/") {
    hits <- refs[x]
    hits[grepl("Escherichia", hits, fixed=TRUE) | grepl("Shigella", hits, fixed=TRUE)] <- "Escherichia/Shigella"
    if(length(unique(hits))<=keep) {
        rval <- do.call(paste, c(as.list(sort(unique(hits))), sep=sep))
    } else { rval <- NA_character_  }
    if(length(rval)==0) rval <- NA_character_
                            rval
}

                                        # Match curated genus names to binomial genus names
                                        # Handles Clostridium groups and split genera names
matchGenera <- function(gen.tax, gen.binom, split.glyph="/") {
    if(is.na(gen.tax) || is.na(gen.binom)) { return(FALSE) }
    if((gen.tax==gen.binom) ||
       grepl(paste0("^", gen.binom, "[ _", split.glyph, "]"), gen.tax) ||
       grepl(paste0(split.glyph, gen.binom, "$"), gen.tax)) {
        return(TRUE)
    } else {
        return(FALSE)
    }

