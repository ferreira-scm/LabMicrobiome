library(DECIPHER)
library(taxonomizr)

#Make training sets for ID taxa
seqs <- readDNAStringSet("/SAN/Susanas_den/AmpMarkers/silva_nr99_v138.1_wSpecies_train_set.fa.gz", format="fasta")

# specific name for the classifier
names(seqs)[500:510]

names(seqs) <- paste("Root;", names(seqs), sep="")

TrainingDB <- function(seqs){
# if they exist, remove any gaps in the sequences:
seqs <- RemoveGaps(seqs)
taxid <- NULL
# ensure that all sequences are in the same orientation:
seqs <- OrientNucleotides(seqs)
groups <- names(seqs)
groupCounts <- table(groups)
u_groups <- names(groupCounts)
length(u_groups)
taxid <- NULL
#subset sequences in large groups
maxGroupSize <- 10
remove <- logical(length(seqs))
for (i in which(groupCounts>maxGroupSize)) {
    index <- which(groups==u_groups[i])
    keep <- sample(length(index),
                   maxGroupSize)
    remove[index[-keep]] <- TRUE
}
sum(remove) 
## the actual training of the database
maxIterations <- 3
allowGroupRemoval <- FALSE
probSeqsPrev <- integer()
for (i in seq_len(maxIterations)) {
    cat("Training iteration: ", i, "\n", sep="")
# train the classifier
    trainingSet <- LearnTaxa(seqs[!remove],
                             names(seqs)[!remove],
                             taxid)
# look for problem sequences
    probSeqs <- trainingSet$problemSequences$Index
    if (length(probSeqs)==0) {
        cat("No problem sequences remaining.\n")
        break
    } else if (length(probSeqs)==length(probSeqsPrev) &&
               all(probSeqsPrev==probSeqs)) {
        cat("Iterations converged.\n")
        break
    }
    if (i==maxIterations)
        break
    probSeqsPrev <- probSeqs
# remove any problem sequences
    index <- which(!remove)[probSeqs]
    remove[index] <- TRUE # remove all problem sequences
    if (!allowGroupRemoval) {
# replace any removed groups
        missing <- !(u_groups %in% groups[!remove])
        missing <- u_groups[missing]
        if (length(missing) > 0) {
            index <- index[groups[index] %in% missing]
            remove[index] <- FALSE # don't remove
        }
    }
}
}
sum(remove) # total number of sequences eliminated
length(probSeqs) # number of remaining problem sequences

trainingSet <- TrainingDB(seqs)
saveRDS(trainingSet, "/SAN/Susanas_den/AmpMarkers/16SSilva138TrainingSet.RDS")

###################### for 18S
seqs <- readDNAStringSet("/SAN/Susanas_den/AmpMarkers/silva132.18Sdada2.fa.gz")
head(names(seqs))
# specific name for the classifier
names(seqs)<- gsub("([A-Z0-9]*\\.[0-9]*\\.[0-9]*;$)", "", names(seqs))
names(seqs) <- paste("Root;", names(seqs), sep="")
trainingSet18S <- TrainingDB(seqs)
saveRDS(trainingSet18S, "/SAN/Susanas_den/AmpMarkers/18SSilva132TrainingSet.RDS")


## for ITS
seqs <- readDNAStringSet("/SAN/Susanas_den/AmpMarkers/UNITE_ITS.fasta")

head(seqs)

# specific name for the classifier
names(seqs) <- gsub("(.*)(k__)", "\\1 Root;k__", names(seqs))
names(seqs) <- gsub("(.*)(Root;)", "\\2", names(seqs))
names(seqs) <- gsub("(\\w__)","", names(seqs))
names(seqs) <- paste(names(seqs), ";", sep="")

# if they exist, remove any gaps in the sequences:
seqs <- RemoveGaps(seqs)
taxid <- NULL

# ensure that all sequences are in the same orientation:
seqs <- OrientNucleotides(seqs)

groups <- names(seqs)
groupCounts <- table(groups)
u_groups <- names(groupCounts)
length(u_groups)

#subset sequences in large groups
maxGroupSize <- 10
remove <- logical(length(seqs))

for (i in which(groupCounts>maxGroupSize)) {
    index <- which(groups==u_groups[i])
    keep <- sample(length(index),
                   maxGroupSize)
    remove[index[-keep]] <- TRUE
}
sum(remove) 

## the actual training of the database
maxIterations <- 3
allowGroupRemoval <- FALSE
probSeqsPrev <- integer()

for (i in seq_len(maxIterations)) {
    cat("Training iteration: ", i, "\n", sep="")
# train the classifier
    trainingSet <- LearnTaxa(seqs[!remove],
                             names(seqs)[!remove],
                             taxid)
# look for problem sequences
    probSeqs <- trainingSet$problemSequences$Index
    if (length(probSeqs)==0) {
        cat("No problem sequences remaining.\n")
        break
    } else if (length(probSeqs)==length(probSeqsPrev) &&
               all(probSeqsPrev==probSeqs)) {
        cat("Iterations converged.\n")
        break
    }
    if (i==maxIterations)
        break
    probSeqsPrev <- probSeqs
# remove any problem sequences
    index <- which(!remove)[probSeqs]
    remove[index] <- TRUE # remove all problem sequences
    if (!allowGroupRemoval) {
# replace any removed groups
        missing <- !(u_groups %in% groups[!remove])
        missing <- u_groups[missing]
        if (length(missing) > 0) {
            index <- index[groups[index] %in% missing]
            remove[index] <- FALSE # don't remove
        }
    }
}


saveRDS(trainingSet,"/SAN/Susanas_den/AmpMarkers/ITS_UNITETrainingSet.RDS")


ITS <- readRDS("/SAN/Susanas_den/AmpMarkers/ITS_UNITETrainingSet.RDS")

### For LSU
seqs <- DNAStringSet(readRNAStringSet("/SAN/Susanas_den/AmpMarkers/SILVAdb/slv_lsu138.1/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz"))

head(names(seqs))

#enaTax <- read.csv("/SAN/Susanas_den/AmpMarkers/SILVAdb/slv_lsu138.1/taxmap_embl-ebi_ena_ssu_ref_nr99_138.1.txt", sep="\t", header=T)

ACC <- read.csv("/SAN/Susanas_den/AmpMarkers/SILVAdb/slv_lsu138.1/tax_slv_lsu_138.1.acc_taxid.gz", sep="\t")

ENA <- read.csv("/SAN/Susanas_den/AmpMarkers/SILVAdb/slv_lsu138.1/taxmap_embl-ebi_ena_ssu_ref_nr99_138.1.txt", sep="\t", header=T)

head(ENA)

head(names(seqs))

length(unique(ENA$ncbi_taxonid))

length(ENA$ncbi_taxonid)

tax <- enaTax$submitted_path

(tax)[100:110]

names(seqs) <- (gsub("(.*?)( .*)", "Root;\\2", names(seqs)))
names(seqs) <- gsub("(Root; )", "Root;", names(seqs)) 
