library(DECIPHER)

library(taxonomizr)


#Make our training sets 16S

seqs <- readDNAStringSet("/SAN/Susanas_den/AmpMarkers/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz")
#trainingSet18S <- readDNAStringSet("/SAN/Susanas_den/AmpMarkers/silva132.18Sdada2.fa.gz")


#TrainingDB <- function(seqs){
# if they exist, remove any gaps in the sequences:
seqs <- RemoveGaps(seqs)
taxid <- NULL

# specific name for the classifier
head(names(seqs))
names(seqs) <- (gsub("(.*?)( .*)", "Root;\\2", names(seqs)))
names(seqs) <- gsub("(Root; )", "Root;", names(seqs)) 

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
#}


saveRDS(trainingSet, "/SAN/Susanas_den/AmpMarkers/LSUSilva138TrainingSet.RDS")


sum(remove) # total number of sequences eliminated
length(probSeqs) # number of remaining problem sequences
