library(DECIPHER)

library(taxonomizr)

## for ITS
seqs <- readDNAStringSet("/SAN/Susanas_den/AmpMarkers/UNITE_ITS.fasta")

seqs <- RemoveGaps(seqs)

taxid <- NULL

                                        # ensure that all sequences are in the same orientation:
seqs <- OrientNucleotides(seqs)

# specific name for the classifier
names(seqs) <- gsub("(.*)(k__)", "\\1 Root; k__", names(seqs))

groups <- names(seqs)

groups <- gsub("(.*)(Root;)", "\\2", groups)

groups

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

print(sum(remove))

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

saveRDS(trainingSet, "/SAN/Susanas_den/AmpMarkers/ITSTrainingSet.RDS")

head(names(seqs))

