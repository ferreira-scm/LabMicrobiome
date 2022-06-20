library(DECIPHER)

library(taxonomizr)




################ New taxonomic annotation
# first we need to know which genes are we targeting
primerL <- read.csv("/SAN/Susanas_den/HMHZ/data/primerInputUnique.csv")
#quick fix here
primerL$Primer_name[122] <- "27M_F_98_F.Klin0341_CR_18_R"

target <- primerL[primerL$Primer_name%in%names(MA@PrimerPairsSet),]

target$Gen

# ok now we make our sequences into DECIPHER format
#seqs <- getSequencesFromTable(MA)
#seqs <- lapply(seq_along(seqs), function(x){
#    DNAStringSet(seqs[[x]])
#    })

#Make our training sets
seqs <- readDNAStringSet("/SAN/Susanas_den/AmpMarkers/silva_nr99_v138.1_wSpecies_train_set.fa.gz")

#trainingSet18S <- readDNAStringSet("/SAN/Susanas_den/AmpMarkers/silva132.18Sdada2.fa.gz")

#trainingSetITS <- readDNAStringSet("/SAN/Susanas_den/AmpMarkers/UNITE_ITS.fasta")

fa16 <- read.table("/SAN/Susanas_den/AmpMarkers/silva_nr99_v138.1_wSpecies_train_set.fa.gz")


# attepting to create my own taxID
ranks <- strsplit(names(seqs), ";", fix=T)
count <- 1L
groups <- "Root"
index <- -1L
level <- 0L
rank <- "rootrank"
pBar <- txtProgressBar(style=3)
taxa <- setNames(c("domain", "phylum", "order", "family", "genus", "species"),
                 c("d__", "p__", "o__", "f__", "g__","s__"))

for (i in seq_along(ranks)) {
      for (j in seq_along(ranks[[i]])) {
          rank_level <- taxa[substring(ranks[[i]][j], 1, 3)]
          group <- substring(ranks[[i]][j], 4)
          w <- which(groups==group & rank==rank_level)
          if (length(w) > 0) {
              parent <- match(substring(ranks[[i]][j - 1], 4),
                              groups)
              if (j==1 || any((parent - 1L)==index[w]))
                  next # already included
          }
          count <- count + 1L
          groups <- c(groups, group)
          if (j==1) {
              index <- c(index, 0)
          } else {
              parent <- match(substring(ranks[[i]][j - 1], 4),
                              groups)
              index <- c(index,
                         parent - 1L)
          }
          level <- c(level, j)
          rank <- c(rank, taxa[j])
      }
      setTxtProgressBar(pBar, i/length(ranks))
      }

level

rank

groups <- gsub("^[ ]+", "", groups)

groups <- gsub("[ ]+$", "", groups)

taxid <- paste(0:(length(index) - 1L), groups, index, level, rank, sep="*")

head(taxid, n=10)

writeLines(taxid,
           con="/SAN/Susanas_den/AmpMarkers/TAXID_Silva138.txt")

saveRDS(taxid, "/SAN/Susanas_den/AmpMarkers/TAXID_Silva138.RDS")
