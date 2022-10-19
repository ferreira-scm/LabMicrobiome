library(Biostrings)
library(ShortRead)

#### geting qiime2 silva DB and fixing names for dada2
Slv138.fa <- readDNAStringSet("/SAN/Susanas_den/AmpMarkers/RESCRIPt/SSURef_NR99/Fastas/silva-138.1-ssu-nr99-seqs.fasta", format="fasta")
Slv138.tax <- read.table("/SAN/Susanas_den/AmpMarkers/RESCRIPt/SSURef_NR99/Fastas/silva-138.1-ssu-nr99-tax.tsv", sep="\t", header=TRUE)

#sanity checks

head(Slv138.fa)

Slv138.tax$Feature.ID

nrow(Slv138.tax)

Slv138.tax$Feature.ID[440464]

names(Slv138.fa)[440464]

grep("CP013077", Slv138.tax$Feature.ID)

Slv138.tax$Feature.ID[70924]

names(Slv138.fa)[1]

# sanity checks
which(!names(Slv138.fa)%in%Slv138.tax$Feature.ID)
which(names(Slv138.fa)==Slv138.tax$Feature.ID)
which(!names(Slv138.fa)==Slv138.tax$Feature.ID)
which(!Slv138.tax$Feature.ID%in%names(Slv138.fa))
# ups, gotta reoder
reorder_idx <- match(names(Slv138.fa), Slv138.tax$Feature.ID)
Slv138.tax <- Slv138.tax[reorder_idx,]
which(!names(Slv138.fa)==Slv138.tax$Feature.ID)

# renaming
names(Slv138.fa) <- gsub(" ", "", Slv138.tax$Taxon, fixed=TRUE)
names(Slv138.fa)

# now saving
writeFasta(Slv138.fa, "/SAN/Susanas_den/AmpMarkers/RESCRIPt/SSURef_NR99/Fastas/Slv138.dada2.fa")

#### LSU
#### geting qiime2 silva DB and fixing names for dada2
Slv138LSU.fa <- readDNAStringSet("/SAN/Susanas_den/AmpMarkers/RESCRIPt/LSURef_NR99/Fastas/silva-138.1-lsu-nr99-seqs.fasta", format="fasta")
Slv138LSU.tax <- read.table("/SAN/Susanas_den/AmpMarkers/RESCRIPt/LSURef_NR99/Fastas/silva-138.1-lsu-nr99-tax.tsv", sep="\t", header=TRUE)

# sanity checks
which(!names(Slv138LSU.fa)%in%Slv138LSU.tax$Feature.ID)
which(names(Slv138LSU.fa)==Slv138LSU.tax$Feature.ID)
which(!names(Slv138LSU.fa)==Slv138LSU.tax$Feature.ID)
which(!Slv138LSU.tax$Feature.ID%in%names(Slv138LSU.fa))

# ups, gotta reoder
reorder_idx <- match(names(Slv138LSU.fa), Slv138LSU.tax$Feature.ID)
Slv138LSU.tax <- Slv138LSU.tax[reorder_idx,]
which(!names(Slv138LSU.fa)==Slv138LSU.tax$Feature.ID)

# renaming
names(Slv138LSU.fa) <- gsub(" ", "", Slv138LSU.tax$Taxon, fixed=TRUE)
names(Slv138LSU.fa)[1:100]

# now saving
writeFasta(Slv138LSU.fa, "/SAN/Susanas_den/AmpMarkers/RESCRIPt/LSURef_NR99/Fastas/Slv138LSU.dada2.fa")

#### NCBI - others
#### geting qiime2 silva DB and fixing names for dada2
other.fa <- readDNAStringSet("/SAN/Susanas_den/AmpMarkers/RESCRIPt/other/Fastas/NCBI_seqs.fasta", format="fasta")
other.tax <- read.table("/SAN/Susanas_den/AmpMarkers/RESCRIPt/other/Fastas/NCBI_tax.tsv", sep="\t", header=TRUE)

head(other.tax)


other.tax$Taxon

# sanity checks
which(!names(Slv138LSU.fa)%in%Slv138LSU.tax$Feature.ID)
which(names(Slv138LSU.fa)==Slv138LSU.tax$Feature.ID)
which(!names(Slv138LSU.fa)==Slv138LSU.tax$Feature.ID)
which(!Slv138LSU.tax$Feature.ID%in%names(Slv138LSU.fa))

# ups, gotta reoder
reorder_idx <- match(names(Slv138LSU.fa), Slv138LSU.tax$Feature.ID)
Slv138LSU.tax <- Slv138LSU.tax[reorder_idx,]
which(!names(Slv138LSU.fa)==Slv138LSU.tax$Feature.ID)

# renaming
names(Slv138LSU.fa) <- gsub(" ", "", Slv138LSU.tax$Taxon, fixed=TRUE)
names(Slv138LSU.fa)[1:100]

# now saving
writeFasta(Slv138LSU.fa, "/SAN/Susanas_den/AmpMarkers/RESCRIPt/LSURef_NR99/Fastas/Slv138LSU.dada2.fa")

