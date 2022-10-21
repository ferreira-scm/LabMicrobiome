library(Biostrings)
library(ShortRead)

#### geting qiime2 silva DB and fixing names for dada2
Slv138.fa <- readDNAStringSet("/SAN/Susanas_den/AmpMarkers/RESCRIPt/SSURef_NR99/Fastas/silva-138.1-ssu-nr99-seqs.fasta", format="fasta")
Slv138.tax <- read.table("/SAN/Susanas_den/AmpMarkers/RESCRIPt/SSURef_NR99/Fastas/silva-138.1-ssu-nr99-tax.tsv", sep='\t', header=TRUE)

nrow(Slv138.tax)
length(Slv138.fa)

# sanity checks
which(!names(Slv138.fa)%in%Slv138.tax$Feature.ID)
which(names(Slv138.fa)==Slv138.tax$Feature.ID)
which(!names(Slv138.fa)==Slv138.tax$Feature.ID)
which(!Slv138.tax$Feature.ID%in%names(Slv138.fa))
# ups, gotta reoder
reorder_idx <- match(names(Slv138.fa), Slv138.tax$Feature.ID)
Slv138.tax <- Slv138.tax[reorder_idx,]
which(!names(Slv138.fa)==Slv138.tax$Feature.ID)

# this is a quick fix until I download the database again without domain
Slv138.tax$Taxon <- (gsub("^d__[A-Za-z]+;", "", Slv138.tax$Taxon))

# renaming
names(Slv138.fa) <- gsub(" ", "", Slv138.tax$Taxon, fixed=TRUE)
names(Slv138.fa)[1:100]

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

# this is a quick fix until I download the database again without domain
Slv138LSU.tax$Taxon <- (gsub("^d__[A-Za-z]+;", "", Slv138LSU.tax$Taxon))

# renaming
names(Slv138LSU.fa) <- gsub(" ", "", Slv138LSU.tax$Taxon, fixed=TRUE)
names(Slv138LSU.fa)[1:100]

# now saving
writeFasta(Slv138LSU.fa, "/SAN/Susanas_den/AmpMarkers/RESCRIPt/LSURef_NR99/Fastas/Slv138LSU.dada2.fa")

#### NCBI - others
#### geting qiime2 silva DB and fixing names for dada2
other.fa <- readDNAStringSet("/SAN/Susanas_den/AmpMarkers/RESCRIPt/other/Fastas/NCBI_seqs.fasta", format="fasta")

#other.fa <- readDNAStringSet("/SAN/Susanas_den/AmpMarkers/RESCRIPt/other/tmp/Fastas/NCBI_seqs.fasta", format="fasta")
#other.tax <- read.table("/SAN/Susanas_den/AmpMarkers/RESCRIPt/other/tmp/Fastas/NCBI_tax.tsv", sep="\t", header=TRUE, quote="")

other.tax <- read.table("/SAN/Susanas_den/AmpMarkers/RESCRIPt/other/Fastas/NCBI_tax.tsv", sep="\t", header=TRUE, quote="")

nrow(other.tax)
length(other.fa)

# sanity checks
which(!names(other.fa)%in%other.tax$Feature.ID)

which(names(other.fa)==other.tax$Feature.ID)

which(!names(other.fa)==other.tax$Feature.ID)

which(!other.tax$Feature.ID%in%names(other.fa))

# ups, gotta reoder
reorder_idx <- match(names(other.fa), other.tax$Feature.ID)
other.tax <- other.tax[reorder_idx,]
which(!names(other.fa)==other.tax$Feature.ID)

# this is a quick fix until I download the database again without domain
#other.tax$Taxon <-gsub(" ", "", other.tax$Taxon, fixed=TRUE)
#other.tax$Taxon <- (gsub("^d__[A-Za-z]+;", "", other.tax$Taxon))

# renaming
names(other.fa) <- gsub(" ", "", other.tax$Taxon, fixed=TRUE)

names(other.fa)[1:100]

# now saving
#writeFasta(other.fa, "/SAN/Susanas_den/AmpMarkers/RESCRIPt/other/tmp/Fastas/other.dada2.fa")

writeFasta(other.fa, "/SAN/Susanas_den/AmpMarkers/RESCRIPt/other/Fastas/other.dada2.fa")

