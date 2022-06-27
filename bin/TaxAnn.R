#!/usr/bin/Rscript

##Sequence cleaning and Multiamplicon pipeline for Access Array of Eimeria infection experiment

library("phyloseq")
library("data.table")
library("taxonomizr")
library(dada2)
library(DECIPHER)
#library("MultiAmplicon", lib.loc="/usr/local/lib/R/site-library") 
## using the devel
devtools::load_all("/SAN/Susanas_den/MultiAmplicon/")

# load MA object from single Amplicon dataset
MA <- readRDS("/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/MATax18S.Rds")

################################################################################################
################ New taxonomic annotation

# first our databases need the ENA taxonomy

dadasilva16 <- readDNAStringSet("/SAN/Susanas_den/AmpMarkers/silva_nr99_v138.1_wSpecies_train_set.fa.gz")

silva16 <- DNAStringSet(readRNAStringSet("/SAN/Susanas_den/AmpMarkers/SILVAdb/slv_ssu138.1/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"))

tax <- read.csv("/SAN/Susanas_den/AmpMarkers/SILVAdb/slv_ssu138.1/taxmap_embl-ebi_ena_ssu_ref_nr99_138.1.txt", sep="\t", stringsAsFactors=FALSE, head=T)

length(rownames(tax))==length(names(silva16))
tax$primaryAccession==accessions

slv <- names(silva16)

accessions <- gsub("([^A-Za-z0-9]*)(\\.)(.*)", "\\1", slv)
#accessions <- gsub("([^A-Za-z0-9]*)(\\.\\d*\\.\\d*)(.*)", "\\1\\2", slv)

# I also want species in here
tax$path <- paste(tax$submitted_path, tax$submitted_name, ";", sep="")

#now I want to replace the names of silva16 with tax$path based on accessions
accessions <- sort(accessions)
tax$primaryAccession <- sort(tax$primaryAccession)

#sanity check
summary((tax$primaryAccession)==(accessions))

names(silva16) <- tax$path

head(names(silva16))

saveRDS(silva16, "/SAN/Susanas_den/AmpMarkers/SILVAdb/silva138_ENA.rds")

# first we need to know which genes are we targeting
primerL <- read.csv("/SAN/Susanas_den/HMHZ/data/primerInputUnique.csv")
#quick fix here
primerL$Primer_name[122] <- "27M_F_98_F.Klin0341_CR_18_R"
target <- primerL[primerL$Primer_name%in%names(MA@PrimerPairsSet),]

##### ok now we make our sequences into DECIPHER format
seqs <- getSequencesFromTable(MA)
#quick sanity check
seqs2 <- getSequenceTableNoChime(MA)
seqs2 <- lapply(seqs2, colnames)
seqs[[3]]==seqs2[[3]]
#reads into DNAstring format
seqs <- lapply(seqs, DNAStringSet)

