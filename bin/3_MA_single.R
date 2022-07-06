#!/usr/bin/Rscript

##Sequence cleaning and Multiamplicon pipeline for Access Array of Eimeria infection experiment

library("lifecycle", lib.loc="/usr/local/lib/R/site-library") 
library("ggplot2")
library("reshape")
library("phyloseq")
library("data.table")
library("taxonomizr")
library("taxize")
library("parallel")
library(dada2)
library(DECIPHER)
#library("MultiAmplicon", lib.loc="/usr/local/lib/R/site-library") 

## using the devel
devtools::load_all("/SAN/Susanas_den/MultiAmplicon/")


## re-run or use pre-computed results for different parts of the pipeline:
## Set to FALSE to use pre-computed and saved results, TRUE to redo analyses.
doFilter <- FALSE
doMultiAmp <- FALSE
doTax <- FALSE

###################Full run Microbiome#######################
#Preparation of files

##These are the same steps that are followed by the DADA2 pipeline
path <- "/SAN/Victors_playground/Eimeria_microbiome/2018_22_Mmb_1/" ## change according to where you downloaded
fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) #take all fastaq files from the folder 
fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE) #separate the forward reads
fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE) #separate the reverse reads 
samples <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", basename(fastqF))
samples<- gsub("S\\d+-", "\\1", basename(samples))

#Extra step in the pipeline: quality plots of the reads 
#plotQualityProfile(fastqF[[1]])
#plotQualityProfile(fastqF[[20]])
#plotQualityProfile(fastqR[[1]])
#plotQualityProfile(fastqR[[20]])

#Creation of a folder for filtrated reads 

filt_path <- "/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/Run18S/filtered18S"

#Pipeline filtration 
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs) <- samples
filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs) <- samples
## some files will be filtered out completely, therefore allowing 50
## files less present and still don't redo filtering

if(doFilter){
  lapply(seq_along(fastqF),  function (i) {
    dada2::filterAndTrim(fastqF[i], filtFs[i], fastqR[i], filtRs[i],
                  truncLen=c(150,170), 
                  maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
                  compress=TRUE, verbose=TRUE, multithread=TRUE)
  })
}

names(filtFs) <- names(filtRs) <- samples
files <- PairedReadFileSet(filtFs, filtRs)

#Preparation of primer file ### Here stats the Multiamplicon pipeline from Emanuel
#Primers used in the arrays 
ptable <- read.csv(file = "/SAN/Victors_playground/Eimeria_microbiome/primer.file.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
primerF <- ptable[, "TS.SequenceF"]
primerR <- ptable[, "TS.SequenceR"]
names(primerF) <- as.character(ptable[, "corrected.NameF"])
names(primerR) <- as.character(ptable[, "corrected.NameR"])
primer <- PrimerPairsSet(primerF, primerR)

##Multi amplicon pipeline
if(doMultiAmp){
  MA <- MultiAmplicon(primer, files)
  filedir <- "/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/Run18S/Stratified_files18S"
  if(dir.exists(filedir)) unlink(filedir, recursive=TRUE)
  MA <- sortAmplicons(MA, n=1e+05, filedir=filedir) ## This step sort the reads into amplicons based on the number of primer pairs
  errF <-  dada2::learnErrors(unlist(getStratifiedFilesF(MA)), nbase=1e8,
                       verbose=0, multithread = 90)
  errR <- dada2::learnErrors(unlist(getStratifiedFilesR(MA)), nbase=1e8,
                      verbose=0, multithread = 90)
  MA <- dadaMulti(MA, Ferr=errF, Rerr=errR,  pool=FALSE,
                  verbose=0, mc.cores=90)
  MA <- mergeMulti(MA, mc.cores=90)
  propMerged <- MultiAmplicon::calcPropMerged(MA)
  MA <- mergeMulti(MA, mc.cores=90) 
  MA <- makeSequenceTableMulti(MA, mc.cores=90) 
  MA <- removeChimeraMulti(MA, mc.cores=90)
  saveRDS(MA, "/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/Run18S/MA_single_18S.RDS")
} else{
    MA <- readRDS("/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/Run18S/MA_single_18S.RDS")
}

#trackingF <- getPipelineSummary(MA) 
## plotPipelineSummary(trackingF) 
## plotPipelineSummary(trackingF) + scale_y_log10()

###New taxonomic assignment 
if(doTax){
MA <- blastTaxAnnot(MA,
                    db = "/SAN/db/blastdb/nt/nt",
                    negative_gilist = "/SAN/db/blastdb/uncultured.gi",
                    infasta = "/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/Run18S/18Sin.fasta",
                    outblast = "/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/Run18S/18Sout.fasta",
                    taxonSQL = "/SAN/db/taxonomy/taxonomizr.sql", 
                    num_threads = 90)
saveRDS(MA, file="/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/MATax18S.Rds")
}else{
    MA <- readRDS("/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/MATax18S.Rds")
}

### Add sample information
if(!exists("sample.data")){
     source("bin/1_Data_preparation.R")
 }

if(!exists("sdt")){
    source("bin/1_qPCR_data_preparation.R")
}

#little fix
rownames(sdt) <- sdt$labels
rownames(sdt)==rownames(sample.data)

##To phyloseq
##Sample data
# temporary fix
source("bin/toPhyloseq.R")
PS <- TMPtoPhyloseq(MA, colnames(MA))

PS@sam_data <- sample_data(sdt)
saveRDS(PS, file="/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/PhyloSeqData18S.Rds")

##Primer data
PS.l <- TMPtoPhyloseq(MA, colnames(MA),  multi2Single=FALSE)

# adding sample data
for (i in 1:3)
{
    sample_data(PS.l[[i]]) <- sdt
}
saveRDS(PS.l, file="/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/PhyloSeqList18S.Rds")

###For Microbiome analysis (Victor)
PS.18S <- PS.l[[2]]
#PS18S <- phyloseq(
#    otu_table(PS.l$wang1141_13_F.Nem_0425_6_3_R), 
#    sample_data(PS.l$wang1141_13_F.Nem_0425_6_3_R), 
#    tax_table(PS.l$wang1141_13_F.Nem_0425_6_3_R))
#sum(otu_table(PS.18S)) ##Total denoised reads = 853,134
saveRDS(PS.18S, file="/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/PS_18Swang.Rds") ###Information from 18S


################################################################################################
################ New taxonomic annotation

# first we need to know which genes are we targeting
#primerL <- read.csv("/SAN/Susanas_den/HMHZ/data/primerInputUnique.csv")
#quick fix here
#primerL$Primer_name[122] <- "27M_F_98_F.Klin0341_CR_18_R"
#target <- primerL[primerL$Primer_name%in%names(MA@PrimerPairsSet),]

##### ok now we make our sequences into DECIPHER format
#seqs <- getSequencesFromTable(MA)
#quick sanity check
#seqs2 <- getSequenceTableNoChime(MA)
#seqs2 <- lapply(seqs2, colnames)
#seqs[[3]]==seqs2[[3]]
#reads into DNAstring format
#seqs <- lapply(seqs, DNAStringSet)

#Load our training sets
#trainingSet16S <- readRDS("/SAN/Susanas_den/AmpMarkers/16SSilva138TrainingSet.RDS")
#trainingSet18S <- readRDS("/SAN/Susanas_den/AmpMarkers/18SSilva132TrainingSet.RDS")

#little fix
#trainingSet16S$ranks <-  c("root", "superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain")

#trainingSet18S$ranks <-  c("root", "superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain")

#ranks <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain") # ranks of interest


#ranks2 <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")

#newTax <- list()
#idtaxa <- list()
#for (i in 1:3){
#    if (target$Gen[i] == "16S") {
#idtaxa <- IdTaxa(seqs[[i]],
#                 trainingSet16S,
#                 strand = "both",
#                 threshold = 40,
#                 bootstraps = 100,
#                 processors = NULL,
#                 verbose = TRUE,
#                 type = "extended")
#newTax[[i]] <- t(sapply(idtaxa, function(x) {
#    m <- match(ranks, x$rank)
#    taxa <- x$taxon[m]
#    taxa[startsWith(taxa, "unclassified_")] <- NA
#    taxa
#}))
#        colnames(newTax[[i]]) <- ranks
#        rownames(newTax[[i]]) <- getSequencesFromTable(MA)[[i]]
#    } else if (target$Gen[i]=="18S"){
#idtaxa <- IdTaxa(seqs[[i]],
#                 trainingSet18S,
#                 strand = "both",
#                 threshold = 40,
#                 bootstraps = 100,
#                 processors = NULL,
#                 verbose = TRUE,
#                 type = "extended")
#newTax[[i]] <- t(sapply(idtaxa, function(x) {
#    m <- match(ranks, x$rank)
#    taxa <- x$taxon[m]
#    taxa[startsWith(taxa, "unclassified_")] <- NA
#    taxa
#}))
#        colnames(newTax[[i]]) <- ranks
#        rownames(newTax[[i]]) <- getSequencesFromTable(MA)[[i]]
#        }
#}

# sanity check
#rownames(MA@taxonTable[[3]])==rownames(newTax[[3]])

#MA1 <- MA

#for (i in (1:3)) {
#    MA1@taxonTable[[i]] <- newTax[[i]]
#    }

# not working
#PS1 <- TMPtoPhyloseq(MA1, colnames(MA1))

#get_taxa_unique(PS1, "species")
# sanity check
#a <- lapply(newTax, rownames)
#a <- unlist(a)
#a==rownames(PS@tax_table)
###################### # we try with dada2 (naive bayesian classifier)
######################################################################

seqs <- getSequencesFromTable(MA)
seqs <- lapply(seqs, DNAStringSet)

primerL <- read.csv("/SAN/Susanas_den/HMHZ/data/primerInputUnique.csv")
#quick fix here
primerL$Primer_name[122] <- "27M_F_98_F.Klin0341_CR_18_R"
target <- primerL[primerL$Primer_name%in%names(MA@PrimerPairsSet),]

target$Gen

taxa2 <- list()

for (i in 1:3){
    if (target$Gen[i]=="16S"){
        taxa2[[i]] <- assignTaxonomy(seqs[[i]],
                                     "/SAN/Susanas_den/AmpMarkers/silva_nr99_v138.1_wSpecies_train_set.fa.gz",
                                    multithread=90,
                                    tryRC = TRUE,
                                    verbose=TRUE)
    }
    else if (target$Gen[i]=="18S"){
     taxa2[[i]] <- assignTaxonomy(seqs[[i]],
                                  "/SAN/Susanas_den/AmpMarkers/silva132.18Sdada2_mod.fa.gz",
                                    multithread=90,
                                    tryRC = TRUE,
                                    verbose=TRUE)
    }   
}

## Little inspection
taxa.print <- taxa2[[2]]
rownames(taxa.print) <- NULL

taxa.print[222:240,]

## assign species to 16S

spec <- list()

for (i in 1:3){
    if (target$Gen[i]=="16S"){
        spec[[i]] <- assignSpecies(seqs[[i]],
                                   "/SAN/Susanas_den/AmpMarkers/silva_species_assignment_v138.1.fa.gz",
                                    allowMultiple=TRUE,
                                    tryRC = TRUE,
                                    verbose=TRUE)
    }}

#for (i in 1:3){
#    if (target$Gen[i]=="18S"){
#        spec[[i]] <- assignSpecies(seqs[[i]],
#                                   "/SAN/Susanas_den/AmpMarkers/silva132.18Sdada2_AssignSpecies.fa",
#                                    allowMultiple=TRUE,
#                                    tryRC = TRUE,
#                                    verbose=TRUE)
#    }}


s.print <- spec[[2]]

rownames(s.print) <- NULL

s.print

##replacing species name
rownames(tax_table(PS.l[[1]]))==rownames(spec[[1]])
rownames(spec[[1]])==rownames(taxa2[[1]])

taxa2[[1]][,7] <- spec[[1]][,2]
taxa2[[3]][,7] <- spec[[3]][,2]

## now saving this taxonomy

MA2 <- MA

MA2@taxonTable <- taxa2


##To phyloseq
##Sample data
# temporary fix
source("bin/toPhyloseq.R")
PS <- TMPtoPhyloseq(MA2, colnames(MA2))

PS@sam_data <- sample_data(sdt)
saveRDS(PS, file="/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/PhyloSeqData18S_SILVA.Rds")

##Primer data
PS.l <- TMPtoPhyloseq(MA2, colnames(MA2),  multi2Single=FALSE)

# adding sample data
for (i in 1:3)
{
    sample_data(PS.l[[i]]) <- sdt
}
saveRDS(PS.l, file="/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/PhyloSeqList18S_SILVA.Rds")

###For Microbiome analysis (Victor)
PS.18S <- PS.l[[2]]
#PS18S <- phyloseq(
#    otu_table(PS.l$wang1141_13_F.Nem_0425_6_3_R), 
#    sample_data(PS.l$wang1141_13_F.Nem_0425_6_3_R), 
#    tax_table(PS.l$wang1141_13_F.Nem_0425_6_3_R))
#sum(otu_table(PS.18S)) ##Total denoised reads = 853,134
saveRDS(PS.18S, file="/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/PS_18Swang_SILVA.Rds") ###Information from 18S



##########################################
########################################
#ID <- taxonomizr::accessionToTaxa(accession, version="base", "/SAN/db/taxonomy/taxonomizr.sql")
#taxonomy <- taxonomizr::getTaxonomy(ID,
#                                    "/SAN/db/taxonomy/taxonomizr.sql",
#                                    desiredTaxa = c("superkingdom", "phylum", "class", "order", "family", "genus", "species"))


