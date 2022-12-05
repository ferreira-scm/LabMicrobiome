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
library(decontam)
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

#there's 2 samples without metadata
rownames(PS@sam_data)[rownames(PS@sam_data)%in%rownames(sdt)==FALSE]

rownames(PS@sam_data)%in%rownames(sample.data)
rownames(sdt)[rownames(sdt)%in%rownames(PS@sam_data)==FALSE]

PS <- subset_samples(PS, rownames(PS@sam_data)%in%rownames(sdt))
sdt <- sdt[match(rownames(PS@sam_data), sdt$labels),]

#sanity check
rownames(PS@otu_table) == rownames(sdt)

# adding sample data slot
PS@sam_data <- sample_data(sdt)
# sanity check
rownames(PS@sam_data)==rownames(PS@otu_table)
saveRDS(PS, file="/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/PhyloSeqData18S.Rds")

##Primer data
PS.l <- TMPtoPhyloseq(MA, colnames(MA),  multi2Single=FALSE)


for (i in 1:3) {
    PS.l[[i]] <- subset_samples(PS.l[[i]], rownames(PS.l[[i]]@sam_data)%in%rownames(sdt))
}

#sanity check
rownames(PS.l[[1]]@otu_table)==rownames(sdt)

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
seqs <- getSequencesFromTable(MA)
seqs <- lapply(seqs, DNAStringSet)

taxa_rescript <- list()

for (i in 1:3){
    taxa_rescript[[i]] <- assignTaxonomy(seqs[[i]],
                                         "/SAN/Susanas_den/AmpMarkers/RESCRIPt/SSURef_NR99/Fastas/Slv138.dada2.fa",
                                    multithread=90,
                                    tryRC = TRUE,
                                    verbose=TRUE)
}

taxa_rescript

MA2 <- MA

MA2@taxonTable <- taxa_rescript

#for (i in 1:3){
#    colnames(MA2@taxonTable[[i]]) <- c("domain", "kingdom", "phylum",  "class", "order", "family", "genus", "species")
#}


##To phyloseq
##Sample data
# temporary fix
source("bin/toPhyloseq.R")
PS2 <- TMPtoPhyloseq(MA2, colnames(MA2))
##Primer data
PS2.l <- TMPtoPhyloseq(MA2, colnames(MA2),  multi2Single=FALSE)

PS2@tax_table[,7]

#there's 2 samples without metadata
PS2 <- subset_samples(PS2, rownames(PS2@sam_data)%in%rownames(sdt))

#sanity check
rownames(PS2@otu_table) == rownames(sdt)

# adding sample data slot
PS2@sam_data <- sample_data(sdt)
# sanity check
rownames(PS2@sam_data)==rownames(PS2@otu_table)
saveRDS(PS2, file="/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/PhyloSeqData18S_SILVA.Rds")

##Primer data
for (i in 1:3) {
    PS2.l[[i]] <- subset_samples(PS2.l[[i]], rownames(PS2.l[[i]]@sam_data)%in%rownames(sdt))
}

#sanity check
rownames(PS2.l[[1]]@otu_table)==rownames(sdt)

# adding sample data
for (i in 1:3)
{
    sample_data(PS2.l[[i]]) <- sdt
}
saveRDS(PS2.l, file="/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/PhyloSeqList18S_SILVA.Rds")



