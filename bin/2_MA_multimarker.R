#!/usr/bin/Rscript

library(ggplot2)
library(dada2)
#library(MultiAmplicon, lib.loc="/usr/local/lib/R/site-library/")
library(reshape)
library(phyloseq)
library(data.table, lib.loc="/usr/local/lib/R/site-library/")
library(taxonomizr)
library(taxize)
library(parallel)

## using the devel
devtools::load_all("/SAN/Susanas_den/MultiAmplicon/")

## re-run or use pre-computed results for different parts of the pipeline:
doFilter <- FALSE
doMultiAmp <- FALSE
doTax <- FALSE

###################Test run Microbiome######################
#Preparation of files
##These are the same steps that are followed by the DADA2 pipeline
path <- "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/2018_22_Eie_TestRun/" ## Test run 24.06.2020
#path <- "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/2018_22_Eie_FullRun_1/" ## Full run 29.06.2020
fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) #take all fastaq files from the folder 
fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE) #separate the forward reads
fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE) #separate the reverse reads 

samples <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", basename(fastqF))
samples<- gsub("S\\d+-", "\\1", basename(samples))
samples<-gsub("-", "_", basename(samples))

#Creation of a folder for filtrated reads 

filt_path <- "/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/TestRun/filtered_TestRun_1/"
#Pipeline filtration 
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs) <- samples
filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs) <- samples

#Quality plots of the reads
pdf("/SAN/Susanas_den/gitProj/LabMicrobiome/fig/quality/qualityProfileF1_testrun1.pdf",
    height = 7, width = 7)
plotQualityProfile(fastqF[[1]])
dev.off()

#Quality plots of the reads
pdf("/SAN/Susanas_den/gitProj/LabMicrobiome/fig/quality/qualityProfileR20_testrun1.pdf",
    height = 7, width = 7)
plotQualityProfile(fastqR[[20]])
dev.off()

## some files will be filtered out completely, therefore allowing 50
## files less present and still don't redo filtering # full run 150 170 test run 150 150
if(doFilter){
  lapply(seq_along(fastqF),  function (i) {
    filterAndTrim(fastqF[i], filtFs[i], fastqR[i], filtRs[i],
                  truncLen=c(150,170), 
                  maxN=0, maxEE=2, truncQ=2, rm.phix = TRUE,
                  compress=TRUE, verbose=TRUE)
  })
}

names(filtFs) <- names(filtRs) <- samples
files <- PairedReadFileSet(filtFs, filtRs)

#Preparation of primer file ### Here stats the Multiamplicon pipeline from Emanuel

#Primers used in the arrays 
ptable <- read.csv(file = "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/primer.file.multi.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
primerF <- ptable[, "TS.SequenceF"]
primerR <- ptable[, "TS.SequenceR"]
names(primerF) <- as.character(ptable[, "corrected.NameF"])
names(primerR) <- as.character(ptable[, "corrected.NameR"])

primer <- PrimerPairsSet(primerF, primerR)

##Multi amplicon pipeline
if(doMultiAmp){
  MA <- MultiAmplicon(primer, files)
  #filedir <- "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/Stratified_files_TestRun"
  #filedir <- "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/Stratified_files_FullRun_1"
  filedir <- "/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/TestRun/stratified_filea_test_run_1"
  if(dir.exists(filedir)) unlink(filedir, recursive=TRUE)
  MA <- sortAmplicons(MA, n=1e+05, filedir=filedir) ## This step sort the reads into amplicons based on the number of primer pairs
  errF <-  learnErrors(unlist(getStratifiedFilesF(MA)), nbase=1e8,
                       verbose=0, multithread = 90)
  errR <- learnErrors(unlist(getStratifiedFilesR(MA)), nbase=1e8,
                      verbose=0, multithread = 90)
  MA <- dadaMulti(MA, Ferr=errF, Rerr=errR,  pool=FALSE,
                  verbose=0, mc.cores=12)
  MA <- mergeMulti(MA, mc.cores=90) 
  propMerged <- MultiAmplicon::calcPropMerged(MA)
  MA <- mergeMulti(MA, mc.cores=90) 
  MA <- makeSequenceTableMulti(MA, mc.cores=90) 
  MA <- removeChimeraMulti(MA, mc.cores=90)
  #saveRDS(MA, "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/MA_Multi_TestRun.RDS")
  saveRDS(MA, "/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/MA_multi_testrun_1.RDS")
} else{
  #MA <- readRDS("/SAN/Victors_playground/Eimeria_microbiome/MA_Multi_TestRun.RDS")
  MA <- readRDS("/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/MA_multi_testrun_1.RDS")
}


err_F <- plotErrors(errF, nominalQ=TRUE)
pdf("fig/quality/TR_Estimeted_error_ratesF_1_1.pdf",
    height = 7, width = 7)
err_F
dev.off()

err_R <- plotErrors(errR, nominalQ=TRUE)
pdf("fig/quality/TR_Estimeted_error_ratesR_1_1.pdf",
    height = 7, width = 7)
err_R
dev.off()

Heatmap <- plotAmpliconNumbers(MA)
pdf("fig/quality/TR_heat_Sequencing_summary_HMHZ_1_1.pdf",
    height = 15, width = 15)
Heatmap
dev.off()

#trackingF <- getPipelineSummary(MA) 
#plotPipelineSummary(trackingF) 
#plotPipelineSummary(trackingF) + scale_y_log10()

#plotAmpliconNumbers(MA)

###New taxonomic assignment 
if(doTax){
MA <- blastTaxAnnot(MA,
                    db = "/SAN/db/blastdb/nt/nt",
                    negative_gilist = "/SAN/db/blastdb/uncultured.gi",
                    #infasta = "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/in_TestRun.fasta",
                    infasta = "/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/TestRun/in_testrun1",
                    #outblast = "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/out_TestRun.fasta",
                    outblast = "/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/TestRun/out_TestRun_1.fasta",
                    taxonSQL = "/SAN/db/taxonomy/taxonomizr.sql", 
                    num_threads = 90) ##Change for use more power!!
#saveRDS(MA, file="/SAN/Victors_playground/Eimeria_microbiome/Multimarker/MATax_TestRun.Rds")
saveRDS(MA, file="/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/MATax_TestRun_1.Rds")
} else {
    MA1 <- readRDS("/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/MATax_TestRun_1.Rds")
}

############################################################
###################Full run Microbiome######################
############################################################

#Preparation of files
##These are the same steps that are followed by the DADA2 pipeline
path <- "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/2018_22_Eie_FullRun_1"
fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) #take all fastaq files from the folder 
fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE) #separate the forward reads
fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE) #separate the reverse reads 

samples <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", basename(fastqF))
samples<- gsub("S\\d+-", "\\1", basename(samples))
samples<-gsub("-", "_", basename(samples))

#Extra step in the pipeline: quality plots of the reads 
## plotQualityProfile(fastqF[[1]])
## plotQualityProfile(fastqF[[2]])
## plotQualityProfile(fastqR[[1]])
## plotQualityProfile(fastqR[[2]])

#Creation of a folder for filtrated reads 
filt_path <- "/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/FullRun/filtered_FullRun_1/"

#Pipeline filtration 
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs) <- samples
filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs) <- samples

#Quality plots of the reads
pdf("/SAN/Susanas_den/gitProj/LabMicrobiome/fig/quality/qualityProfileF1_fullrun1.pdf",
    height = 7, width = 7)
plotQualityProfile(fastqF[[1]])
dev.off()

#Quality plots of the reads
pdf("/SAN/Susanas_den/gitProj/LabMicrobiome/fig/quality/qualityProfileR20_fullrun1.pdf",
    height = 7, width = 7)
plotQualityProfile(fastqR[[20]])
dev.off()

## some files will be filtered out completely, therefore allowing 50
## files less present and still don't redo filtering # full run 150 170 test run 150 150
if(doFilter){
  lapply(seq_along(fastqF),  function (i) {
    filterAndTrim(fastqF[i], filtFs[i], fastqR[i], filtRs[i],
                  truncLen=c(150,170), 
                  maxN=0, maxEE=2, truncQ=2, rm.phix = TRUE,
                  compress=TRUE, verbose=TRUE)
  })
}

names(filtFs) <- names(filtRs) <- samples
files <- PairedReadFileSet(filtFs, filtRs)

#Preparation of primer file ### Here stats the Multiamplicon pipeline from Emanuel

#Primers used in the arrays 
ptable <- read.csv(file = "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/primer.file.multi.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
primerF <- ptable[, "TS.SequenceF"]
primerR <- ptable[, "TS.SequenceR"]
names(primerF) <- as.character(ptable[, "corrected.NameF"])
names(primerR) <- as.character(ptable[, "corrected.NameR"])

primer <- PrimerPairsSet(primerF, primerR)

##Multi amplicon pipeline
if(doMultiAmp){
  MA <- MultiAmplicon(primer, files)
  filedir <- "/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/FullRun/stratified_files_full_run_1"
  if(dir.exists(filedir)) unlink(filedir, recursive=TRUE)
  MA <- sortAmplicons(MA, n=1e+05, filedir=filedir) ## This step sort the reads into amplicons based on the number of primer pairs
  errF <-  learnErrors(unlist(getStratifiedFilesF(MA)), nbase=1e8,
                       verbose=0, multithread = 90)
  errR <- learnErrors(unlist(getStratifiedFilesR(MA)), nbase=1e8,
                      verbose=0, multithread = 90)
  MA <- dadaMulti(MA, Ferr=errF, Rerr=errR,  pool=FALSE,
                  verbose=0, mc.cores=12)
  MA <- mergeMulti(MA, mc.cores=90) 
  propMerged <- MultiAmplicon::calcPropMerged(MA)
  MA <- mergeMulti(MA, mc.cores=90) 
  MA <- makeSequenceTableMulti(MA, mc.cores=90) 
  MA <- removeChimeraMulti(MA, mc.cores=90)
  #saveRDS(MA, "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/MA_Multi_TestRun.RDS")
  saveRDS(MA, "/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/MA_multi_fullrun_1.RDS")
} else{
  #MA <- readRDS("/SAN/Victors_playground/Eimeria_microbiome/MA_Multi_TestRun.RDS")
  MA <- readRDS("/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/MA_multi_fullrun_1.RDS")
}

err_F <- plotErrors(errF, nominalQ=TRUE)
pdf("fig/quality/FR_Estimeted_error_ratesF_1_1.pdf",
    height = 7, width = 7)
err_F
dev.off()

err_R <- plotErrors(errR, nominalQ=TRUE)
pdf("fig/quality/FR_Estimeted_error_ratesR_1_1.pdf",
    height = 7, width = 7)
err_R
dev.off()

Heatmap <- plotAmpliconNumbers(MA)
pdf("fig/quality/FR_heat_Sequencing_summary_HMHZ_1_1.pdf",
    height = 15, width = 15)
Heatmap
dev.off()

###New taxonomic assignment 
if(doTax){
MA <- blastTaxAnnot(MA,
                    db = "/SAN/db/blastdb/nt/nt",
                    negative_gilist = "/SAN/db/blastdb/uncultured.gi",
                    infasta = "/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/FullRun/in_fullrun1",
                    outblast = "/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/FullRun/out_FullRun_1.fasta",
                    taxonSQL = "/SAN/db/taxonomy/taxonomizr.sql", 
                    num_threads = 90) ##Change for use more power!!
saveRDS(MA, file="/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/MATax_FullRun_1.Rds")
} else {
    MA1 <- readRDS("/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/MATax_TestRun_1.Rds")
    MA2 <- readRDS("/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/MATax_FullRun_1.Rds")  
}

### Add sample information
if(!exists("sample.data")){
    source("/SAN/Susanas_den/gitProj/LabMicrobiome/bin/1_Data_preparation.R")
}

if(!exists("sdt")){
    source("/SAN/Susanas_den/gitProj/LabMicrobiome/bin/1_qPCR_data_preparation.R")
}

#little fix
rownames(sdt) <- sdt$labels
rownames(sdt)==rownames(sample.data)

##To phyloseq

source("bin/toPhyloseq.R")

PS1 <- TMPtoPhyloseq(MA1, colnames(MA1))
PS2<- TMPtoPhyloseq(MA2, colnames(MA2))

PS <- merge_phyloseq(PS1, PS2) ###Works!

PS@sam_data <- sample_data(sdt)

saveRDS(PS, file="/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/PhyloSeqData_All.Rds") ###Results from full + test run 


##Primer data
## just sorting out primers whithout any taxannot

##Makethe next function work 
#MA1 <- MA1[which(!unlist(lapply(MA1@taxonTable, is.null))), ] 
PS1.l <- TMPtoPhyloseq(MA1, colnames(MA1),  multi2Single=FALSE)
#MA2 <- MA2[which( !unlist(lapply(MA2@taxonTable, is.null))), ]
PS2.l <- TMPtoPhyloseq(MA2, colnames(MA2),  multi2Single=FALSE) 


names(PS1.l)== names(PS2.l)
# 7 amplicons are empty.
PS.l <- lapply(along, function(i) try(merge_phyloseq(PS1.l[[i]], PS2.l[[i]]))) ##Merge all the information from both experiments

names(PS.l) <- names(PS2.l) ###Use the names from test list


row.names(sdt) <- sdt$labels

row.names(sdt)==row.names(sample.data)

length(PS.l)

# adding sample data
for (i in 1:48)
{
    try(sam_data(PS.l[[i]]) <- sdt)
}

head(sam_data(PS.l[[1]]))

saveRDS(PS.l, file="/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/PhyloSeqList_All.Rds") ###For primer analysis

