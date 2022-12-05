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
                  truncLen=c(200,200), 
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
                  verbose=0, mc.cores=90)
  MA <- mergeMulti(MA, mc.cores=90) 
  propMerged <- MultiAmplicon::calcPropMerged(MA)
  MA <- makeSequenceTableMulti(MA, mc.cores=90)
  MA <- removeChimeraMulti(MA, mc.cores=90)  
  saveRDS(MA, "/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/MA_multi_testrun_1.RDS")
} else{
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
                  truncLen=c(200,200), 
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


primerL <- read.table("/SAN/Susanas_den/gitProj/HMHZ/data/primerInputUnique.csv", head=T, sep=",")
ptable$Primer_name <- paste(ptable$corrected.NameF, ptable$corrected.NameR, sep=".")

which(!ptable$Primer_name %in% primerL$Primer_name)

## some primers are not on the list 
Primer_name <- c(ptable$Primer_name[6],
ptable$Primer_name[24],
ptable$Primer_name[40],
ptable$Primer_name[44],
ptable$Primer_name[45])
Target <- c("Eimeria", "Bacteria", "Eimeria", "Eimeria", "Eimeria")
Gen <- c("COX0", "16S", "COX1", "EfaB_31746", "EfaB_31746")

Pdf <- data.frame(Primer_name, Gen, Target)

# manual correction
#ptable$Primer_name[16] <- "16S.1100.F16_100_F.1492R_100_R"
#ptable$Primer_name[29] <- "NLF184cw_74_F.NL818cw_74_R"

primerL$Primer_name[4] <- "16S.1100.F16_100_for.1492R_100_rev"
primerL$Primer_name[122] <- "27M_F_98_F.Klin0341_CR_18_R"
primerL$Primer_name[70] <-  "LSU_Fwd_2_3Mod_55_F.LSU_R_4_54_R"
primerL$Primer_name[6] <- "18S_0067a_deg_3Mod_53_F.NSR399_3Mod_53_R"
primerL$Primer_name[7] <-  "18S_0067a_deg_5Mod_52_F.NSR399_5Mod_52_R"
primerL$Primer_name[120] <- "Bgf_132_F.Bgr_132_R"
primerL$Primer_name[86] <- "NLF184cw _74_F.NL818cw_74_R"


p.df <- primerL[which(primerL$Primer_name %in%ptable$Primer_name),c(2, 10, 11)]
p.df <- rbind(p.df, Pdf)
p.df <- p.df[match(names(MA1@PrimerPairsSet), p.df$Primer_name),]

which(ptable$Primer_name == p.df$Primer_name)
names(MA1@PrimerPairsSet)==p.df$Primer_name
rownames(p.df) <- seq(1,48,1)
table(p.df$Gen, p.df$Target)
p.df$Gen

taxT1 <- list()
seqs <- getSequencesFromTable(MA1)
seqs <- lapply(seqs, DNAStringSet)
for (i in 1:48){
    if (p.df$Gen[i]=="16S"){
        try(taxT1[[i]] <- assignTaxonomy(seqs[[i]],
   "/SAN/Susanas_den/AmpMarkers/RESCRIPt/SSURef_NR99/Fastas/Slv138.dada2.fa",
          multithread=90,
                                    tryRC = TRUE,
                                   verbose=TRUE))
    }
    else if (p.df$Gen[i]=="18S"){
        try(taxT1[[i]] <- assignTaxonomy(seqs[[i]],
   "/SAN/Susanas_den/AmpMarkers/RESCRIPt/SSURef_NR99/Fastas/Slv138.dada2.fa",
                                     multithread=90,
                                    tryRC = TRUE,
                                    verbose=TRUE))
    }
    else if (p.df$Gen[i]=="28S"){
        try(taxT1[[i]] <- assignTaxonomy(seqs[[i]],
 "/SAN/Susanas_den/AmpMarkers/RESCRIPt/LSURef_NR99/Fastas/Slv138LSU.dada2.fa",
                                     multithread=90,
                                    tryRC = TRUE,
                                    verbose=TRUE))
    }   
    else if (p.df$Gen[i]=="ITS"){
     try(taxT1[[i]] <- assignTaxonomy(seqs[[i]],
  "/SAN/Susanas_den/AmpMarkers/UNITE/sh_general_release_s_all_10.05.2021/sh_general_release_dynamic_s_all_10.05.2021.fasta",
                                     multithread=90,
                                    tryRC = TRUE,
                                    verbose=TRUE))
    }
    else {
     try(taxT1[[i]] <- assignTaxonomy(seqs[[i]],
       "/SAN/Susanas_den/AmpMarkers/RESCRIPt/other/Fastas/other.dada2.fa",
                                     multithread=90,
                                    tryRC = TRUE,
                                    verbose=TRUE))
    }   
}
taxT2 <- list()
seqs <- getSequencesFromTable(MA2)
seqs <- lapply(seqs, DNAStringSet)
for (i in 1:48){
    if (p.df$Gen[i]=="16S"){
        try(taxT2[[i]] <- assignTaxonomy(seqs[[i]],
   "/SAN/Susanas_den/AmpMarkers/RESCRIPt/SSURef_NR99/Fastas/Slv138.dada2.fa",
          multithread=90,
                                    tryRC = TRUE,
                                   verbose=TRUE))
    }
    else if (p.df$Gen[i]=="18S"){
        try(taxT2[[i]] <- assignTaxonomy(seqs[[i]],
   "/SAN/Susanas_den/AmpMarkers/RESCRIPt/SSURef_NR99/Fastas/Slv138.dada2.fa",
                                    multithread=90,
                                    tryRC = TRUE,
                                    verbose=TRUE))
    }
    else if (p.df$Gen[i]=="28S"){
        try(taxT2[[i]] <- assignTaxonomy(seqs[[i]],
 "/SAN/Susanas_den/AmpMarkers/RESCRIPt/LSURef_NR99/Fastas/Slv138LSU.dada2.fa",
                                    multithread=90,
                                    tryRC = TRUE,
                                    verbose=TRUE))
    }   
    else if (p.df$Gen[i]=="ITS"){
     try(taxT2[[i]] <- assignTaxonomy(seqs[[i]],
  "/SAN/Susanas_den/AmpMarkers/UNITE/sh_general_release_s_all_10.05.2021/sh_general_release_dynamic_s_all_10.05.2021.fasta",
                                    multithread=90,
                                    tryRC = TRUE,
                                    verbose=TRUE))
    }
    else {
     try(taxT2[[i]] <- assignTaxonomy(seqs[[i]],
       "/SAN/Susanas_den/AmpMarkers/RESCRIPt/other/Fastas/other.dada2.fa",
                                    multithread=90,
                                    tryRC = TRUE,
                                    verbose=TRUE))
    }   
}


taxT2[[44]]

## Little inspection
taxa.print1 <- taxT1[[1]]
taxa.print2 <- taxT1[[5]]
taxa.print4 <- taxT1[[7]]

taxa.print5 <- taxT2[[1]]
taxa.print6 <- taxT2[[5]]
taxa.print7 <- taxT2[[7]]


rownames(taxa.print1) <- NULL
rownames(taxa.print2) <- NULL
rownames(taxa.print4) <- NULL
rownames(taxa.print5) <- NULL
rownames(taxa.print6) <- NULL
rownames(taxa.print7) <- NULL



taxa.print1[1:10,]
taxa.print2
taxa.print4[1:10,]

taxa.print5[1:10,]
taxa.print6[1:10,]
taxa.print7[1:10,]

MA1@taxonTable <- taxT1
MA2@taxonTable <- taxT2

saveRDS(MA1, "/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/MATax_TestRun_1_NewTax.Rds")
saveRDS(MA2, "/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/MATax_FullRun_1_NewTax.Rds")  



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

PS1.l <- TMPtoPhyloseq(MA1, colnames(MA1),  multi2Single=FALSE)
PS2.l <- TMPtoPhyloseq(MA2, colnames(MA2),  multi2Single=FALSE) 

#reordering metadata
sdt <- sdt[match(rownames(PS1@sam_data), sdt$labels),]

#sanity check
rownames(PS1@otu_table) == rownames(sdt)

# adding sample data slot
PS1@sam_data <- sample_data(sdt)

# sanity check
rownames(PS1@sam_data)==rownames(PS1@otu_table)
rownames(PS1@otu_table)==sample_names(PS1)
PS1@sam_data[which(!rownames(PS1@sam_data)==rownames(PS1@otu_table))]

rownames(PS1@sam_data) <- rownames(PS1@otu_table)

PS_neg <- subset_samples(PS1, grepl("NEGATIVE",rownames(PS1@otu_table)))

PS1@sam_data$Control <- FALSE
PS1@sam_data$Control[which(sample_names(PS1)%in%sample_names(PS_neg))] <- TRUE
# sanity check
PS1@sam_data$labels[PS1@sam_data$Control==FALSE]
rownames(PS1@sam_data)[PS1@sam_data$Control==TRUE]
library("decontam")
###### removing contaminants
## assuming that negative controls have 0 DNA
PS1@sam_data$Total_DNA[PS1@sam_data$Control==TRUE] <- 0.0001
## ----see-depths---------------------------------------------------------------
#df <- as.data.frame(sample_data(PS1)) # Put sample_data into a ggplot-friendly data.frame
#df$LibrarySize <- sample_sums(PS1)
#df <- df[order(df$LibrarySize),]
#df$Index <- seq(nrow(df))
#ggplot(data=df, aes(x=Index, y=LibrarySize, color=Control)) + geom_point()
ps <- phyloseq::prune_samples(sample_sums(PS1)>0, PS1)
contamdf.freq <- isContaminant(ps, method="either", conc="Total_DNA", neg="Control", threshold=c(0.1,0.5), normalize=TRUE)
table(contamdf.freq$contaminant)
### taxa to remove
ps@tax_table[rownames(contamdf.freq[contamdf.freq$contaminant==TRUE,]),5]
## let's remove them now and negative controls
Keep <- rownames(contamdf.freq[contamdf.freq$contaminant==FALSE,])
PS1 <- prune_samples(sample_data(PS1)$Control == FALSE, PS1)
PS1 <- prune_taxa(Keep, PS1)  
#plot_bar(PS_neg, fill="phylum")

## adding metadata, removing contaminants and controls
pos <- sample_names(subset_samples(PS1.l[[1]], !grepl("NEGATIVE",rownames(PS1.l[[1]]@otu_table))))

for (i in 1:48) {
    try(PS1.l[[i]]@sam_data <- sample_data(sdt), silent=TRUE)
    try(rownames(PS1.l[[i]]@sam_data) <- rownames(PS1.l[[i]]@otu_table), silent=TRUE)
}

for (i in 1:48) {
    try(PS1.l[[i]] <- prune_taxa(Keep, PS1.l[[i]]), silent=TRUE)
    try(PS1.l[[i]] <- prune_samples(pos, PS1.l[[i]]), silent=TRUE)
}

##### for PS2
# I know it's annoying, but I want to remove contaminants before merging

#reordering metadata
sdt <- sdt[match(rownames(PS2@sam_data), sdt$labels),]

#sanity check
rownames(PS2@otu_table) == rownames(sdt)

# adding sample data slot
PS2@sam_data <- sample_data(sdt)

# sanity check
rownames(PS2@sam_data)==rownames(PS2@otu_table)
rownames(PS2@otu_table)==sample_names(PS2)
PS2@sam_data[which(!rownames(PS2@sam_data)==rownames(PS2@otu_table))]
rownames(PS2@sam_data) <- rownames(PS2@otu_table)

PS_neg <- subset_samples(PS2, grepl("NEGATIVE",rownames(PS2@otu_table)))

PS2@sam_data$Control <- FALSE
PS2@sam_data$Control[which(sample_names(PS2)%in%sample_names(PS_neg))] <- TRUE
# sanity check
PS2@sam_data$labels[PS2@sam_data$Control==FALSE]
rownames(PS2@sam_data)[PS2@sam_data$Control==TRUE]


###### removing contaminants
## assuming that negative controls have 0 DNA
PS2@sam_data$Total_DNA[PS2@sam_data$Control==TRUE] <- 0.0001

## ----see-depths---------------------------------------------------------------
#df <- as.data.frame(sample_data(PS2)) # Put sample_data into a ggplot-friendly data.frame
#df$LibrarySize <- sample_sums(PS2)
#df <- df[order(df$LibrarySize),]
#df$Index <- seq(nrow(df))
#ggplot(data=df, aes(x=Index, y=LibrarySize, color=Control)) + geom_point()

ps <- phyloseq::prune_samples(sample_sums(PS2)>0, PS2)
contamdf.freq <- isContaminant(ps, method="either", conc="Total_DNA", neg="Control", threshold=c(0.1,0.5), normalize=TRUE)
table(contamdf.freq$contaminant)
### taxa to remove
ps@tax_table[rownames(contamdf.freq[contamdf.freq$contaminant==TRUE,]),5]
## let's remove them now and negative controls
Keep <- rownames(contamdf.freq[contamdf.freq$contaminant==FALSE,])
PS2 <- prune_samples(sample_data(PS2)$Control == FALSE, PS2)
PS2 <- prune_taxa(Keep, PS2)  
#plot_bar(PS_neg, fill="phylum")

## adding metadata, removing contaminants and controls
pos <- sample_names(subset_samples(PS2.l[[1]], !grepl("NEGATIVE",rownames(PS2.l[[1]]@otu_table))))

for (i in 1:48) {
    try(PS2.l[[i]]@sam_data <- sample_data(sdt), silent=TRUE)
    try(rownames(PS2.l[[i]]@sam_data) <- rownames(PS2.l[[i]]@otu_table), silent=TRUE)
}


for (i in 1:48) {
    try(PS2.l[[i]] <- prune_taxa(Keep, PS2.l[[i]]), silent=TRUE)
    try(PS2.l[[i]] <- prune_samples(pos, PS2.l[[i]]), silent=TRUE)
}

PS <- merge_phyloseq(PS1, PS2) ###Works!

#saveRDS(PS, file="/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/PhyloSeqData_All.Rds") ###Results from full + test run 

saveRDS(PS, file="/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/PhyloSeqData_All_Tax_New.Rds") ###Results from full + test run 


names(PS1.l)== names(PS2.l)

along<- names(PS2.l) ## Run with less primers working
# 7 amplicons are empty.

PS.l <- lapply(along, function(i) try(merge_phyloseq(PS1.l[[i]], PS2.l[[i]]))) ##Merge all the information from both experiments
names(PS.l) <- names(PS2.l) ###Use the names from test list

#sanity check
rownames(PS.l[[1]]@otu_table)%in%rownames(sdt)

saveRDS(PS.l, file="/SAN/Susanas_den/gitProj/LabMicrobiome/tmp/PhyloSeqList_All_Tax_New.Rds") ###For primer analysis

