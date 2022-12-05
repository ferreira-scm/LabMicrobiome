#!/usr/bin/Rscript

library(ggplot2)
#library(dada2)
#library(MultiAmplicon, lib.loc="/usr/local/lib/R/site-library/")
#library(reshape)
#remotes::install_github("vmikk/metagMisc")
library(metagMisc)
require(grid)
library(Hmisc)
library(phyloseq)
#library(grid_extra)
require(cowplot)
library(RColorBrewer)
library(dplyr)
## using the devel
#devtools::load_all("/SAN/Susanas_den/MultiAmplicon/")
source("bin/PlottingCor.R")

## All runs pooled
all.PS <- readRDS("tmp/PhyloSeqData_All.Rds")

all.PS.slv <- readRDS("tmp/PhyloSeqData_All_Tax_New.Rds")

all.PS.l <- readRDS("tmp/PhyloSeqList_All.Rds")
all.PS.l.slv <- readRDS("tmp/PhyloSeqList_All_Tax_New.Rds")

# single amplicon from MA run
all.PSwang <- all.PS.l[[37]]

## Single amplicon 18S
sin.PS18S <- readRDS("tmp/PhyloSeqList18S.Rds")
sin.PS18S.slv <- readRDS("tmp/PhyloSeqList18S_SILVA.Rds")
sin.PS18S <- sin.PS18S[[2]]
sin.PS18S.slv <- sin.PS18S.slv[[2]]

## Single amplicon "pooled"
sin.PS <- readRDS("tmp/PhyloSeqData18S.Rds")
sin.PS.slv <- readRDS("tmp/PhyloSeqData18S_SILVA.Rds")




# let's filter
f.sin18 <- fil(sin.PS18S)
f.sin18.slv <- fil(sin.PS18S.slv)
f.all <- fil(all.PS)
f.all.slv <- fil(all.PS.slv)



## filtering MA by amplicon
f.all.l <- list()
for (i in 1:48) {
    try(f.all.l[[i]] <- fil(all.PS.l[[i]]), silent=TRUE)
}

f.all.l

f.all.lp <- f.all.l[[1]]
for (i in 2:47){
    f.all.lp <- try(merge_phyloseq(f.all.lp,f.all.l[[i]]))
    }

f.all.l.slv <- list()
for (i in 1:48) {
    try(f.all.l.slv[[i]] <- fil(all.PS.l.slv[[i]]), silent=TRUE)
}

f.all.lp.slv <- f.all.l.slv[[1]]
for (i in 2:47){
    f.all.lp.slv <- try(merge_phyloseq(f.all.lp.slv,f.all.l.slv[[i]]))
    }



# sanity check
#merge_phyloseq(f.all.l[[45]], f.all.l[[44]], f.all.l[[47]], f.all.l[[37]],f.all.l[[1]], f.all.l[[7]],
#               f.all.l[[12]], f.all.l[[13]], f.all.l[[13]], f.all.l[[22]], f.all.l[[23]], f.all.l[[28]],
#               f.all.l[[33]], f.all.l[[35]])
## Single amplicon "pooled"
f.sin <- fil(sin.PS)
f.sin.slv <- fil(sin.PS.slv)

# and transform
sin18TSS  <- transform_sample_counts(f.sin18, function(x) x / sum(x))
#sinTSS <- transform_sample_counts(f.sin, function(x) x / sum(x))
allTSS <-  transform_sample_counts(f.all.lp, function(x) x / sum(x))
sin18TSS.slv  <- transform_sample_counts(f.sin18.slv, function(x) x / sum(x))
#sinTSS <- transform_sample_counts(f.sin, function(x) x / sum(x))
allTSS.slv <-  transform_sample_counts(f.all.lp.slv, function(x) x / sum(x))

T.sin18 <- sin18TSS
T.all <- allTSS
T.sin18.slv <- sin18TSS.slv
T.all.slv <- allTSS.slv
otu_table(T.sin18) <- otu_table(T.sin18)*sample_data(T.sin18)$Total_DNA
otu_table(T.all) <- otu_table(T.all)*sample_data(T.all)$Total_DNA
otu_table(T.sin18.slv) <- otu_table(T.sin18.slv)*sample_data(T.sin18.slv)$Total_DNA
otu_table(T.all.slv) <- otu_table(T.all.slv)*sample_data(T.all.slv)$Total_DNA


#### exploring MA
#for (i in 1:48) {
#    nm <- names(all.PS.l)[i]
#    ps <- all.PS.l[[i]]
#    print(nm)
#    try(Plotting_cor(ps, name=nm, dir="fig/MA/"))
#}

#for (i in 1:48) {
#    nm <- names(all.PS.l)[i]
#    ps <- all.PS.l[[i]]
#    try(NoFilPlotting_cor(ps, name=nm, dir="fig/MA/NoFil/"))
#}

rank_names(all.PS.l.slv[[1]])

rank_names(T.sin18.slv)

get_taxa_unique(T.sin18.slv, "Phylum")

get_taxa_unique(T.sin18, "genus")

rank_names(T.sin18.slv)

get_taxa_unique(f.all.lp.slv, "Phylum")




#########################################################
#how many primers amplify Apicomplexa and which families?
for (i in 1:48) {
#    print(names(all.PS.l)[[i]])
    try(p <- subset_taxa(all.PS.l.slv[[i]],Phylum=="p__Apicomplexa"), silent=TRUE)
#    try(get_taxa_unique(p, "family"), silent=TRUE)
    if (exists("p")) {
        a <- get_taxa_unique(p, "Family")
        print(paste(i, "- ", names(all.PS.l.slv[i]), ": ", length(a), sep=""))
        print(a)
}
    rm(p)
}

##### and how many amplicons have eimeria?
for (i in 1:48) {
#    print(names(all.PS.l)[i])
    try(p <- subset_taxa(all.PS.l.slv[[i]],Genus=="g__Eimeria"), silent=TRUE)
#    try(get_taxa_unique(p, "genus"), silent=TRUE)
    if (exists("p")) {
        a <- get_taxa_unique(p, "Genus")
        print(paste(i, "- ", names(all.PS.l.slv[i]), ": ", nrow(p@tax_table), sep=""))
        print(a)
}
    rm(p)
}


### and what happens when we filter?
for (i in 1:48) {
#    print(names(all.PS.l)[i])
    try(p <- subset_taxa(f.all.l.slv[[i]],Genus=="g__Eimeria"), silent=TRUE)
#    try(get_taxa_unique(p, "genus"), silent=TRUE)
    if (exists("p")) {
        a <- get_taxa_unique(p, "Genus")
        print(paste(i, "- ", names(f.all.l.slv[i]), ": ", nrow(p@tax_table), sep=""))
        print(a)
}
    rm(p)
}

names(all.PS.l.slv)

## how many eimeria ASV reads do we have in single amplicon before and after filtering
subset_taxa(sin.PS.slv, Genus=="g__Eimeria")

# sanity check
subset_taxa(sin.PS18S.slv, Genus=="g__Eimeria")

# after filtering
subset_taxa(f.sin.slv, Genus=="g__Eimeria")

#sanity check
subset_taxa(f.sin18.slv, Genus=="g__Eimeria")

## OK, now we want all the Eimeria sequences
#Eim2 <- subset_taxa(T.sin18, genus%in%"Eimeria")
#Eim <- subset_taxa(T.all, genus%in%"Eimeria")

Eim2 <- subset_taxa(T.sin18.slv, Genus%in%"g__Eimeria")
Eim <- subset_taxa(T.all.slv, Genus%in%"g__Eimeria")
