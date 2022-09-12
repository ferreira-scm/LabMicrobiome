#!/usr/bin/Rscript

library(ggplot2)
#library(dada2)
#library(MultiAmplicon, lib.loc="/usr/local/lib/R/site-library/")
#library(reshape)
#remotes::install_github("vmikk/metagMisc")
library(metagMisc)
require(grid)
library(Hmisc)

#library(grid_extra)
require(cowplot)
library(RColorBrewer)
library(dplyr)
## using the devel
#devtools::load_all("/SAN/Susanas_den/MultiAmplicon/")
source("bin/PlottingCor.R")

## All runs pooled
all.PS <- readRDS("tmp/PhyloSeqData_All.Rds")
all.PS.l <- readRDS("tmp/PhyloSeqList_All.Rds")
# single amplicon from MA run
all.PSwang <- all.PS.l[[37]]

## Single amplicon 18S
sin.PS18S <- readRDS("tmp/PS_18Swang.Rds")
sin.PS18S.slv <- readRDS("tmp/PS_18Swang_SILVA.Rds")

## Single amplicon "pooled"
sin.PS <- readRDS("tmp/PhyloSeqData18S.Rds")
sin.PS.slv <- readRDS("tmp/PhyloSeqData18S_SILVA.Rds")

# let's filter
f.sin18 <- fil(sin.PS18S)
f.all <- fil(all.PS)
f.allwang <- fil(all.PSwang)
f.sin18.slv <- fil(sin.PS18S.slv)
## filtering MA by amplicon
f.all.l <- list()
for (i in 1:48) {
    try(f.all.l[[i]] <- fil(all.PS.l[[i]]), silent=TRUE)
}

f.all.lp <- f.all.l[[1]]
for (i in 2:47){
    f.all.lp <- try(merge_phyloseq(f.all.lp,f.all.l[[i]]))
    }


# sanity check
#merge_phyloseq(f.all.l[[45]], f.all.l[[44]], f.all.l[[47]], f.all.l[[37]],f.all.l[[1]], f.all.l[[7]],
#               f.all.l[[12]], f.all.l[[13]], f.all.l[[13]], f.all.l[[22]], f.all.l[[23]], f.all.l[[28]],
#               f.all.l[[33]], f.all.l[[35]])
## Single amplicon "pooled"
f.sin <- fil(sin.PS)
f.sin.slv <- fil(sin.PS.slv)
# and transform
T.sin18 <- f.sin18
T.all <- f.all
T.all.l <- f.all.lp
T.allwang <-f.allwang
T.sin18.slv <- f.sin18.slv
T.sin <-f.sin
T.sin.slv <-f.sin.slv
otu_table(T.sin18) <- otu_table(T.sin18)*sample_data(T.sin18)$DNA_g_feces
otu_table(T.all) <- otu_table(T.all)*sample_data(T.all)$DNA_g_feces
otu_table(T.all.l) <- otu_table(T.all.l)*sample_data(T.all.l)$DNA_g_feces
otu_table(T.allwang) <- otu_table(T.allwang)*sample_data(T.allwang)$DNA_g_feces
otu_table(T.sin18.slv) <- otu_table(T.sin18.slv)*sample_data(T.sin18.slv)$DNA_g_feces
otu_table(T.sin) <- otu_table(T.sin)*sample_data(T.sin)$DNA_g_feces
otu_table(T.sin.slv) <- otu_table(T.sin.slv)*sample_data(T.sin.slv)$DNA_g_feces

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


#########################################################
#how many primers amplify Apicomplexa and which families?
for (i in 1:48) {
#    print(names(all.PS.l)[[i]])
    try(p <- subset_taxa(all.PS.l[[i]],phylum=="Apicomplexa"), silent=TRUE)
    try(get_taxa_unique(p, "family"), silent=TRUE)
    if (exists("p")) {
        a <- get_taxa_unique(p, "family")
        print(paste(i, "- ", names(all.PS.l[i]), ": ", length(a), sep=""))
        print(a)
}
    rm(p)
}

##### and how many amplicons have eimeria?
for (i in 1:48) {
#    print(names(all.PS.l)[i])
    try(p <- subset_taxa(all.PS.l[[i]],genus=="Eimeria"), silent=TRUE)
#    try(get_taxa_unique(p, "genus"), silent=TRUE)
    if (exists("p")) {
        a <- get_taxa_unique(p, "genus")
        print(paste(i, "- ", names(all.PS.l[i]), ": ", nrow(p@tax_table), sep=""))
        print(a)
}
    rm(p)
}


##### and how many amplicons have eimeria?
for (i in 1:48) {
#    print(names(all.PS.l)[i])
    try(p <- subset_taxa(all.PS.l[[i]],genus=="Eimeria"), silent=TRUE)
#    try(get_taxa_unique(p, "genus"), silent=TRUE)
    if (exists("p")) {
        a <- get_taxa_unique(p, "genus")
        print(paste(i, "- ", names(all.PS.l[i]), ": ", nrow(p@tax_table), sep=""))
        print(a)
}
    rm(p)
}


### and what happens when we filter?
for (i in 1:48) {
#    print(names(all.PS.l)[i])
    try(p <- subset_taxa(f.all.l[[i]],genus=="Eimeria"), silent=TRUE)
#    try(get_taxa_unique(p, "genus"), silent=TRUE)
    if (exists("p")) {
        a <- get_taxa_unique(p, "genus")
        print(paste(i, "- ", names(f.all.l[i]), ": ", nrow(p@tax_table), sep=""))
        print(a)
}
    rm(p)
}


## how many eimeria ASV reads do we have in single amplicon before and after filtering
subset_taxa(sin.PS, genus=="Eimeria")
# sanity check
subset_taxa(sin.PS18S, genus=="Eimeria")
# after filtering
subset_taxa(f.sin, genus=="Eimeria")
#sanity check
subset_taxa(f.sin18, genus=="Eimeria")

## OK, now we want all the Eimeria sequences
Eim <- subset_taxa(T.all.l, genus%in%"Eimeria")
Eim2 <- subset_taxa(T.sin18, genus%in%"Eimeria")
#Eim.slv <- subset_taxa(T.sin18.slv, Family%in%"Eimeriorina")
Eim_nf <- subset_taxa(all.PS, genus%in%"Eimeria")
Eim_nf_wang <- subset_taxa(all.PS.l[[37]], genus%in%"Eimeria")
Eim2_nf <- subset_taxa(sin.PS18S, genus%in%"Eimeria")

