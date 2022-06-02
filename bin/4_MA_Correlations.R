#!/usr/bin/Rscript

library(ggplot2)
library(dada2)
#library(MultiAmplicon, lib.loc="/usr/local/lib/R/site-library/")
library(reshape)
library(phyloseq)
#remotes::install_github("vmikk/metagMisc")
library(metagMisc)
require(grid)
require(grid_extra)
require(cowplot)
library(microbiome)
library(DESeq2)
## using the devel
#devtools::load_all("/SAN/Susanas_den/MultiAmplicon/")

PS <- readRDS("tmp/PhyloSeqData_All.Rds")
PS.l <- readRDS("tmp/PhyloSeqList_All.Rds")
PS18S <- readRDS("tmp/PhyloSeqData18S.Rds")
PSwang <- PS.l[[31]]

sam <- data.frame(sample_data(PS))

#let's remove samples with no reads
PS <- prune_samples(sample_sums(PS)>0, PS)
PS18S <- prune_samples(sample_sums(PS18S)>0, PS18S)
PSwang <- prune_samples(sample_sums(PSwang)>0, PSwang)
source("bin/PlottingCor.R")

Plotting_cor(PS, "MA", dir="fig/MA/")

Plotting_cor(PS18S, "SA", dir="fig/SA/")

for (i in 1:38) {
    nm <- names(PS.l)[i]
    ps <- PS.l[[i]]
    try(Plotting_cor(PS=ps, name=nm, dir="fig/MA/"))
}

for (i in 1:38) {
    nm <- names(PS.l)[i]
    ps <- PS.l[[i]]
    try(NoFilPlotting_cor(PS=ps, name=nm, dir="fig/MA/NoFil/"))
}

for (i in 1:38) {
#    print(names(PS.l)[i])
    try(p <- subset_taxa(PS.l[[i]],phylum=="Apicomplexa"), silent=TRUE)
    try(get_taxa_unique(p, "family"), silent=TRUE)
    if (exists("p")) {
        a <- get_taxa_unique(p, "family")
        print(paste(i, "- ", names(PS.l[i]), ": ", length(a), sep=""))
        print(a)
}
    rm(p)
}

#######################
