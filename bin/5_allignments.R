#!/usr/bin/Rscript

library(ggplot2)
library(dada2)
#library(MultiAmplicon, lib.loc="/usr/local/lib/R/site-library/")
library(reshape)
library(phyloseq)
library(DECIPHER)
library(phangorn)

#remotes::install_github("vmikk/metagMisc")
## using the devel
devtools::load_all("/SAN/Susanas_den/MultiAmplicon/")

PS <- readRDS("tmp/PhyloSeqData_All.Rds")
PS.l <- readRDS("tmp/PhyloSeqList_All.Rds")
PS18S <- readRDS("tmp/PS_18Swang.Rds")
PSwang <- PS.l[[37]]

trainingSet <- readRDS("/SAN/Susanas_den/AmpMarkers/wildEimeria18S/EimrefTrainingSet.RDS")

trainingSet

source("bin/PlottingCor.R")
# let's filter
fPS18S <- fil(PS18S)
fPS <- fil(PS)
fPSwang <- fil(PSwang)

## OK, now we want all the Eimeria sequences

Eim <- subset_taxa(fPS, family%in%"Eimeriidae")
Eim2 <- subset_taxa(fPS18S, family%in%"Eimeriidae")

seqs <- DNAStringSet(getSequences(colnames(Eim@otu_table)))
seqs2 <- DNAStringSet(getSequences(colnames(Eim2@otu_table)))

MA_dada2Sp <- dada2::assignSpecies(seqs2, "/SAN/Susanas_den/AmpMarkers/wildEimeria18S/Eim_ref.fa")
SA_dada2Sp <- dada2::assignSpecies(seqs, "/SAN/Susanas_den/AmpMarkers/wildEimeria18S/Eim_ref.fa", allowMultiple=TRUE)

SA_dada2Tx <- dada2::assignTaxonomy(seqs2, "/SAN/Susanas_den/AmpMarkers/wildEimeria18S/Eim_refAssignTaxonomy.fa")
MA_dada2Tx <- dada2::assignTaxonomy(seqs, "/SAN/Susanas_den/AmpMarkers/wildEimeria18S/Eim_refAssignTaxonomy.fa")

MA_IDtax <- IdTaxa(seqs,
                 trainingSet,
                 strand = "both",
                 threshold = 40,
                 bootstraps = 100,
                 processors = NULL,
                 verbose = TRUE,
                 type = "extended")

SA_IDtax <- IdTaxa(seqs2,
                 trainingSet,
                 strand = "both",
                 threshold = 40,
                 bootstraps = 100,
                 processors = NULL,
                 verbose = TRUE,
                 type = "extended")


### check out all taxonomic annotations

rownames(MA_dada2Sp) <- NULL
rownames(SA_dada2Sp) <- NULL
rownames(MA_dada2Tx) <- NULL
rownames(SA_dada2Tx) <- NULL

MA_dada2Sp
SA_dada2Sp

MA_dada2Tx
SA_dada2Tx

sapply(SA_IDtax, function(x) {
    x$taxon[8]})

sapply(MA_IDtax, function(x) {
    x$taxon[8]})

##aligments
refEim <- readDNAStringSet("/SAN/Susanas_den/AmpMarkers/wildEimeria18S/Eim_ref.fa")
poolSeqs <- c(seqs, seqs2)
alignment <- AlignSeqs(poolSeqs, anchor=NA, verbose=FALSE)
alignmentdb <- AlignSeqs(refEim, anchor=NA, verbose=FALSE)
#BrowseSeqs(alignment, highlight=0)

max(DistanceMatrix(alignment))

alignment

# making a phylo tree
PhanAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(PhanAlign)
treeUPGMA  <- upgma(dm)
treeNJ <- NJ(dm)

fun <- function(x) upgma(dist.ml(x))
bs_upgma <- bootstrap.phyDat(PhanAlign, fun)
plotBS(treeUPGMA, bs_upgma, main="UPGMA")
#ml
fit <- pml(treeNJ, data=PhanAlign)

# now for fun, our reference eimeria db

Peim <- phyDat(as(alignmentdb, "matrix"), type="DNA")
dm <- dist.ml(Peim)
treeUPGMA  <- upgma(dm)
treeNJ <- NJ(dm)

fun <- function(x) upgma(dist.ml(x))
bs_upgma <- bootstrap.phyDat(Peim, fun)
plotBS(treeUPGMA, bs_upgma, main="UPGMA")

                                        # for NJ too
fun <- function(x) NJ(dist.ml(x))
bs_NJ <- bootstrap.phyDat(Peim, fun)
plotBS(treeNJ, bs_NJ, main="NJ")

### ML

fit <- pml(treeNJ, data=Peim)
fit
