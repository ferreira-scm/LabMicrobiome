#!/usr/bin/Rscript

library(ggplot2)
library(dada2)
#library(MultiAmplicon, lib.loc="/usr/local/lib/R/site-library/")
library(reshape)
library(phyloseq)
library(DECIPHER)
library(phangorn)
#library(BiocManager, lib.loc="/usr/local/lib/R/site-library/")

                                        #remotes::install_github("vmikk/metagMisc")
## using the devel
devtools::load_all("/SAN/Susanas_den/MultiAmplicon/")

devtools::load_all("/SAN/Susanas_den/phangorn/")


PS <- readRDS("tmp/PhyloSeqData_All.Rds")
PS.l <- readRDS("tmp/PhyloSeqList_All.Rds")
PS18S <- readRDS("tmp/PS_18Swang.Rds")
PSwang <- PS.l[[37]]

PSslv <- readRDS("tmp/PhyloSeqData18S.Rds")

trainingSet <- readRDS("/SAN/Susanas_den/AmpMarkers/wildEimeria18S/EimrefTrainingSet.RDS")

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

MA_dada2Sp <- dada2::assignSpecies(seqs, "/SAN/Susanas_den/AmpMarkers/wildEimeria18S/Eim_ref.fa")
SA_dada2Sp <- dada2::assignSpecies(seqs2, "/SAN/Susanas_den/AmpMarkers/wildEimeria18S/Eim_ref.fa", allowMultiple=TRUE)

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
#rownames(MA_dada2Sp) <- NULL
#rownames(SA_dada2Sp) <- NULL
#rownames(MA_dada2Tx) <- NULL
#rownames(SA_dada2Tx) <- NULL

MA_dada2Sp
SA_dada2Sp

MA_dada2Tx
SA_dada2Tx

sapply(SA_IDtax, function(x) {
    x$taxon[8]})

sapply(MA_IDtax, function(x) {
    x$taxon[8]})

SA_IDtax

MA_IDtax

#quick inspection here
rownames(MA_dada2Sp)==rownames(Eim@tax_table)

Eim@tax_table[,7] <- MA_dada2Sp[,2]
Eimp <- psmelt(Eim)

Eimeria_reads <- ggplot(Eimp, aes(x=Sample, y=Abundance, fill=species))+
    geom_bar(stat="identity", position="stack", color="black", size=0.02)+
    theme_classic()

Eimeria_reads2

ggplot2::ggsave(file="fig/MA_Eimeria_reads2.pdf", Eimeria_reads2, width = 5, height = 3, dpi = 300)

rownames(SA_dada2Sp)==rownames(Eim2@tax_table)
Eim2@tax_table[,7] <- SA_dada2Sp[,2]
Eimp2 <- psmelt(Eim2)

Eimeria_reads2 <- ggplot(Eimp2, aes(x=Sample, y=Abundance, fill=species))+
    geom_bar(stat="identity", position="stack", color="black", size=0.02)+
    theme_classic()

ggplot2::ggsave(file="fig/SA_Eimeria_reads.pdf", Eimeria_reads2, width = 5, height = 5, dpi = 300)

##aligments
#first better name reads
names(seqs) <- MA_dada2Sp[,2]
names(seqs) <- paste(names(seqs),c("MA1", "MA2", "MA3"))
names(seqs2) <- SA_dada2Sp[,2]
names(seqs2) <- paste(names(seqs2), c("SA1", "SA2", "SA3", "SA4", "SA5"))

refEim <- readDNAStringSet("/SAN/Susanas_den/AmpMarkers/wildEimeria18S/Eim_ref.fa")

poolSeqs <- c(seqs, seqs2)
alignment <- AlignSeqs(poolSeqs, anchor=NA, verbose=FALSE)
alignmentdb <- AlignSeqs(refEim, anchor=NA, verbose=FALSE)

# we need an outgroup

Teur <- c(readDNAStringSet("/SAN/Susanas_den/AmpMarkers/wildEimeria18S/Teur1.fasta"),
          readDNAStringSet("/SAN/Susanas_den/AmpMarkers/wildEimeria18S/Teur2.fasta"),
          readDNAStringSet("/SAN/Susanas_den/AmpMarkers/wildEimeria18S/Teur3.fasta"))

names(Teur) <- gsub("(.*sp.)(.*)","\\1", names(Teur))

allal <- c(poolSeqs, refEim, Teur)
Allal <- AlignSeqs(allal, anchor=NA)

library(ShortRead)

writeFasta(allal, "/SAN/Susanas_den/AmpMarkers/wildEimeria18S/EimreadsRef.fa")

writeFasta(Allal, "/SAN/Susanas_den/AmpMarkers/wildEimeria18S/EimreadsRef_allign.fa")

#library(msa)

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


############
Pall <- phyDat(as(Allal, "matrix"), type="DNA")
dm <- dist.ml(Pall)
treeUPGMA  <- upgma(dm)
treeNJ <- NJ(dm)

plot(treeNJ)

plot(treeUPGMA)

fun <- function(x) upgma(dist.ml(x))
bs_upgma <- bootstrap.phyDat(Pall, fun)
plotBS(treeUPGMA, bs_upgma, main="UPGMA")

fit <- pml(treeNJ, data=Pall)

fit2 <- pml(treeUPGMA, data=Pall)

fitJC  <- optim.pml(fit, rearrangement="NNI")

fitJC2  <- optim.pml(fit2, rearrangement="NNI")

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,optNni=TRUE,
                    optBf=TRUE, optQ=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))


fitGTR2 <- update(fit2, k=4, inv=0.2)
fitGTR2 <- optim.pml(fitGTR2, model="GTR", optInv=TRUE, optGamma=TRUE,optNni=TRUE,
                    optBf=TRUE, optQ=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))


#anova(fitJC, fitGTR)
#bs1 <- bootstrap.pml(fitJC, bs=100, optNni=TRUE,
#                    control = pml.control(trace = 0))
#plotBS(midpoint(fitJC$tree), bs1, p = 50, type="p")
#bs <- bootstrap.pml(fitGTR, bs=100, optNni=TRUE,
#                    control = pml.control(trace = 0))
#plotBS(midpoint(fitGTR$tree), bs, p = 50, type="p")

library(ggtree)

#pat <- c("Eimeria ferrisi", "Eimeria vermiformis", "Eimeria falciformis", "SA", "MA")

unloadNamespace("microbiome")
unloadNamespace("phyloseq")

group

group <- fitGTR$tree$tip.label

ferrisi <- group[grepl("Eimeria ferrisi", fitGTR$tree$tip.label)]

falciformis <- group[grepl("Eimeria falciformis", fitGTR$tree$tip.label)]

vermiformis <- group[grepl("Eimeria vermiformis", fitGTR$tree$tip.label)]

Isospora <- group[grepl("Isospora", fitGTR$tree$tip.label)]

#MA <- group[grepl("MA", fitGTR$tree$tip.label)]
#SA <- group[grepl("SA", fitGTR$tree$tip.label)]


cls <- list(ferrisi, falciformis, vermiformis, Isospora)


mytree <- fitGTR$tree

mytree2 <- fitGTR2$tree

ggtree(mytree2)+geom_text2(aes(subset=!isTip, label=node), hjust=-.3)

ggtree(mytree)+geom_text2(aes(subset=!isTip, label=node), hjust=-.3)

###

tree1 <- groupOTU(mytree, cls)

tree2 <- groupOTU(mytree2, cls)

tree1

library("colorspace")

t1 <- ggtree(tree1, aes(color=group)) +
    scale_color_manual(values=c("black", rainbow_hcl(4))) + theme(legend.position="right")+
    geom_nodepoint(aes(color=group), size=3, alpha=.8) 


t1

t2 <- ggtree(tree2, aes(color=group)) +
    scale_color_manual(values=c("black", rainbow_hcl(4))) + theme(legend.position="right")+
        geom_nodepoint(aes(color=group), size=3, alpha=.8) 

t2

ggplot2::ggsave(file="fig/Eimeria_NJtree.pdf", t1, width = 5, height = 5, dpi = 300)
ggplot2::ggsave(file="fig/Eimeria_NJtree.png", t1, width = 10, height = 10, dpi = 300)
ggplot2::ggsave(file="fig/Eimeria_Utree.pdf", t2, width = 5, height = 5, dpi = 300)


t1