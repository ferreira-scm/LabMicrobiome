#!/usr/bin/Rscript

library(ggplot2)
library(dada2)
#library(MultiAmplicon, lib.loc="/usr/local/lib/R/site-library/")
library(reshape)
library(phyloseq)
library(DECIPHER)
library(phangorn)
library(ShortRead)
library(RColorBrewer)
library("seqinr")
library(treeio)
library(ggtree)
library(magrittr)

#library(BiocManager, lib.loc="/usr/local/lib/R/site-library/")
## All runs pooled
all.PS <- readRDS("tmp/PhyloSeqData_All.Rds")
all.PS.l <- readRDS("tmp/PhyloSeqList_All.Rds")
# single amplicon from MA run
all.PSwang <- all.PS.l[[37]]

## Single amplicon 18S
sin.PS18S <- readRDS("tmp/PS_18Swang.Rds")
sin.PS18S.slv <- readRDS("tmp/PS_18Swang_SILVA.Rds")

## Single amplicon "pooled"
#sin.PS <- readRDS("tmp/PhyloSeqData18S.Rds")
#sin.PS.slv <- readRDS("tmp/PhyloSeqData18S_SILVA.Rds")

trainingSet <- readRDS("/SAN/Susanas_den/AmpMarkers/wildEimeria18S/EimrefTrainingSet.RDS")

source("bin/PlottingCor.R")
# let's filter
f.sin18 <- fil(sin.PS18S)
f.all <- fil(all.PS)
f.allwang <- fil(all.PSwang)
f.sin18.slv <- fil(sin.PS18S.slv)
# and transform
T.sin18 <- f.sin18
T.all <- f.all
T.allwang <-f.allwang
T.sin18.slv <- f.sin18.slv
otu_table(T.sin18) <- otu_table(T.sin18)*sample_data(T.sin18)$DNA_g_feces
otu_table(T.all) <- otu_table(T.all)*sample_data(T.all)$DNA_g_feces
otu_table(T.allwang) <- otu_table(T.allwang)*sample_data(T.allwang)$DNA_g_feces
otu_table(T.sin18.slv) <- otu_table(T.sin18.slv)*sample_data(T.sin18.slv)$DNA_g_feces

#all.TSS <- transform_sample_counts(f.all, function(x) x / sum(x))

## OK, now we want all the Eimeria sequences
Eim <- subset_taxa(T.all, family%in%"Eimeriidae")
Eim2 <- subset_taxa(T.sin18, family%in%"Eimeriidae")
Eim.slv <- subset_taxa(T.sin18.slv, Family%in%"Eimeriorina")
Eim.w <- subset_taxa(T.allwang, family%in%"Eimeriidae")

nfEim <- subset_taxa(all.PS, family%in%"Eimeriidae")
nfEim2 <- subset_taxa(sin.PS18S, family%in%"Eimeriidae")
nfEim.w <- subset_taxa(all.PSwang, family%in%"Eimeriidae")

get_taxa_unique(Eim, "genus")

# Do ASV match beween MA and SA?
rownames(Eim2@tax_table) %in% rownames(Eim@tax_table)

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

rownames(Eim.slv@tax_table)==rownames(Eim2@tax_table)

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


Eim2@tax_table[,7] <- paste("sa.read", seq(1:nrow(Eim2@tax_table)), sep="")

Eim2@tax_table[1,7] <- SA_dada2Sp[1,2]

Eimp.sa <- psmelt(Eim2)

Eimp.sa.slv <- psmelt(Eim.slv)

plot_bar(Eim2, "species")

Eimeria_reads <- ggplot(Eimp, aes(x=Sample, y=Abundance, fill=species))+
    geom_bar(stat="identity", position="stack", color="black", size=0.02)+
    theme_classic()

Eimeria_reads.sa <- ggplot(Eimp.sa, aes(x=Sample, y=Abundance, fill=species))+
    geom_bar(stat="identity", position="stack", color="black", size=0.02)+
    theme_classic()

Eimeria_reads.sa

Eimeria_reads.sa.slv <- ggplot(Eimp.sa.slv, aes(x=Sample, y=Abundance, fill=Genus))+
    geom_bar(stat="identity", position="stack", color="black", size=0.02)+
    theme_classic()

Eimeria_reads.sa.2 <- ggplot(Eimp.sa, aes(x=Sample, y=Abundance, fill=OTU))+
    geom_bar(stat="identity", position="stack", color="black", size=0.02)+
    theme_classic()

ggplot2::ggsave(file="fig/SA_Eimeria_reads.sa.pdf", Eimeria_reads.sa, width = 5, height = 3, dpi = 300)
ggplot2::ggsave(file="fig/SA_Eimeria_reads.sa.png", Eimeria_reads.sa, width = 5, height = 3, dpi = 300)

##aligments
#first better name reads
names(seqs2) <- c("ASV1", "ASV2", "ASV3", "ASV4", "ASV5")

## let's get all the sequences used in Jarquín-Díaz et al. 2019


access.l <- c("JQ993669", "JQ993670","JQ993671","JQ993665", "JQ993659", "KU174470", "KU174461", "KU174464", "KU174480", "KU174483", "KU174462", "KU174463", "KU174468", "KU174479", "KU174449", "KU174450", "KU174465", "KU174467", "KU174474", "KU174451", "KU174459", "KU174485", "KU174487", "KU174472", "JQ993666","JQ993649", "JQ993650", "AF080614", "MH751998", "KT360995", "JF304148", "U40263","JQ993651", "AF311643","AF311644", "JQ993654", "JQ993655", "JQ993656", "JQ993657", "JQ993658", "JQ993660", "KU174475", "JQ993661", "JQ993662", "JQ993663", "JQ993664", "JQ993652", "JQ993667", "KU174454", "KU174469", "KU174481", "KU174456","KU174484", "KU174478", "KU174473", "KU174471",  "KU174455", "KU174457", "KU174476","KU174466", "KU174486", "KU174453","KU174458", "KU174460", "KU174452","KU174482", "AF246717", "KT184355","JQ993653", "AF307880", "AF339489","AF307878", "AF307876", "AF324214","AF339490", "AF339491", "AF307879", "AF339492", "AF307877", "AF311642", "KU192965", "KU192958", "KU192936", "KU192961", "KU192956", "KU192931", "KU192938", "KU192916")
              
eim.db <- read.GenBank(access.l)

eim.db.ID <- paste(names(eim.db), attr(eim.db, "species"), sep="_")

eim.db.ID

# convert to DNAStringset

eim.DB <- eim.db %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
names(eim.DB) <- eim.db.ID

refEim <- readDNAStringSet("/SAN/Susanas_den/AmpMarkers/wildEimeria18S/Eim_ref.fa")
names(refEim) <- gsub("(\\s)", "_", names(refEim))
names(refEim)

allSeqs <- c(eim.DB,refEim, seqs2)

alignment <- AlignSeqs(seqs2, anchor=NA, verbose=FALSE, iterations=10, refinements=10, processors=90)
Allal <- AlignSeqs(allSeqs, anchor=NA, verbose=FALSE, iterations=10, refinements=10, processors=90)
Allal <- AdjustAlignment(Allal)

writeFasta(allSeqs, "tmp/Eimeria_seqs.fa")
writeFasta(Allal, "tmp/Eimeria_alignment.fa")
writeFasta(alignment, "tmp/Eimeria_reads_alignment.fa")

#BrowseSeqs(Allal, highlight=0)

# read our tree constructed with iqtree2
#~/iqtree-2.2.0-Linux/bin/iqtree2 -s tmp/Eimeria_alignment.fa -m TN+F+R2 -B 5000 -redo

my_tree <- read.iqtree('tmp/Eimeria_alignment.fa.iqtree') error



Peim <- phyDat(as(Allal, "matrix"), type="DNA")
dm <- dist.ml(Peim)
treeUPGMA  <- upgma(dm)

fun <- function(x) upgma(dist.ml(x))
bs_upgma <- bootstrap.phyDat(Peim, fun)
plotBS(treeUPGMA, bs_upgma, main="UPGMA")

# for NJ too
fun <- function(x) NJ(dist.ml(x))
bs_NJ <- bootstrap.phyDat(Peim, fun)
plotBS(treeNJ, bs_NJ, main="NJ")

### ML

fit <- pml(treeNJ, data=Peim)
#fit2 <- pml(treeUPGMA, data=Peim)
fitJC  <- optim.pml(fit, rearrangement="NNI")
#fitJC2  <- optim.pml(fit2, rearrangement="NNI")
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,optNni=TRUE, optBf=TRUE, optQ=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))

#fitGTR2 <- update(fit2, k=4, inv=0.2)
#fitGTR2 <- optim.pml(fitGTR2, model="GTR", optInv=TRUE, optGamma=TRUE,optNni=TRUE, optBf=TRUE, optQ=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))


anova(fitJC, fitGTR)
#bs1 <- bootstrap.pml(fitJC, bs=100, optNni=TRUE,
#                    control = pml.control(trace = 0))
#plotBS(midpoint(fitJC$tree), bs1, p = 50, type="p")
#bs <- bootstrap.pml(fitGTR, bs=100, optNni=TRUE,
#                    control = pml.control(trace = 0))
#plotBS(midpoint(fitGTR$tree), bs, p = 50, type="p")

library(ggtree)

#pat <- c("Eimeria ferrisi", "Eimeria vermiformis", "Eimeria falciformis", "SA", "MA")

#unloadNamespace("microbiome")
#unloadNamespace("phyloseq")

group <- fitGTR$tree$tip.label

ferrisi <- group[grepl("Eimeria ferrisi", fitGTR$tree$tip.label)]

falciformis <- group[grepl("Eimeria falciformis", fitGTR$tree$tip.label)]

vermiformis <- group[grepl("Eimeria vermiformis", fitGTR$tree$tip.label)]

Isospora <- group[grepl("Isospora", fitGTR$tree$tip.label)]

#MA <- group[grepl("MA", fitGTR$tree$tip.label)]
#SA <- group[grepl("SA", fitGTR$tree$tip.label)]


cls <- list(ferrisi, falciformis, vermiformis, Isospora)


mytree <- fitGTR$tree

#mytree2 <- fitGTR2$tree

ggtree(mytree)+geom_tiplab()

#+geom_text2(aes(subset=!isTip, label=node), hjust=-.3)

#ggtree(mytree)+geom_text2(aes(subset=!isTip, label=node), hjust=-.3)

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
#####################################################
########################### Eimeria ASVs
## Are ASV's the same between SA and MA?
colnames(Eim@otu_table)[1] == colnames(Eim2@otu_table)[1]
colnames(Eim@otu_table)[2] == colnames(Eim2@otu_table)[2]
colnames(Eim@otu_table)[3] == colnames(Eim2@otu_table)[3]

#plotting individual ASV's
MA.e <- psmelt(Eim)
SA.e <- psmelt(Eim2)

Eim2.g <- tax_glom(Eim2, taxrank="family")
SA.e.g <- psmelt(Eim2.g)
Eim.g <- tax_glom(Eim, taxrank="family")
MA.e.g <- psmelt(Eim.g)
Eim.g.w <- tax_glom(Eim.w, taxrank="family")
MA.e.g.w <- psmelt(Eim.g.w)

nfEim2.g <- tax_glom(nfEim2, taxrank="family")
nfSA.e.g <- psmelt(nfEim2.g)
nfEim.g <- tax_glom(nfEim, taxrank="family")
nfMA.e.g <- psmelt(nfEim.g)
nfEim.g.w <- tax_glom(nfEim.w, taxrank="family")
nfMA.e.g.w <- psmelt(nfEim.g.w)

SA.e$ASV <- "ASV"
SA.e$ASV[which(SA.e$OTU==colnames(Eim2@otu_table)[5])] <- "ASV5"
SA.e$ASV[which(SA.e$OTU==colnames(Eim2@otu_table)[4])] <- "ASV4"
SA.e$ASV[which(SA.e$OTU==colnames(Eim2@otu_table)[3])] <- "ASV3"
SA.e$ASV[which(SA.e$OTU==colnames(Eim2@otu_table)[2])] <- "ASV2"
SA.e$ASV[which(SA.e$OTU==colnames(Eim2@otu_table)[1])] <- "ASV1"

MA.e$ASV <- "ASV"
MA.e$ASV[which(MA.e$OTU==colnames(Eim@otu_table)[3])] <- "ASV3"
MA.e$ASV[which(MA.e$OTU==colnames(Eim@otu_table)[2])] <- "ASV2"
MA.e$ASV[which(MA.e$OTU==colnames(Eim@otu_table)[1])] <- "ASV1"


SA.e5 <- SA.e[which(SA.e$OTU==colnames(Eim2@otu_table)[5]),]
SA.e4 <- SA.e[which(SA.e$OTU==colnames(Eim2@otu_table)[4]),]
SA.e3 <- SA.e[which(SA.e$OTU==colnames(Eim2@otu_table)[3]),]
SA.e2 <- SA.e[which(SA.e$OTU==colnames(Eim2@otu_table)[2]),]
SA.e1 <- SA.e[which(SA.e$OTU==colnames(Eim2@otu_table)[1]),]

MA.e3 <- MA.e[which(MA.e$OTU==colnames(Eim@otu_table)[3]),]
MA.e2 <- MA.e[which(MA.e$OTU==colnames(Eim@otu_table)[2]),]
MA.e1 <- MA.e[which(MA.e$OTU==colnames(Eim@otu_table)[1]),]

nb.col=length(levels(as.factor(SA.e5$EH_ID)))+1
coul <- colorRampPalette(brewer.pal(8, "Accent"))(nb.col)

library(MASS)

SA.eim.m <- glm.nb(Eim2@sam_data$Genome_copies_gFaeces~Eim2@otu_table[,1]+Eim2@otu_table[,2]+Eim2@otu_table[,3]+Eim2@otu_table[,4]+Eim2@otu_table[,5])

gau.m <- lm(log(1+Eim2@sam_data$Genome_copies_gFaeces)~log(1+Eim2@otu_table[,1])+log(1+Eim2@otu_table[,2])+log(1+Eim2@otu_table[,3])+log(Eim2@otu_table[,4]+1)+log(1+Eim2@otu_table[,5]))

gau.m1 <- lm(log(1+Eim2@sam_data$Genome_copies_gFaeces)~log(1+Eim2@otu_table[,2])+log(1+Eim2@otu_table[,3])+log(Eim2@otu_table[,4]+1)+log(1+Eim2@otu_table[,5]))
gau.m2 <- lm(log(1+Eim2@sam_data$Genome_copies_gFaeces)~log(1+Eim2@otu_table[,1])+log(1+Eim2@otu_table[,3])+log(Eim2@otu_table[,4]+1)+log(1+Eim2@otu_table[,5]))
gau.m3 <- lm(log(1+Eim2@sam_data$Genome_copies_gFaeces)~log(1+Eim2@otu_table[,1])+log(1+Eim2@otu_table[,2])+log(Eim2@otu_table[,4]+1)+log(1+Eim2@otu_table[,5]))
gau.m4 <- lm(log(1+Eim2@sam_data$Genome_copies_gFaeces)~log(1+Eim2@otu_table[,1])+log(1+Eim2@otu_table[,2])+log(1+Eim2@otu_table[,3])+log(1+Eim2@otu_table[,5]))
gau.m5 <- lm(log(1+Eim2@sam_data$Genome_copies_gFaeces)~log(1+Eim2@otu_table[,1])+log(1+Eim2@otu_table[,2])+log(1+Eim2@otu_table[,3])+log(Eim2@otu_table[,4]+1))
gau.m0 <- lm(log(1+Eim2@sam_data$Genome_copies_gFaeces)~1)

anova(gau.m, test="LRT")
summary(gau.m)
#plot(gau.m) # not great
#summary(SA.eim.m) # terrible

### testing Victor's hypothesis for oocysts
gau.o <- lm(log(1+Eim2@sam_data$OPG)~log(1+Eim2@otu_table[,1])+log(1+Eim2@otu_table[,2])+log(1+Eim2@otu_table[,3])+log(Eim2@otu_table[,4]+1)+log(1+Eim2@otu_table[,5]))

cor.test(log(1+SA.e2$Abundance), log(1+SA.e2$OPG))
cor.test(log(1+SA.e2$Abundance), log(1+SA.e2$Genome_copies_gFaeces))

#ggplot(SA.e1, aes(log(1+Abundance), log(1+OPG)))+
#    geom_jitter()

#ggplot(SA.e1, aes(log(1+Abundance), log(1+Genome_copies_gFaeces)))+
#    geom_jitter()

summary(gau.o)


MA.gau <- lm(log(1+Eim@sam_data$Genome_copies_gFaeces)~log(1+Eim@otu_table[,1])+log(1+Eim@otu_table[,2])+log(1+Eim@otu_table[,3]))
summary(MA.gau)

calc.relimp(MA.gau)

anova(MA.gau, test="LRT")

#plot(MA.gau) # uff, not great either!

SA_Eimeiria.ASVs <- ggplot(SA.e, aes(x=log(1+Genome_copies_gFaeces), y=log(1+Abundance), fill=ASV))+
    geom_point(shape=21, size=4, alpha=0.7)+
    scale_fill_manual(values=c("#009E73", "#F0E442", "#0072B2",
                               "#D55E00", "#E63C9A"))

ggplot2::ggsave(file="fig/SA/SA_EimeriaASVs_qPCR.pdf", SA_Eimeiria.ASVs, width = 5, height = 3, dpi = 300)

SA_Eimeiria.ASVs

### make glm qPCR~ASV1+ASV2+...
## ASV5+ASV5 useful?
#### Plot by EH_ID, only positive animals
keep <- SA.e5$EH_ID[SA.e5$Abundance>0]
SA.e5 <- SA.e5[which(SA.e5$EH_ID%in%keep),]
SA5 <- ggplot(SA.e5, aes(x=dpi, y=log(1+Abundance), fill=EH_ID))+
    geom_point(shape=21, position=position_jitter(0.2), size=2, alpha=0.7)+
    geom_line(aes(group=EH_ID), alpha=0.2)+
    scale_fill_manual(values=coul)+
    labs(y="ASV5 abundance, log(+1)", x="Days post infection")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          legend.position = "none",
          axis.line = element_line(colour = "black"))


keep <- SA.e4$EH_ID[SA.e4$Abundance>0]
SA.e4 <- SA.e4[which(SA.e4$EH_ID%in%keep),]
SA4 <- ggplot(SA.e4, aes(x=dpi, y=log(1+Abundance), fill=EH_ID))+
    geom_point(shape=21, position=position_jitter(0.2), size=2, alpha=0.7)+
    geom_line(aes(group=EH_ID), alpha=0.2)+
    labs(y="ASV4 abundance, log(+1)", x="Days post infection")+
    scale_fill_manual(values=coul)+    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          legend.position = "none",
          axis.line = element_line(colour = "black"))


keep <- SA.e3$EH_ID[SA.e3$Abundance>0]
SA.e3 <- SA.e3[which(SA.e3$EH_ID%in%keep),]
SA3 <- ggplot(SA.e3, aes(x=dpi, y=log(1+Abundance), fill=EH_ID))+
    geom_point(shape=21, position=position_jitter(0.2), size=2, alpha=0.7)+
    geom_line(aes(group=EH_ID), alpha=0.2)+
    scale_fill_manual(values=coul)+
    labs(y="ASV3 abundance, log(+1)", x="Days post infection")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          legend.position = "none",
          axis.line = element_line(colour = "black"))


keep <- SA.e2$EH_ID[SA.e2$Abundance>0]
SA.e2 <- SA.e2[which(SA.e2$EH_ID%in%keep),]
SA2 <- ggplot(SA.e2, aes(x=dpi, y=log(1+Abundance), fill=EH_ID))+
    geom_point(shape=21, position=position_jitter(0.2), size=2, alpha=0.7)+
    geom_line(aes(group=EH_ID), alpha=0.2)+
    scale_fill_manual(values=coul)+
        labs(y="ASV2 abundance, log(+1)", x="Days post infection")+
        theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          legend.position = "none",
          axis.line = element_line(colour = "black"))


keep <- SA.e1$EH_ID[SA.e1$Abundance>0]
SA.e1 <- SA.e1[which(SA.e1$EH_ID%in%keep),]
SA1 <- ggplot(SA.e1, aes(x=dpi, y=log(1+Abundance), fill=EH_ID))+
    geom_point(shape=21, position=position_jitter(0.2), size=2, alpha=0.7)+
    geom_line(aes(group=EH_ID), alpha=0.2)+
    scale_fill_manual(values=coul)+
    labs(y="ASV1 abundance, log(+1)", x="Days post infection")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          legend.position = "none",
          axis.line = element_line(colour = "black"))


### ASV from MA run plotted
keep <- MA.e3$EH_ID[MA.e3$Abundance>0]
MA.e3 <- MA.e3[which(MA.e3$EH_ID%in%keep),]
MA3 <- ggplot(MA.e3, aes(x=dpi, y=log(1+Abundance), fill=EH_ID))+
    geom_point(shape=21, position=position_jitter(0.2), size=2, alpha=0.7)+
    geom_line(aes(group=EH_ID), alpha=0.2)+
    scale_fill_manual(values=coul)+
    labs(y="ASV3 abundance, log(+1)", x="Days post infection")+
        theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          legend.position = "none",
          axis.line = element_line(colour = "black"))


keep <- MA.e2$EH_ID[MA.e2$Abundance>0]
MA.e2 <- MA.e2[which(MA.e2$EH_ID%in%keep),]
MA2 <- ggplot(MA.e2, aes(x=dpi, y=log(1+Abundance), fill=EH_ID))+
    geom_point(shape=21, position=position_jitter(0.2), size=2, alpha=0.7)+
    geom_line(aes(group=EH_ID), alpha=0.2)+
    scale_fill_manual(values=coul)+
    labs(y="ASV2 abundance, log(+1)", x="Days post infection")+
        theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          legend.position = "none",
          axis.line = element_line(colour = "black"))


keep <- MA.e1$EH_ID[MA.e1$Abundance>0]
MA.e1 <- MA.e1[which(MA.e1$EH_ID%in%keep),]
MA1 <- ggplot(MA.e1, aes(x=dpi, y=log(1+Abundance), fill=EH_ID))+
    geom_point(shape=21, position=position_jitter(0.2), size=2, alpha=0.7)+
    geom_line(aes(group=EH_ID), alpha=0.2)+
    scale_fill_manual(values=coul)+
    labs(y="ASV1 abundance, log(+1)", x="Days post infection")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          legend.position = "none",
          axis.line = element_line(colour = "black"))


qpcr <- ggplot(MA.e1, aes(x=dpi, y=log(1+Genome_copies_gFaeces), fill=EH_ID))+
    geom_point(shape=21, position=position_jitter(0.2), size=2, alpha=0.7)+
    geom_line(aes(group=EH_ID), alpha=0.2)+
    scale_fill_manual(values=coul)+
    labs(y="Genome copies log(1+)", x="Days post infection")+
        theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          legend.position = "none",
          axis.line = element_line(colour = "black"))

ma.sa <- SA.e[, c("ASV", "Abundance", "EH_ID", "dpi", "labels")]
ma.sa.t <- MA.e[, c("ASV", "Abundance", "labels")]
ma.sa <- merge(ma.sa, ma.sa.t, by=c("ASV", "labels"))
# remove zeros

t <- ma.sa[ma.sa$Abundance.x>0,]
t <- t[t$Abundance.y>0,]

cor.test(log(1+t$Abundance.x), log(1+t$Abundance.y), method="pearson")

ASV.c <- ggplot(ma.sa, aes(x=log(1+Abundance.x), y=log(1+Abundance.y), fill=ASV))+
    geom_point(size=2, shape=21, alpha=0.5)+
    scale_fill_manual(values=c("#009E73", "#F0E442", "#0072B2"))+
    xlab("SA - ASV abundance (log1+)")+
    ylab("MA - ASV abundance (log1+)")+
    annotate(geom="text", x=15, y=19, label="Pearson rho=0.92, p<0.001", size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=11),
          legend.position = "top",
          axis.line = element_line(colour = "black"))

ASV.c

plot_grid(SA1,SA2,SA3,SA4,SA5) -> SA.asv
plot_grid(MA1,MA2,MA3, nrow=2) -> MA.asv
plot_grid(SA1,MA1,SA2,MA2,SA3,MA3,qpcr, ASV.c,  nrow=4, labels="auto") -> SAMA.asv

ggplot2::ggsave(file="fig/Eimeria_ASVs_dpi.pdf", SAMA.asv, width = 10, height = 15, dpi = 300)
ggplot2::ggsave(file="fig/Eimeria_ASVs_dpi.png", SAMA.asv, width = 10, height = 15, dpi = 300)
ggplot2::ggsave(file="fig/SA/Eimeria_ASVs_dpi.pdf", SA.asv, width = 15, height = 10, dpi = 300)
ggplot2::ggsave(file="fig/SA/Eimeria_ASVs_dpi.png", SA.asv, width = 15, height = 10, dpi = 300)
ggplot2::ggsave(file="fig/MA/Eimeria_ASVs_dpi.pdf", MA.asv, width = 10, height = 8, dpi = 300)
ggplot2::ggsave(file="fig/MA/Eimeria_ASVs_dpi.png", MA.asv, width = 10, height = 8, dpi = 300)

ggplot2::ggsave(file="fig/MA_SA_Eimeria_ASVs.pdf", ASV.c, width = 5, height = 5, dpi = 300)
ggplot2::ggsave(file="fig/MA_SA_Eimeria_ASVs.png", ASV.c, width = 5, height = 5, dpi = 300)


#################### plotting individuals by ASV
library(cowplot) # to plot a list of plots

head(SA.e5$EH_ID)

cl <- colorRampPalette(brewer.pal(8, "Accent"))(5)

length(levels(SA.e$EH_ID))

SA.e$EH_ID

p.ID <- function(i){
    SA.e%>%
    dplyr::filter(EH_ID%in%SA.e$EH_ID[i])%>%
    dplyr::select(EH_ID, dpi, ASV, Genome_copies_gFaeces, Abundance)%>%
    ggplot(aes(x=dpi, Abundance+1, fill=ASV))+
    geom_point(shape=21, position=position_jitter(0.2), size=4, alpha=0.7, aes(fill=ASV), color="black")+
    geom_line(aes(group=ASV), color="gray", alpha=0.5)+
    scale_fill_manual(values=c("#009E73", "#F0E442", "#0072B2", "#D55E00", "#E63C9A"))+
    scale_y_log10("log10 (Eimeira /gFaeces + 1) (qPCR)",
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    theme_bw()+
    theme(text = element_text(size=16), axis.title.x = element_blank(), legend.position = "top") -> nm
    nm
}

p.EH <- list()
for (i in 1:22){
    p <- p.ID(i)
    p.EH[[i]] <- p
}

pdf("fig/SA/Eimeria_ASVs_ID.pdf")
for (i in 1:22) {
    print(p.EH[[i]])
}
dev.off()

p.ID <- function(i){
    SA.e%>%
    dplyr::filter(EH_ID%in%SA.e$EH_ID[i])%>%
    dplyr::select(EH_ID, dpi, ASV, Genome_copies_gFaeces, Abundance)%>%
    ggplot(aes(x=dpi, Abundance+1, fill=ASV))+
    geom_point(shape=21, position=position_jitter(0.2), size=4, alpha=0.7, aes(fill=ASV), color="black")+
    geom_line(aes(group=ASV), color="gray", alpha=0.5)+
    scale_fill_manual(values=c("#009E73", "#F0E442", "#0072B2", "#D55E00", "#E63C9A"))+
    scale_y_log10("log10 (Eimeira /gFaeces + 1) (qPCR)",
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    theme_bw()+
    theme(text = element_text(size=16), axis.title.x = element_blank(), legend.position = "top") -> nm
    nm
}

p.EH <- list()
for (i in 1:22){
    p <- p.ID(i)
    p.EH[[i]] <- p
}

pdf("fig/SA/Eimeria_ASVs_ID.pdf")
for (i in 1:22) {
    print(p.EH[[i]])
}
dev.off()

head(MA.e)

p.ID.MA <- function(i){
    MA.e%>%
    dplyr::filter(EH_ID%in%MA.e$EH_ID[i])%>%
    dplyr::select(EH_ID, dpi, ASV, Genome_copies_gFaeces, Abundance)%>%
    ggplot(aes(x=dpi, Abundance+1, fill=ASV))+
    geom_point(shape=21, position=position_jitter(0.2), size=4, alpha=0.7, aes(fill=ASV), color="black")+
    geom_line(aes(group=ASV), color="gray", alpha=0.5)+
    scale_fill_manual(values=c("#009E73", "#F0E442", "#0072B2"))+
    scale_y_log10("log10 (Eimeira /gFaeces + 1) (qPCR)",
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    theme_bw()+
    theme(text = element_text(size=16), axis.title.x = element_blank(), legend.position = "top") -> nm
    nm
}

p.EH.M <- list()
for (i in 1:22){
    p <- p.ID.MA(i)
    p.EH.M[[i]] <- p
}

p.EH.M[1]

pdf("fig/MA/Eimeria_ASVs_ID.pdf")
for (i in 1:22) {
    print(p.EH.M[[i]])
}
dev.off()

################################



ggplot(SA.e, aes(x=log(1+Abundance), y=ASV))+
    geom_point(aes(fill=ASV), size=4, shape=21, apha=0.5, position=position_jitter(height=0.2))+
#      geom_density(aes(color = ASV), size = 1) +
    scale_fill_manual(values=c("#009E73", "#F0E442", "#0072B2", "#D55E00", "#E63C9A"))+
#        scale_colour_manual(values=c("#009E73", "#F0E442", "#0072B2", "#D55E00", "#E63C9A"))+
#    ylab("SA - ASV abundance (log1+)")+
#    xlab("MA - ASV abundance (log1+)")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
           legend.position = "top",
          axis.line = element_line(colour = "black"))


ASV1 <- ma.sa[ma.sa$ASV=="ASV1",]
ASV2 <- ma.sa[ma.sa$ASV=="ASV2",]
ASV3 <- ma.sa[ma.sa$ASV=="ASV3",]

cor.test(ASV1$Abundance.x, ASV1$Abundance.y)
cor.test(ASV2$Abundance.x, ASV2$Abundance.y)
cor.test(log(1+ASV1$Abundance.x), log(1+ASV1$Abundance.y))
cor.test(log(1+ASV2$Abundance.x), log(1+ASV2$Abundance.y))
cor.test(log(1+ASV3$Abundance.x), log(1+ASV3$Abundance.y))

ASV1.m <- ggplot(ASV1, aes(x=log(1+Abundance.x), y=log(1+Abundance.y)))+
    geom_point(size=4, shape=21, alpha=0.5, fill="#009E73")+
    ylab("SA - ASV3 abundance (log1+)")+
    xlab("MA - ASV3 abundance (log1+)")+
        annotate(geom="text", x=5, y=18, label="Pearson rho=0.78, p<0.001")+
        theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
           legend.position = "none",
          axis.line = element_line(colour = "black"))

ASV2.m <- ggplot(ASV2, aes(x=log(1+Abundance.x), y=log(1+Abundance.y)))+
    geom_point(size=4, shape=21, alpha=0.5, fill= "#F0E442")+
    ylab("SA - ASV3 abundance (log1+)")+
    xlab("MA - ASV3 abundance (log1+)")+
        annotate(geom="text", x=5, y=16, label="Pearson rho=0.67, p<0.001")+
        theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
           legend.position = "none",
          axis.line = element_line(colour = "black"))


ASV3.m <- ggplot(ASV3, aes(x=log(1+Abundance.x), y=log(1+Abundance.y)))+
    geom_point(size=4, shape=21, alpha=0.5, fill="#0072B2")+
    ylab("SA - ASV3 abundance (log1+)")+
    xlab("MA - ASV3 abundance (log1+)")+
        annotate(geom="text", x=5, y=13, label="Pearson rho=0.57, p<0.001")+
        theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
           legend.position = "none",
          axis.line = element_line(colour = "black"))

ASV3.m

bottom_row <- plot_grid(ASV1.m,ASV2.m,ASV3.m, labels=c("b", "c", "d"), nrow=1, label_size=12)

ASV.cor <- plot_grid(ASV.c,bottom_row, labels=c("a", ""), label_size=12, ncol=1)


ASV.cor

ggplot2::ggsave(file="fig/ASV_CORRELATION.pdf", ASV.cor, width = 10, height = 10, dpi = 600)



############# sensitivity and specificity
#remove NA row
SA.e.g <- SA.e.g[-which(is.na(SA.e.g$Genome_copies_gFaeces)),]
nrow(SA.e.g)

sensit <- function(SA.e.g){
SA.eim <- matrix(as.numeric(c(summary(SA.e.g$Abundance>0 & SA.e.g$Genome_copies_gFaeces>0)[3],
    summary(SA.e.g$Abundance==0 & SA.e.g$Genome_copies_gFaeces>0)[3],
    summary(SA.e.g$Abundance>0 & SA.e.g$Genome_copies_gFaeces==0)[3],
    summary(SA.e.g$Abundance==0 & SA.e.g$Genome_copies_gFaeces==0)[3])), ncol=2, byrow=TRUE)
margin1 <- margin.table(SA.eim, margin=1)
margin2 <- margin.table(SA.eim, margin=2)
eimN <- margin1[2]
eimP <- margin1[1]
testP <- margin2[1]
testN <- margin2[2]
tN <- margin2[2]
tP <- margin2[1]
tP <- SA.eim[1,1]
fP <- SA.eim[2,1]
tN <- SA.eim[2,2]
fN <- SA.eim[1,2]
print("Sensitivity:")
print(tP/eimP)
print("Specificity")
print(tN/eimN)
print("positive predictive value")
print(tP/testP)
print("negative predictive value")
print(tN/testN)
}

#sensitivity for the 3 datasets filtered and not filtered
sensit(SA.e.g)

sensit(nfSA.e.g)

sensit(MA.e.g)
sensit(nfMA.e.g)

sensit(MA.e.g.w)
sensit(nfMA.e.g.w)

tP/testP

testP
