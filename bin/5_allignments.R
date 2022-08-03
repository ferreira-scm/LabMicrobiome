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

#library(BiocManager, lib.loc="/usr/local/lib/R/site-library/")

#remotes::install_github("vmikk/metagMisc")
## using the devel
devtools::load_all("/SAN/Susanas_den/MultiAmplicon/")

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
names(seqs) <- MA_dada2Sp[,2]
names(seqs) <- paste(names(seqs),c("MA1", "MA2", "MA3"),sep="")
names(seqs2) <- SA_dada2Sp[,2]
names(seqs2) <- paste(names(seqs2), c("SA1", "SA2", "SA3", "SA4", "SA5"), sep="")


refEim <- readDNAStringSet("/SAN/Susanas_den/AmpMarkers/wildEimeria18S/Eim_ref.fa")

names(refEim) <- gsub("(\\s)", "_", names(refEim))

poolSeqs <- c(seqs, seqs2)

allSeqs <- c(poolSeqs, refEim)

alignment <- AlignSeqs(poolSeqs, anchor=NA, verbose=FALSE)

alignmentdb <- AlignSeqs(refEim, anchor=NA, verbose=FALSE)


# we need an outgroup

Teur <- c(readDNAStringSet("/SAN/Susanas_den/AmpMarkers/wildEimeria18S/Teur1.fasta"),
          readDNAStringSet("/SAN/Susanas_den/AmpMarkers/wildEimeria18S/Teur2.fasta"),
          readDNAStringSet("/SAN/Susanas_den/AmpMarkers/wildEimeria18S/Teur3.fasta"))

names(Teur) <- gsub("(.*sp.)(.*)","\\1", names(Teur))

alignmentAll <- AlignSeqs(c(poolSeqs, refEim, Teur), anchor=NA)

writeFasta(alignmentAll, "tmp/Eimeria_alignment.fa")

writeFasta(alignment, "tmp/Eimeria_reads_alignment.fa")

#library(msa)

#BrowseSeqs(alignment, highlight=0)

max(DistanceMatrix(alignment))

DistanceMatrix(alignment)

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

SA.e$ASV <- "ASV"
SA.e$ASV[which(SA.e$OTU==colnames(Eim2@otu_table)[5])] <- "ASV5"
SA.e$ASV[which(SA.e$OTU==colnames(Eim2@otu_table)[4])] <- "ASV4"
SA.e$ASV[which(SA.e$OTU==colnames(Eim2@otu_table)[3])] <- "ASV3"
SA.e$ASV[which(SA.e$OTU==colnames(Eim2@otu_table)[2])] <- "ASV2"
SA.e$ASV[which(SA.e$OTU==colnames(Eim2@otu_table)[1])] <- "ASV1"

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


lrtest(gau.m)

lrtest(gau.m, gau.m1)
lrtest(gau.m, gau.m2)
lrtest(gau.m, gau.m3)
lrtest(gau.m, gau.m4)
lrtest(gau.m, gau.m5)

anova(gau.m, test="LRT")

summary(gau.m)

plot(gau.m) # not great

summary(SA.eim.m) # terrible

library(lmtest)
lrtest(gau.m, gau.m0)

library(relaimpo)

calc.relimp(gau.m)

MA.gau <- lm(log(1+Eim@sam_data$Genome_copies_gFaeces)~log(1+Eim@otu_table[,1])+log(1+Eim@otu_table[,2])+log(1+Eim@otu_table[,3]))

summary(MA.gau)

plot(MA.gau) # uff, not great either!

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
SA5 <- ggplot(SA.e5, aes(x=dpi, y=Abundance, fill=EH_ID))+
    geom_point(shape=21, position=position_jitter(0.2), size=4, alpha=0.7)+
    geom_line(aes(group=EH_ID), alpha=0.2)+
    scale_fill_manual(values=coul)

keep <- SA.e4$EH_ID[SA.e4$Abundance>0]
SA.e4 <- SA.e4[which(SA.e4$EH_ID%in%keep),]
SA4 <- ggplot(SA.e4, aes(x=dpi, y=Abundance, fill=EH_ID))+
    geom_point(shape=21, position=position_jitter(0.2), size=4, alpha=0.7)+
    geom_line(aes(group=EH_ID), alpha=0.2)+
    scale_fill_manual(values=coul)

keep <- SA.e3$EH_ID[SA.e3$Abundance>0]
SA.e3 <- SA.e3[which(SA.e3$EH_ID%in%keep),]
SA3 <- ggplot(SA.e3, aes(x=dpi, y=Abundance, fill=EH_ID))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, alpha=0.7)+
    geom_line(aes(group=EH_ID), alpha=0.2)+
    scale_fill_manual(values=coul)

keep <- SA.e2$EH_ID[SA.e2$Abundance>0]
SA.e2 <- SA.e2[which(SA.e2$EH_ID%in%keep),]
SA2 <- ggplot(SA.e2, aes(x=dpi, y=Abundance, fill=EH_ID))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, alpha=0.7)+
    geom_line(aes(group=EH_ID, alpha=0.2))+
    scale_fill_manual(values=coul)

keep <- SA.e1$EH_ID[SA.e1$Abundance>0]
SA.e1 <- SA.e1[which(SA.e1$EH_ID%in%keep),]
SA1 <- ggplot(SA.e1, aes(x=dpi, y=Abundance, fill=EH_ID))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, alpha=0.7)+
    geom_line(aes(group=EH_ID), alpha=0.2)+
    scale_fill_manual(values=coul)

### ASV from MA run plotted
keep <- MA.e3$EH_ID[MA.e3$Abundance>0]
MA.e3 <- MA.e3[which(MA.e3$EH_ID%in%keep),]
MA3 <- ggplot(MA.e3, aes(x=dpi, y=Abundance, fill=EH_ID))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, alpha=0.7)+
    geom_line(aes(group=EH_ID), alpha=0.2)+
    scale_fill_manual(values=coul)

keep <- MA.e2$EH_ID[MA.e2$Abundance>0]
MA.e2 <- MA.e2[which(MA.e2$EH_ID%in%keep),]
MA2 <- ggplot(MA.e2, aes(x=dpi, y=Abundance, fill=EH_ID))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, alpha=0.7)+
    geom_line(aes(group=EH_ID, alpha=0.2))+
    scale_fill_manual(values=coul)

keep <- MA.e1$EH_ID[MA.e1$Abundance>0]
MA.e1 <- MA.e1[which(MA.e1$EH_ID%in%keep),]
MA1 <- ggplot(MA.e1, aes(x=dpi, y=Abundance, fill=EH_ID))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, alpha=0.7)+
    geom_line(aes(group=EH_ID), alpha=0.2)+
    scale_fill_manual(values=coul)


plot_grid(SA1,SA2,SA3,SA4,SA5) -> SA.asv

plot_grid(MA1,MA2,MA3) -> MA.asv

SA.asv

ggplot2::ggsave(file="fig/SA/Eimeria_ASVs_dpi.pdf", SA.asv, width = 10, height = 5, dpi = 300)

ggplot2::ggsave(file="fig/MA/Eimeria_ASVs_dpi.pdf", MA.asv, width = 10, height = 5, dpi = 300)

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

