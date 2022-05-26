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
PSwang <- readRDS("tmp/PS_wang.Rds")

sam <- data.frame(sample_data(PS))

#let's remove samples with no reads
PS <- prune_samples(sample_sums(PS)>0, PS)
PS18S <- prune_samples(sample_sums(PS18S)>0, PS18S)
PSwang <- prune_samples(sample_sums(PSwang)>0, PSwang)

### here we start with the single amplicon analysis (wang)
#No filters
PSeimf <- subset_taxa(PSwang, family%in%"Eimeriidae")

#create total sums and Eimeria sums data frame
df <- data.frame(sample_sums(otu_table(PSwang)))
df$labels <- rownames(df)
eimf <-as.data.frame(sample_sums(PSeimf))
eimf$labels <- rownames(eimf)
names(eimf) <- c("EimeriaSums", "labels")
names(df) <- c("TotalSums", "labels")

#merge
df <- merge(df,eimf, by="labels", all=FALSE) 
sdt <- merge(df,sam, by="labels")

#correlation tests
cor.test(sdt$Genome_copies_gFaeces, sdt$EimeriaSums, method="spearman")
cor.test(sdt$OPG, sdt$EimeriaSums, method="spearman")
cor.test(sdt$TotalSums, sdt$EimeriaSums, method="spearman")

# Linear models
Alm <- lm(Genome_copies_gFaeces ~ EimeriaSums + TotalSums, sdt)
summary(Alm)

# plotting no filter correlation
a <-ggplot(sdt, aes(log(1+Genome_copies_gFaeces), log(1+EimeriaSums)))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
        xlab("Genome copies gFaeces(log 1+)")+
    ylab("SA Eimeriidae (log1+)")+
    ggtitle("No filter, raw counts")+
        labs(tag= "a)")+
        annotate(geom="text", x=12, y=7, label="Spearman rho=0.92, p<0.001")+
        theme_bw()+
    theme(text = element_text(size=16))

## Now we filter for abundance 0.01% and prevalence 1%
#PS <- prune_samples(sample_sums(PS)>0, PS)
#PS.l <- lapply(PS.l, function(x){

# abundance filtering to 0.001%?
x = taxa_sums(PSwang)
keepTaxa = (x / sum(x) > 0.00001)
summary(keepTaxa)
pPS = prune_taxa(keepTaxa, PSwang)

# plus prevalnce filter at 1%
ppPS <- phyloseq_filter_prevalence(pPS, prev.trh=0.01)

KeepTaxap <- prevalence(pPS)>0.01
ppPS <- prune_taxa(KeepTaxap, pPS)
ppPS <- prune_samples(sample_sums(ppPS)>0, ppPS)
# subset samples based on total read count (500 reads)
ppPS <- phyloseq::subset_samples(ppPS, phyloseq::sample_sums(PS) > 500)
ppPS <- prune_samples(sample_sums(ppPS)>0, ppPS)
ppPS

#now we make the data frame
bPSeimf <- subset_taxa(ppPS, family%in%"Eimeriidae")

#create total sums and Eimeria sums data frame

bdf <- data.frame(sample_sums(otu_table(ppPS)))
bdf$labels <- rownames(bdf)
bdf$eim <-(sample_sums(otu_table(bPSeimf)))
names(bdf) <- c("Fil_TotalSums", "labels", "FilEimeriaSums")

#merge
sdt <- merge(bdf,sdt, by="labels", all=TRUE) 

nrow(sdt)

#correlation tests
cor.test(sdt$Genome_copies_gFaeces, sdt$FilEimeriaSums, method="spearman")
cor.test(sdt$FilEimeriaSums, sdt$Fil_TotalSums, method="spearman")
cor.test(sdt$OPG, sdt$FilEimeriaSums, method="spearman")

# Linear models
Blm <- lm(Genome_copies_gFaeces ~ Fil_EimeriaSums, bsdt)

# plotting with abundance and prevalence filter correlation
b <-ggplot(sdt, aes(log(1+Genome_copies_gFaeces), log(1+FilEimeriaSums)))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
        xlab("Genome copies gFaeces(log 1+)")+
    ylab("multiamplicon Eimeriidae (log1+)")+
    ggtitle("prev 1%, ab 0.001%, 500 sample")+
        labs(tag= "b)")+
        annotate(geom="text", x=12, y=7, label="Spearman rho=0.93, p<0.001")+
        theme_bw()+
    theme(text = element_text(size=16))

#### using relative abundance
PSTSS = transform_sample_counts(ppPS, function(x) x / sum(x))
cPSeimf <- subset_taxa(PSTSS, family%in%"Eimeriidae")

#create total sums and Eimeria sums data frame
bdf <- data.frame(sample_sums(otu_table(ppPS)))
bdf$labels <- rownames(bdf)
bdf$eim <-(sample_sums(otu_table(bPSeimf)))
names(bdf) <- c("Fil_TotalSums", "labels", "FilEimeriaSums")


df <-data.frame(sample_sums(otu_table(cPSeimf)))
df$labels <- rownames(df)
names(df) <- c("TSS_Eim", "labels")

#merge
sdt <- merge(df, sdt, by="labels", all=TRUE) 

#correlation tests
cor.test(sdt$Genome_copies_gFaeces, sdt$TSS_Eim, method="spearman")
cor.test(sdt$OPG, sdt$TSS_Eim, method="spearman")

# Linear models
Clm <- lm(Genome_copies_gFaeces ~ TSS_Eims, sdt)
summary(Clm)

#library(lmtest)
#coxtest(Alm, Blm)
#jtest(Alm, Blm)

# plotting TSS correlation
c <-ggplot(sdt, aes(log(1+Genome_copies_gFaeces), log(1+TSS_Eim)))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
        xlab("Genome copies gFaeces log(1+)")+
    ylab("SA Eimeriidae log(1+)")+
    ggtitle("TSS")+
        labs(tag= "c)")+
        annotate(geom="text", x=13, y=0.5, label="Spearman rho=0.94, p<0.001")+
        theme_bw()+
    theme(text = element_text(size=16))

c

Sys.setenv("DISPLAY"=":10.0")

#### using Relative log expression
library(edgeR)
edgePS <- phyloseq_to_edgeR(bPSeimf)
edgePS$samples$labels <- rownames(edgePS$samples)
df <- (edgePS$samples)

head(df)

df$group <- NULL
df$lib.size <- NULL

names(df) <- c("REL_Eim", "labels")

#merge
sdt <- merge(df, sdt, by="labels", all=TRUE) 

#hist(log10(apply(otu_table(edgePS),1,var)),
#     breaks=87,
#     xlab="log10(variance)")

#correlation tests
cor.test(sdt$Genome_copies_gFaeces, sdt$REL_Eim, method="spearman")
cor.test(sdt$OPG, sdt$REL_Eim, method="spearman")

# Linear models
Dlm <- lm(Genome_copies_gFaeces ~ REL_Eim, sdt)
summary(Dlm)

# plotting REL correlation
d <-ggplot(sdt, aes(log(1+Genome_copies_gFaeces), REL_Eim))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
        xlab("Genome copies gFaeces log(1+)")+
    ylab("SA Eimeriidae")+
    ggtitle("REL")+
        labs(tag= "c)")+
        annotate(geom="text", x=20, y=5, label="Spearman rho=-0.90, p<0.001")+
        theme_bw()+
    theme(text = element_text(size=16))

d

# plotting VST correlation
# I think this is not correct. Come back to this later
DPS <- bPSeimf
otu_table(DPS) <- otu_table(DPS)+1
deseq <- phyloseq_to_deseq2(DPS, ~Strain)
deseq <- estimateSizeFactors(deseq)
deseq <- estimateDispersions(deseq)
deseq.vst <- getVarianceStabilizedData(deseq)
deseq.vst <- data.frame(deseq.vst)
deseq <- t(deseq.vst)
deseq.vst$vst_Eim <- rowSums(deseq.vst)

# save plots of what we have so far
plot_grid(a,b,c,d) -> fCor

ggplot2::ggsave(file="fig/SACOrrs.pdf", fCor, width = 10, height = 10, dpi = 600)

fCor

png(filename="fig/SA_cor.png",
    width =10, height = 10, units = "in", res= 300)
fCor
dev.off()


########### we start with our "biological" normalization

# what is the most abundant taxa?

#remotes::install_github("gmteunisse/Fantaxtic")

library("fantaxtic")

top <- get_top_taxa(ppPS, 10, relative=TRUE, discard_other=FALSE, other_label="Other")

ps_tmp <- name_taxa(top, label="Unknown", species = T, other_label="Other")

fantaxtic_bar(ps_tmp, color_by="Phylum", label_by="family", other_label="Other")

get_taxa_unique(top, "family")

most_abundant_taxa <- sort(taxa_sums(ppPS), TRUE)[1:20]

PS20 <- prune_taxa(names(most_abundant_taxa), ppPS)

PS20

get_taxa_unique(PS20, "family")


library(MASS)

