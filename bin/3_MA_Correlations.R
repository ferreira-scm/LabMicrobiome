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

## using the devel
#devtools::load <- all("/SAN/Susanas_den/MultiAmplicon/")

PS <- readRDS("tmp/PhyloSeqData_All.Rds")
PS.l <- readRDS("tmp/PhyloSeqList_All.Rds")

#quick look of what's in there
get_taxa_unique(PS, "family")
sum(otu_table(subset_taxa(PS, superkingdom%in%"Eukaryota")))/sum(otu_table(PS))
sum(otu_table(subset_taxa(PS, family%in%"Eimeriidae")))/sum(otu_table(PS))

sam <- data.frame(sample_data(PS))

#let's remove samples with no reads
PS <- prune_samples(sample_sums(PS)>0, PS)

PS

#No filters
PSeimf <- subset_taxa(PS, family%in%"Eimeriidae")

#create total sums and Eimeria sums data frame
df <- data.frame(sample_sums(otu_table(PS)))
df$labels <- rownames(df)
eimf <-as.data.frame(sample_sums(PSeimf))
eimf$labels <- rownames(eimf)
names(eimf) <- c("EimeriaSums", "labels")
names(df) <- c("TotalSums", "labels")

#merge
df <- merge(df,eimf, by="labels", all=FALSE) 
sdt <- merge(df,sam, by="labels") 

#correlation tests
cor.test(sdt$Genome_copies_gFaeces, sdt$EimeriaSums)
cor.test(sdt$Genome_copies_gFaeces, sdt$TotalSums)
cor.test(sdt$OPG, sdt$EimeriaSums)

# plotting no filter correlation
a <-ggplot(sdt, aes(log(1+Genome_copies_gFaeces), log(1+EimeriaSums)))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
        xlab("Genome copies gFaeces(log 1+)")+
    ylab("multiamplicon Eimeriidae (log1+)")+
    ggtitle("No filter, raw counts")+
        labs(tag= "a)")+
        annotate(geom="text", x=5, y=2.5, label="Spearman rho=0.53, p<0.001")+
        theme_bw()+
    theme(text = element_text(size=16))


a

## Now we filter for abundance 0.01% and prevalence 1%
#PS <- prune_samples(sample_sums(PS)>0, PS)
#PS.l <- lapply(PS.l, function(x){

# abundance filtering to 0.001%?
x = taxa_sums(PS)
keepTaxa = (x / sum(x) > 0.00001)
summary(keepTaxa)
pPS = prune_taxa(keepTaxa, PS)

# plus prevalnce filter at 1%
ppPS <- phyloseq_filter_prevalence(pPS, prev.trh=0.01)

# subset samples based on total read count (500 reads)
summary(phyloseq::sample_sums(ppPS))

ppPS <- prune_samples(sample_sums(ppPS)>0, ppPS)

ppPS <- phyloseq::subset_samples(ppPS, phyloseq::sample_sums(PS) > 500)

ppPS

#now we make the data frame
bPSeimf <- subset_taxa(ppPS, family%in%"Eimeriidae")

#create total sums and Eimeria sums data frame
bdf <- data.frame(sample_sums(otu_table(ppPS)))
bdf$labels <- rownames(bdf)
beimf <-as.data.frame(sample_sums(bPSeimf))
beimf$labels <- rownames(beimf)
names(beimf) <- c("EimeriaSums", "labels")
names(bdf) <- c("TotalSums", "labels")

#merge
bdf <- merge(bdf,beimf, by="labels", all=FALSE) 
bsdt <- merge(bdf,sam, by="labels") 

#correlation tests
cor.test(bsdt$Genome_copies_gFaeces, bsdt$EimeriaSums)
cor.test(bsdt$EimeriaSums          , bsdt$TotalSums)
cor.test(bsdt$OPG, bsdt$EimeriaSums)


# plotting with abundance and prevalence filter correlation
b <-ggplot(bsdt, aes(log(1+Genome_copies_gFaeces), log(1+EimeriaSums)))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
        xlab("Genome copies gFaeces(log 1+)")+
    ylab("multiamplicon Eimeriidae (log1+)")+
    ggtitle("prev 1%, ab 0.001%, 500 sample")+
        labs(tag= "b)")+
        annotate(geom="text", x=5, y=2.5, label="Spearman rho=0.53, p<0.001")+
        theme_bw()+
    theme(text = element_text(size=16))


b

### OK now only positive samples in MA

Csdt <- sdt[sdt$EimeriaSums>0,]

#correlation tests
cor.test(Csdt$Genome_copies_gFaeces, Csdt$EimeriaSums)
cor.test(Csdt$EimeriaSums, Csdt$TotalSums)
cor.test(Csdt$OPG, Csdt$EimeriaSums)


# plotting with abundance and prevalence filter correlation
c <-ggplot(Csdt, aes(log(1+Genome_copies_gFaeces), log(1+EimeriaSums)))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    xlab("Genome copies gFaeces(log 1+)")+
    ylab("multiamplicon Eimeriidae (log1+)")+
    ggtitle("no filter, raw, only positive")+
        labs(tag= "c)")+
        annotate(geom="text", x=5, y=2.5, label="Spearman rho=0.56, p<0.001")+
        theme_bw()+
    theme(text = element_text(size=16))


c

#### using relative abundance

PSTSS = transform_sample_counts(PS, function(x) x / sum(x))

#No filters
dPSeimf <- subset_taxa(PSTSS, family%in%"Eimeriidae")

#create total sums and Eimeria sums data frame
ddf <- data.frame(sample_sums(otu_table(PSTSS)))
ddf$labels <- rownames(ddf)
deimf <-as.data.frame(sample_sums(dPSeimf))
deimf$labels <- rownames(deimf)
names(deimf) <- c("EimeriaSums", "labels")
names(ddf) <- c("TotalSums", "labels")

#merge
ddf <- merge(ddf,deimf, by="labels", all=FALSE) 
dsdt <- merge(ddf,sam, by="labels") 

#correlation tests
cor.test(dsdt$Genome_copies_gFaeces, dsdt$EimeriaSums)
cor.test(dsdt$Genome_copies_gFaeces, dsdt$TotalSums)
cor.test(dsdt$OPG, dsdt$EimeriaSums)

# plotting no filter correlation
d <-ggplot(dsdt, aes(Genome_copies_gFaeces, EimeriaSums))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
        xlab("Genome copies gFaeces")+
    ylab("multiamplicon Eimeriidae")+
    ggtitle("Relative abundance, no filter")+
        labs(tag= "d)")+
        annotate(geom="text", x=10, y=0.5, label="Spearman rho=0.52, p<0.001")+
        theme_bw()+
    theme(text = element_text(size=16))


d

plot_grid(a,b,c,d) -> fCor

ggplot2::ggsave(file="fig/FirstCOrrs.pdf", fCor, width = 10, height = 10, dpi = 600)
