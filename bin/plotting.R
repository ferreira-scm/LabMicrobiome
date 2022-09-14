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

remotes::install_github("microsud/microbiomeutilities")

library(microbiomeutilities) # messes up with ggplot2::alpha

## using the devel
#devtools::load_all("/SAN/Susanas_den/MultiAmplicon/")

source("bin/4_MA_SA_filtering.R")

allTSS

T.all

## going to plot with ACS
tall <- psmelt(T.all)

head(tall)
f.all.lp

phyTSS <- aggregate_top_taxa2(allTSS, level="phylum", 11)

genTSS <- aggregate_top_taxa2(allTSS, level="genus", 11)

genPS <- aggregate_top_taxa2(T.all, level="genus", 11)

gen <- psmelt(genPS)

genTSS <- psmelt(genTSS)

ggplot(gen, aes(x=Sample, y=Abundance, fill=fct_reorder(genus, Abundance)))+
    geom_bar(position="stack", stat="identity")+
    scale_x_discrete(labels=gen$EH_ID, breaks=gen$Sample)+
    scale_fill_brewer("genus", palette="Paired")+
    facet_grid(~dpi, scales="free")+
#    rremove("x.text")+
    theme_bw()

gen$Abundance

gen$Abundance

heat_gen <- ggplot(gen, aes(x=Sample, y=genus))+
    geom_tile(aes(fill=Abundance))+
    scale_fill_distiller("Abundance/ngDNA", palette = "Spectral") + theme_bw()+
    facet_grid(~dpi, scales="free")+
    theme()+
    theme(axis.text.y = element_text(colour = 'black', size = 10, face = 'italic'),
        axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.key = element_blank(),
          strip.background = element_rect(colour="black", fill="white")) 

heat_genTSS <- ggplot(genTSS, aes(x=Sample, y=genus))+
    geom_tile(aes(fill=Abundance))+
    scale_fill_distiller("Abundance", palette = "Spectral") + theme_bw()+
    facet_grid(~dpi, scales="free")+
    theme()+
    theme(axis.text.y = element_text(colour = 'black', size = 10, face = 'italic'),
        axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.key = element_blank(),
          strip.background = element_rect(colour="black", fill="white")) 

heat_genTSS

ggsave("fig/Genus_heat_composition.pdf", heat_gen, height=10, width=10, dpi=400)
ggsave("fig/Genus_heat_compositionTSS.pdf", heat_genTSS, height=10, width=10, dpi=400)

# 20 most abundant ASV's

phylum_comp <- plot_composition(phyTSS, sample.sort="dpi", x.label="dpi", otu.sort="abundance")+
    scale_fill_brewer("phylum", palette="Paired")+theme_bw()

ggsave("fig/Phylum_composition.pdf", phylum_comp, height=6, width=10, dpi=400)

gen_comp <- plot_composition(genTSS, sample.sort="dpi", x.label="dpi", otu.sort="abundance")+
    scale_fill_brewer("genus", palette="Paired")+theme_bw()
ggsave("fig/Genus_composition.pdf", gen_comp, height=6, width=10, dpi=400)

gen_comp$data

#T.all@otu_table

ggplot(tall, aes(x=reorder(labels, -as.numeric(tall$dpi)), y=Abundance, fill=phylum))+
#    coord_flip()+
#    geom_jitter(shape=21, alpha = 0.7, position=position_jitter(0.2), size=4, aes(fill=phylum))+
    geom_bar(stat="identity", position="stack", color="black", size=0.02)+
    scale_fill_manual(values=mycolors)+
    labs(fill="Phylum", x="Samples", y="Abundance (per g of faeces)")+
    theme_bw()+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
#          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
#          axis.line = element_line(colour = "black"))



## plotting ordinations
ps.ord <- ordinate(T.all, "NMDS", "bray")
p1 <- plot_ordination(T.all, ps.ord, type="samples", color="dpi")
print(p1)
p3 <- plot_ordination(allTSS, ps.ord, type="taxa", color="phylum")
print(p3)

p1

#plot(sample_data(T.all)$DNA_g_feces ~ as.numeric(sample_data(T.all)$dpi))

dpiPS <- merge_samples(T.all, "dpi")
phyDPI <- aggregate_taxa(dpiPS, level = "phylum")

phy.melt <- psmelt(phyDPI)
dpi.melt <- psmelt(dpiPS)

library(forcats)

ggplot(dpi.melt, aes(x=dpi, y=Abundance, fill=phylum))+
#    coord_flip()+
#    geom_jitter(shape=21, alpha = 0.7, position=position_jitter(0.2), size=4, aes(fill=phylum))+
    geom_bar(stat="identity", position="stack")+
    scale_color_brewer(palette="Dark2")+
    theme_minimal()

phy.melt$phylum <- factor(phy.melt$phylum)

library(RColorBrewer)

nb.c <-length(levels(phy.melt$phylum))

mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.c)

PHYdpi <- ggplot(phy.melt, aes(x=dpi, y=Abundance, fill=fct_reorder(phylum, Abundance)))+
#    coord_flip()+
#    geom_jitter(shape=21, alpha = 0.7, position=position_jitter(0.2), size=4, aes(fill=phylum))+
    geom_bar(stat="identity", position="stack", color="black", size=0.02)+
    scale_fill_manual(values=mycolors)+
    labs(fill="Phylum")+
#    scale_color_brewer(palette="Set2")+
    theme_minimal()

PHYdpi

##### intra individual abundance

sample_data(T.all)$EH_ID <- as.factor(sample_data(T.all)$EH_ID)

ID <- merge_samples(T.all, "EH_ID")

phyID <- aggregate_taxa(ID, level = "phylum")

IDmelt <- psmelt(phyID)

IDmelt$phylum <- factor(IDmelt$phylum)

nb.c <-length(levels(IDmelt$phylum))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.c)
ID <- ggplot(IDmelt, aes(x=EH_ID, y=Abundance, fill=fct_reorder(phylum, Abundance)))+
#    coord_flip()+
#    geom_jitter(shape=21, alpha = 0.7, position=position_jitter(0.2), size=4, aes(fill=phylum))+
    geom_bar(stat="identity", position="stack", color="black", size=0.02)+
    scale_fill_manual(values=mycolors)+
    labs(fill="Phylum")+
#    scale_color_brewer(palette="Set2")+
    theme_minimal()

ID

pTps <- aggregate_taxa(T.all, level = "phylum")
pTps.mel <- psmelt(pTps)

pTps.mel$phylum <- as.factor(pTps.mel$phylum)
nb.c <-length(levels(pTps.mel$phylum))
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.c)

abps <- ggplot(pTps.mel, aes(x=labels, y=Abundance, fill=fct_reorder(phylum, Abundance)))+
#    coord_flip()+
#    geom_jitter(shape=21, alpha = 0.7, position=position_jitter(0.2), size=4, aes(fill=phylum))+
    geom_bar(stat="identity", position="stack", color="black", size=0.02)+
    scale_fill_manual(values=mycolors)+
    labs(fill="Phylum", x="Samples", y="Abundance (per g of faeces)")+
    theme_bw()+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
#          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
#          axis.line = element_line(colour = "black"))

abps

