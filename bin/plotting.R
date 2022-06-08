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

source("bin/4_MA_Correlations.R")
source("bin/PlottingCor.R")


ps.ord <- ordinate(Tps, "NMDS", "bray")

p1 <- plot_ordination(Tps, ps.ord, type="samples", color="dpi", shape="HybridStatus")
print(p1)

p3 <- plot_ordination(Tps, ps.ord, type="taxa", color="phylum")
print(p3)

plot(sample_data(Tps)$DNA_g_feces ~ as.numeric(sample_data(Tps)$dpi))

ps18S.ord <- ordinate(Tps18S, "NMDS", "bray")
p2 <- plot_ordination(Tps18S, ps18S.ord, type="samples", color="dpi")
print(p2)


dpiPS <- merge_samples(Tps, "dpi")

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

#quick dirty fix
sample_data(Tps18S) <- sample_data(Tps)

dpi18 <- merge_samples(Tps18S, "dpi")

phyDPI18 <- aggregate_taxa(dpi18, level = "phylum")

dpi.melt18 <- psmelt(phyDPI18)

dpi.melt18$phylum <- factor(dpi.melt18$phylum)

nb.c <-length(levels(dpi.melt18$phylum))

mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.c)

PHYdpi18 <- ggplot(dpi.melt18, aes(x=dpi, y=Abundance, fill=fct_reorder(phylum, Abundance)))+
#    coord_flip()+
#    geom_jitter(shape=21, alpha = 0.7, position=position_jitter(0.2), size=4, aes(fill=phylum))+
    geom_bar(stat="identity", position="stack", color="black", size=0.02)+
    scale_fill_manual(values=mycolors)+
    labs(fill="Phylum")+
#    scale_color_brewer(palette="Set2")+
    theme_minimal()

PHYdpi18

# save

ggplot2::ggsave(file="fig/phylum_abundance_DPI.pdf", PHYdpi, width = 5, height = 5, dpi = 300)
ggplot2::ggsave(file="fig/phylum_abundance_DPI.png", PHYdpi, width = 5, height = 5, dpi = 300)

ggplot2::ggsave(file="fig/phylum_abundance_DPI18S.pdf", PHYdpi18, width = 5, height = 5, dpi = 300)
ggplot2::ggsave(file="fig/phylum_abundance_DPI18S.png", PHYdpi18, width = 5, height = 5, dpi = 300)


##### intra individual abundance

sample_data(Tps)$EH_ID <- as.factor(sample_data(Tps)$EH_ID)
sample_data(Tps18S)$EH_ID <- as.factor(sample_data(Tps18S)$EH_ID)

ID <- merge_samples(Tps, "EH_ID")
ID18 <- merge_samples(Tps18S, "EH_ID")

phyID <- aggregate_taxa(ID, level = "phylum")
phyID18 <- aggregate_taxa(ID18, level = "phylum")

IDmelt18 <- psmelt(phyID18)
IDmelt <- psmelt(phyID)

IDmelt18$phylum <- factor(IDmelt18$phylum)
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


nb.c <-length(levels(IDmelt18$phylum))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.c)
ID18 <- ggplot(IDmelt18, aes(x=EH_ID, y=Abundance, fill=fct_reorder(phylum, Abundance)))+
#    coord_flip()+
#    geom_jitter(shape=21, alpha = 0.7, position=position_jitter(0.2), size=4, aes(fill=phylum))+
    geom_bar(stat="identity", position="stack", color="black", size=0.02)+
    scale_fill_manual(values=mycolors)+
    labs(fill="Phylum")+
#    scale_color_brewer(palette="Set2")+
    theme_minimal()

ID18

pTps <- aggregate_taxa(Tps, level = "phylum")
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


pTps18 <- aggregate_taxa(Tps18S, level = "phylum")
pTps18.mel <- psmelt(pTps18)

pTps18.mel$phylum <- as.factor(pTps18.mel$phylum)
nb.c <-length(levels(pTps18.mel$phylum))
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.c)

abps18 <- ggplot(pTps18.mel, aes(x=labels, y=Abundance, fill=fct_reorder(phylum, Abundance)))+
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


abps18

ggplot2::ggsave(file="fig/phylum_abundance_18S.pdf", abps18, width = 7, height = 5, dpi = 300)
ggplot2::ggsave(file="fig/phylum_abundance_all.pdf", abps, width = 7, height = 5, dpi = 300)

ggplot2::ggsave(file="fig/phylum_abundance_18S.png", abps18, width = 7, height = 5, dpi = 300)
ggplot2::ggsave(file="fig/phylum_abundance_all.png", abps, width = 7, height = 5, dpi = 300)
