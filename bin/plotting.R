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
    geom_bar(stat="identity", position="stack", color="gray")+
    scale_fill_manual(values=mycolors)+
#    scale_color_brewer(palette="Set2")+
    theme_minimal()

PHYdpi


sample_data(Tps18S)$dpi <- factor(sample_data(Tps18S)$dpi)

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
    geom_bar(stat="identity", position="stack", color="gray")+
    scale_fill_manual(values=mycolors)+
#    scale_color_brewer(palette="Set2")+
    theme_minimal()

PHYdpi18

# save

ggplot2::ggsave(file="fig/phylum_abundance_DPI.pdf", PHYdpi, width = 5, height = 5, dpi = 300)
ggplot2::ggsave(file="fig/phylum_abundance_DPI.png", PHYdpi, width = 5, height = 5, dpi = 300)

ggplot2::ggsave(file="fig/phylum_abundance_DPI18S.pdf", PHYdpi18, width = 5, height = 5, dpi = 300)
ggplot2::ggsave(file="fig/phylum_abundance_DPI18S.png", PHYdpi18, width = 5, height = 5, dpi = 300)


names(prevalencedf)

sample_data(Tps18S)$dpi <- as.factor(sample_data(Tps18S)$dpi)

library("fantaxtic")

top <- get_top_taxa(fPS, 10, relative=TRUE, discard_other=FALSE, other_label="Other")
ps_tmp <- name_taxa(top, label="Unknown", species = T, other_label="Other")
all <- fantaxtic_bar(ps_tmp, color_by="phylum", label_by="family", other_label="Other")

all

top <- get_top_taxa(fPS18S, 10, relative=FALSE, discard_other=FALSE, other_label="Other")
ps_tmp <- name_taxa(top, label="Unknown", species = T, other_label="Other")
wang <- fantaxtic_bar(ps_tmp, color_by="family", label_by="genus", other_label="Other")

png(filename="fig/abundance_all.png",
    width =7, height = 5, units = "in", res= 300)
all
dev.off()
png(filename="fig/abundance_wang.png",
    width =7, height = 5, units = "in", res= 300)
wang
dev.off()

##########################################
Tps18S

plot_bar(Tps18S, fill="phylum", x="dpi")

plot_bar(Tps, fill="phylum", x="dpi")


#plot_heatmap(Tps, sample.label="dpi", sample.order="dpi")

#plot_heatmap(Tps18S, sample.label="dpi", sample.order="dpi", taxa.label="family")

