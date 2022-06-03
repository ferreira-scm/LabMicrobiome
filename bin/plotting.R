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

#plot_heatmap(Tps, sample.label="dpi", sample.order="dpi")

#plot_heatmap(Tps18S, sample.label="dpi", sample.order="dpi", taxa.label="family")

ps.ord <- ordinate(Tps, "NMDS", "bray")

p1 <- plot_ordination(Tps, ps.ord, type="samples", color="dpi", shape="HybridStatus")
print(p1)

p3 <- plot_ordination(Tps, ps.ord, type="taxa", color="phylum")
print(p3)


plot(sample_data(Tps)$DNA_g_feces ~ as.numeric(sample_data(Tps)$dpi))

ps18S.ord <- ordinate(Tps18S, "NMDS", "bray")
p2 <- plot_ordination(Tps18S, ps18S.ord, type="samples", color="dpi")
print(p2)



sample_data(Tps18S)$dpi <- as.factor(sample_data(Tps18S)$dpi)

library("fantaxtic")

top <- get_top_taxa(fPS, 10, relative=TRUE, discard_other=FALSE, other_label="Other")
ps_tmp <- name_taxa(top, label="Unknown", species = T, other_label="Other")
all <- fantaxtic_bar(ps_tmp, color_by="phylum", label_by="family", other_label="Other")

all

top <- get_top_taxa(fPS18S, 10, relative=TRUE, discard_other=FALSE, other_label="Other")
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

