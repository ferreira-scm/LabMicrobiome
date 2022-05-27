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


# abundance filtering to 0.001%?
x = taxa_sums(PS)
keepTaxa = (x / sum(x) > 0.00001)
summary(keepTaxa)
filPS = prune_taxa(keepTaxa, PS)

# plus prevalnce filter at 1%
filPS <- phyloseq_filter_prevalence(filPS, prev.trh=0.01)

KeepTaxap <- prevalence(filPS)>0.01
filPS <- prune_taxa(KeepTaxap, filPS)
# subset samples based on total read count (500 reads)
filPS <- phyloseq::subset_samples(filPS, phyloseq::sample_sums(PS) > 500)
filPS <- prune_samples(sample_sums(filPS)>0, filPS)
filPS

library("fantaxtic")

top <- get_top_taxa(filPS, 10, relative=TRUE, discard_other=FALSE, other_label="Other")
ps_tmp <- name_taxa(top, label="Unknown", species = T, other_label="Other")
all <- fantaxtic_bar(ps_tmp, color_by="phylum", label_by="family", other_label="Other")

top <- get_top_taxa(ppPS, 10, relative=TRUE, discard_other=FALSE, other_label="Other")
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

