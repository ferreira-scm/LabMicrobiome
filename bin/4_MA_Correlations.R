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

source("bin/PlottingCor.R")

PS <- readRDS("tmp/PhyloSeqData_All.Rds")
PS.l <- readRDS("tmp/PhyloSeqList_All.Rds")
PS18S <- readRDS("tmp/PS_18Swang.Rds")
PSwang <- PS.l[[37]]


# let's filter
fPS18S <-fil(PS18S)
fPS <- fil(PS)

# and transform
Tps18S <- fPS18S
Tps <- fPS
otu_table(Tps18S) <- otu_table(fPS18S)*sample_data(fPS18S)$DNA_g_feces
otu_table(Tps) <- otu_table(fPS)*sample_data(fPS)$DNA_g_feces

Plotting_cor(ps=PS, "MA", dir="fig/MA/")
Plotting_cor(ps=PS18S, "SA", dir="fig/SA/")

#Plotting_cor(ps=PSwang, name="wang1141_13_F.Nem_0425_6_3_R", dir="fig/MA/")

for (i in 1:48) {
    nm <- names(PS.l)[i]
    ps <- PS.l[[i]]
    try(Plotting_cor(ps, name=nm, dir="fig/MA/"))
}

for (i in 1:48) {
    nm <- names(PS.l)[i]
    ps <- PS.l[[i]]
    try(NoFilPlotting_cor(PS=ps, name=nm, dir="fig/MA/NoFil/"))
}

for (i in 1:48) {
#    print(names(PS.l)[i])
    try(p <- subset_taxa(PS.l[[i]],phylum=="Apicomplexa"), silent=TRUE)
    try(get_taxa_unique(p, "family"), silent=TRUE)
    if (exists("p")) {
        a <- get_taxa_unique(p, "family")
        print(paste(i, "- ", names(PS.l[i]), ": ", length(a), sep=""))
        print(a)
}
    rm(p)
}

#######################
# how does host DNA change with infection
Mus <- subset_taxa(Tps, genus%in%"Mus")
Mus <-aggregate_taxa(Mus, level="genus")

# no Mus reads in SA
#Mus18S <- subset_taxa(Tps18S, genus%in%"Mus")
m <- psmelt(Mus)

mp <- ggplot(m, aes(dpi, Abundance, color=EH_ID))+
    geom_point()+
    geom_line(aes(group=EH_ID))+
    labs(x="Days post infection", y="Host reads (per g of faeces)")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none")
#          axis.line = element_line(colour = "black"))


# how does DNA g/faeces change with infection
dnap <- ggplot(m, aes(y=DNA_g_feces, x=dpi, color=EH_ID))+
    geom_point()+
    geom_line(aes(group=EH_ID))+
    labs(x="Days post infection", y="DNA (per g of faeces)")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none")


names(m)

wl <- ggplot(m, aes(y=weightloss, x=dpi, color=EH_ID))+
    geom_point()+
    geom_line(aes(group=EH_ID))+
    labs(x="Days post infection", y="DNA (per g of faeces)")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none")

ggplot2::ggsave(file="fig/Mus_abundance.pdf", mp, width = 7, height = 5, dpi = 300)
ggplot2::ggsave(file="fig/DNA_abundance_all.pdf", dnap, width = 7, height = 5, dpi = 300)
ggplot2::ggsave(file="fig/Mus_abundance.png", mp, width = 7, height = 5, dpi = 300)
ggplot2::ggsave(file="fig/DNA_abundance_all.png", dnap, width = 7, height = 5, dpi = 300)

plot_grid(mp, dnap, wl, nrow=3) -> mplot

ggplot2::ggsave(file="fig/Mus_WeighLoss.png", mplot, width = 7, height = 15, dpi = 300)

ggplot2::ggsave(file="fig/Mus_WeighLoss.pdf", mplot, width = 7, height = 15, dpi = 300)


# and does it correlate with weight loss?
cor.test(m$Abundance, m$weightloss, method="pearson")

### OK, so let's remove food

rank_names(fPS)

get_taxa_unique(PS18S, "phylum")


get_taxa_unique(fPS, "phylum")

plant18S <- subset_taxa(PS18S, !phylum%in%"Streptophyta")

plant <-  subset_taxa(PS, !phylum%in%"Streptophyta")

Plotting_cor(ps=plant18S, "SA_no_plants", dir="fig/SA/")

Plotting_cor(ps=plant, "MA_no_plants", dir="fig/MA/")


plant

PS

PlantMusWorms

PlantMusWorms <-  subset_taxa(PS, !(phylum%in%"Streptophyta"| phylum%in%"Nematoda" | phylum%in%"Chordata"))

PlantMusWorms

Plotting_cor(psb=PlantMusWorms, "MA_no_plants_Mus_Nematodes", dir="fig/MA/")



