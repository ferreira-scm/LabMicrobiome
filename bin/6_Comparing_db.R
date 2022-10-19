#!/usr/bin/Rscript

library(ggplot2)
library(dada2)
#library(MultiAmplicon, lib.loc="/usr/local/lib/R/site-library/")
library(reshape)
library(phyloseq)

library(RColorBrewer)
library(forcats)

require(grid)
#require(grid_extra)
require(cowplot)
                                        #library(DECIPHER)
#library(phangorn)
#library(ShortRead)

#library(BiocManager, lib.loc="/usr/local/lib/R/site-library/")

                                        #remotes::install_github("vmikk/metagMisc")
## using the devel
devtools::load_all("/SAN/Susanas_den/MultiAmplicon/")

## Single amplicon 18S
#sin.PS18S <- readRDS("tmp/PS_18Swang.Rds")
#sin.PS18S.slv <- readRDS("tmp/PS_18Swang_SILVA.Rds")

## Single amplicon "pooled"
sin.PS <- readRDS("tmp/PhyloSeqData18S.Rds")
sin.PS.slv <- readRDS("tmp/PhyloSeqData18S_SILVA.Rds")

## Single amplicon "listed"
sin.PS.l <- readRDS("tmp/PhyloSeqList18S.Rds")
sin.PS.l.slv <- readRDS("tmp/PhyloSeqList18S_SILVA.Rds")

source("bin/PlottingCor.R")

# let's filter
#f.sin18 <- fil(sin.PS18S)
#f.sin18.slv <- fil(sin.PS18S.slv)

f.sin <- list()

for (i in 1:3){
    try(f.sin[[i]] <- fil(sin.PS.l[[i]]), silent=TRUE)
}


f.sin <- merge_phyloseq(f.sin[[3]], f.sin[[2]])

f.sin.l.slv <- list()
for (i in 1:3){
    try(f.sin.l.slv[[i]] <- fil(sin.PS.l.slv[[i]]), silent=TRUE)
}
f.sin.slv <- merge_phyloseq(f.sin.l.slv[[3]], f.sin.l.slv[[2]])


#get_taxa_unique(f.sin.PS.l.slv[[1]], "domain")
get_taxa_unique(f.sin.l.slv[[2]], "phylum")

get_taxa_unique(f.sin.l.slv[[3]], "domain")
get_taxa_unique(f.sin, "superkingdom")

#dirty little fix for downstream plotting
f.sin@sam_data$ave <- "average"
f.sin.slv@sam_data$ave <- "average"

rank_names(f.sin.slv)

Es <- subset_taxa(f.sin, superkingdom%in%"Eukaryota")
Bs <- subset_taxa(f.sin, superkingdom%in%"Bacteria")
#As <- subset_taxa(f.sin, superkingdom%in%"Archaea")
#ps18 <-aggregate_taxa(ps18, level="family")
#ps18 <- psmelt(ps18)
#ps18$logA <- log(1+ps18$Abundance)
#ps18$logGC <- log(1+ps18$Genome_copies_gFaeces)

Es.slv <- subset_taxa(f.sin.slv, domain%in%"d__Eukaryota")
Bs.slv <- subset_taxa(f.sin.slv, domain%in%"d__Bacteria")
As.slv <- subset_taxa(f.sin.slv, domain%in%"d__Archaea")
#ps18s <-aggregate_taxa(ps18s, level="Family")
#ps18s <- psmelt(ps18s)
#ps18s$logA <- log(1+ps18s$Abundance)
#ps18s$logGC <- log(1+ps18s$Genome_copies_gFaeces)



Esm <- psmelt(Es)
Esm.slv <- psmelt(Es.slv)
Bsm <- psmelt(Bs)
Bsm.slv <- psmelt(Bs.slv)
Asm.slv <- psmelt(As.slv)
#Esm$ave <- "average"
#Esm.slv$ave <- "average"

#nb.c <-length(levels(as.factor(Esm$phylum)))
#mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(length(levels(as.factor(Esm$phylum))))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(length(levels(as.factor(Esm.slv$phylum))))

############Eukaryotes
Esm_plot <- ggplot(Esm, aes(x=ave, y=Abundance, fill=phylum))+
    geom_bar(stat="identity", position="stack", size=0.02)+
    scale_fill_manual(values=mycolors)+
    labs(fill="Phylum")+
    theme_bw()+
    theme(text = element_text(size=16))


Esm.slv_plot <- ggplot(Esm.slv, aes(x=ave, y=Abundance, fill=phylum))+
    geom_bar(stat="identity", position="stack", size=0.02)+
    scale_fill_manual(values=mycolors)+
    labs(fill="Phylum")+
    theme_bw()+
    theme(text = element_text(size=16))

Comp_Euk <- plot_grid(Esm_plot, Esm.slv_plot)

ggplot2::ggsave(file="fig/TaxComp_Eukaryotes.png",Comp_Euk, width = 10, height = 5, dpi = 300)

#########Bacteria

mycolour <- colorRampPalette(brewer.pal(8, "Set2"))(length(levels(as.factor(Bsm.slv$phylum))))

Bs.slv_plot <- ggplot(Bsm.slv, aes(x=ave, y=Abundance, fill=phylum))+
    geom_bar(stat="identity", position="stack", size=0.02)+
    scale_fill_manual(values=mycolour)+
    labs(fill="Phylum")+
    theme_bw()+
    theme(text = element_text(size=16))

Bs_plot <- ggplot(Bsm, aes(x=ave, y=Abundance, fill=phylum))+
    geom_bar(stat="identity", position="stack", size=0.02)+
    scale_fill_manual(values=mycolour)+
    labs(fill="Phylum")+
    theme_bw()+
    theme(text = element_text(size=16))

Comp_bac <- plot_grid(Bs_plot, Bs.slv_plot, ncol=2)

Comp_bac

ggplot2::ggsave(file="fig/TaxComp_Bacteria.pdf",Comp_bac, width = 10, height = 5, dpi = 300)

As_plot <- ggplot(Asm.slv, aes(x=ave, y=Abundance, fill=phylum))+
    geom_bar(stat="identity", position="stack", size=0.02)+
    scale_fill_manual(values=mycolour)+
    labs(fill="Phylum")+
    theme_bw()+
    theme(text = element_text(size=16))

As_plot

################### for Nematodes
Nem <- subset_taxa(Es, phylum%in%"Nematoda")

nem <- psmelt(Nem)

Nem_plot_F <- ggplot(nem, aes(x=ave, y=Abundance, fill=family))+
    geom_bar(stat="identity", position="stack", size=0.02)+
    scale_fill_manual(values=mycolors)+
    labs(fill="Family")+
    theme_bw()+
    theme(text = element_text(size=16))

Nem_plot_G <-  ggplot(nem, aes(x=ave, y=Abundance, fill=genus))+
    geom_bar(stat="identity", position="stack", size=0.02)+
    scale_fill_manual(values=mycolors)+
    labs(fill="Genus")+
    theme_bw()+
    theme(text = element_text(size=16))

Nem_plot_S <-  ggplot(nem, aes(x=ave, y=Abundance, fill=species))+
    geom_bar(stat="identity", position="stack", size=0.02)+
    scale_fill_manual(values=mycolors)+
    labs(fill="Species")+
    theme_bw()+
    theme(text = element_text(size=16))

get_taxa_unique(Es.slv, "phylum")

Nem.slv <- subset_taxa(Es.slv, phylum%in%"p__Nematozoa")

nem.slv <- psmelt(Nem.slv)

Nem.slv_plot_F <- ggplot(nem.slv, aes(x=ave, y=Abundance, fill=family))+
    geom_bar(stat="identity", position="stack", size=0.02)+
    scale_fill_manual(values=mycolors)+
    labs(fill="Family")+
    theme_bw()+
    theme(text = element_text(size=16))

Nem.slv_plot_G <- ggplot(nem.slv, aes(x=ave, y=Abundance, fill=genus))+
    geom_bar(stat="identity", position="stack", size=0.02)+
    scale_fill_manual(values=mycolors)+
    labs(fill="Genus")+
    theme_bw()+
    theme(text = element_text(size=16))

Nem.slv_plot_S <- ggplot(nem.slv, aes(x=ave, y=Abundance, fill=species))+
    geom_bar(stat="identity", position="stack", size=0.02)+
    scale_fill_manual(values=mycolors)+
    labs(fill="Species")+
    theme_bw()+
    theme(text = element_text(size=16))

# save plots of what we have so far
plot_grid(Nem.slv_plot_F, Nem.slv_plot_G, Nem.slv_plot_S, ncol=3) -> NemSilva

plot_grid(Nem_plot_F, Nem_plot_G, Nem_plot_S, ncol=3) -> NemB

ggplot2::ggsave(file="fig/Tax_Comp_NemSilva.png", NemSilva, width = 12, height = 5, dpi = 600)
ggplot2::ggsave(file="fig/Tax_Comp_Nem.png", NemB, width = 12, height = 5, dpi = 600)


#### Is it a spider or is it a worm?

get_taxa_unique(Es.slv, "genus")

################# Eimeria

Eim.g <- subset_taxa(Es, genus%in%"Eimeria")
Eim.f <- subset_taxa(Es, family%in%"Eimeriidae")

Eim.s.g <- subset_taxa(Es.slv, genus%in%"g__Eimeria")
Eim.s.f <- subset_taxa(Es.slv, family%in%"f__Eimeriorina")

Eim.s.g
Eim.g

Eim.s.f
Eim.f


library(microbiome)

Eim.g <-aggregate_taxa(Eim.g, level="genus")
Eim.f <-aggregate_taxa(Eim.f, level="family")

Eim.s.g <-aggregate_taxa(Eim.s.g, level="genus")
Eim.s.f <-aggregate_taxa(Eim.s.f, level="family")

eim.g <- psmelt(Eim.g)
eim.f <- psmelt(Eim.f)
eim.s.g <- psmelt(Eim.s.g)
eim.s.f<- psmelt(Eim.s.f)

mycolors <- colorRampPalette(brewer.pal(8, "Spectral"))(11)

mycolors

library(ggpmisc)

head(eim.f)

Eim.f.p <- ggplot(eim.f, aes(y=log(1+Genome_copies_ngDNA), x=log(1+Abundance)))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    scale_fill_manual(values=mycolors)+
            geom_smooth(method = "lm", se=FALSE) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                                          parse = TRUE) +  
    ylab("Genome copies gFaeces log(1+)")+
    xlab("Blast:NCBI Eimeria log(1+)")+
    labs(fill="dpi")+
    theme_bw()+
    theme(text = element_text(size=16))

Eim.f.p


Eim.s.g.p <- ggplot(eim.s.g, aes(y=log(1+Genome_copies_ngDNA), x=log(1+Abundance)))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    scale_fill_manual(values=mycolors)+
            geom_smooth(method = "lm", se=FALSE) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                                          parse = TRUE) +  
    ylab("Genome copies gFaeces log(1+)")+
    xlab("Silva:AssignTaxonomy Eimeria log(1+)")+
    labs(fill="dpi")+
    theme_bw()+
    theme(text = element_text(size=16))


Silva_Blast_Eimeria <- plot_grid(Eim.s.g.p, Eim.f.p)

Silva_Blast_Eimeria


ggplot2::ggsave(file="fig/Silva_Blast_Eimeria.png", Silva_Blast_Eimeria, width = 8, height = 4, dpi = 600)

# this bit was moved from 4_* script it might break here on its own #################################################
################################################
#############################################
################## comparing taxonomies
Tslv <- psmelt(TPSslv)
TSA <- psmelt(TPS_SA)

library(forcats)

TPSslv

TPSslv@sam_data$glom <- "a"
TPS_SA@sam_data$glom <- "a"

Avgslv <- merge_samples(TPSslv, "glom")
Avgbla <- merge_samples(TPS_SA, "glom")

Tslv <- psmelt(Avgslv)
TSA <- psmelt(Avgbla)


silva_phy <- ggplot(Tslv, aes(dpi, Abundance, fill=fct_reorder(Phylum, Abundance)))+
    geom_bar(stat="identity", position="stack")+
    labs(y="ASV reads (per g of faeces)")+
    scale_color_brewer(palette="Dark2")+
    theme_minimal()+
        theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
#          legend.position = "none",
          axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

silva_phy

blast_phy <- ggplot(TSA, aes(dpi, Abundance, fill=fct_reorder(phylum, Abundance)))+
    geom_bar(stat="identity", position="stack")+
    labs(y="ASV reads (per g of faeces)")+
    scale_color_brewer(palette="Dark2")+
    theme_minimal()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
#          legend.position = "none",
          axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())


silva_gen <- ggplot(Tslv, aes(dpi, Abundance, fill=fct_reorderGenus, Abundance)))+
    geom_bar(stat="identity", position="stack")+
    labs(y="ASV reads (per g of faeces)")+
    scale_color_brewer(palette="Dark2")+
    theme_minimal()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

silva_gen

blast_gen <- ggplot(TSA, aes(dpi, Abundance, fill=fct_reorder(genus, Abundance)))+
    geom_bar(stat="identity", position="stack")+
    labs(y="ASV reads (per g of faeces)")+
    scale_color_brewer(palette="Dark2")+
    theme_minimal()+
        theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())



ggplot2::ggsave(file="fig/phy_silva_Abundace.pdf", silva_phy, width = 5, height = 5, dpi = 600)
ggplot2::ggsave(file="fig/phy_silva_Abundace.png", silva_phy, width = 5, height = 5, dpi = 600)
ggplot2::ggsave(file="fig/phy_blast_Abundace.pdf", blast_phy, width = 5, height = 5, dpi = 600)
ggplot2::ggsave(file="fig/phy_blast_Abundace.png", blast_phy, width = 5, height = 5, dpi = 600)

ggplot2::ggsave(file="fig/gen_silva_Abundace.pdf", silva_gen, width = 15, height = 10, dpi = 600)
ggplot2::ggsave(file="fig/gen_silva_Abundace.png", silva_gen, width = 15, height = 10, dpi = 600)
ggplot2::ggsave(file="fig/gen_blast_Abundace.pdf", blast_gen, width = 8, height = 10, dpi = 600)
ggplot2::ggsave(file="fig/gen_blast_Abundace.png", blast_gen, width = 8, height = 10, dpi = 600)

silva_gen

Tpsa@sam_data$glom <- "a"
TPSslv@sam_data$glom <- "a"

Tpsa@sam_data$glom

colnames(TPSslv@tax_table)

psa <- merge_samples(Tpsa, "glom")
PSslv <- merge_samples(TPSslv, "glom")

Ebla <- subset_taxa(psa, superkingdom%in%"Eukaryota") 
Bbla <- subset_taxa(psa, superkingdom%in%"Bacteria")

Bslv <- subset_taxa(PSslv, Kingdom%in%"Bacteria")
Eslv <- subset_taxa(PSslv, Kingdom%in%"Eukaryota")

Ebla <-  psmelt(Ebla)
Eslv <- psmelt(Eslv)

Bbla <- psmelt(Bbla)
Bslv <- psmelt(Bslv)

library(RColorBrewer)

Ebla_phy <- ggplot(Ebla, aes(dpi, Abundance, fill=phylum))+
    geom_bar(stat="identity", position="stack")+
    labs(y="ASV reads (per g of faeces)")+
    scale_fill_brewer(palette="Set3")+
    theme_minimal()+
        theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
#          legend.position = "none",
          axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

Eslv_phy <- ggplot(Eslv, aes(dpi, Abundance,  fill=Order))+
    geom_bar(stat="identity", position="stack")+
    labs(y="ASV reads (per g of faeces)")+
    scale_fill_brewer(palette="Dark2")+
    theme_minimal()+
        theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
#          legend.position = "none",
          axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())


Bbla_phy <- ggplot(Bbla, aes(dpi, Abundance, fill=phylum))+
    geom_bar(stat="identity", position="stack")+
    labs(y="ASV reads (per g of faeces)")+
    scale_fill_brewer(palette="Set3")+
    theme_minimal()+
        theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
#          legend.position = "none",
          axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

Bslv_phy <- ggplot(Bslv, aes(dpi, Abundance, fill=Phylum))+
    geom_bar(stat="identity", position="stack")+
    labs(y="ASV reads (per g of faeces)")+
#    scale_fill_brewer(palette="Dark2")+
    theme_minimal()+
        theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
#          legend.position = "none",
          axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

Ebla_phy

plot_grid(Ebla_phy, Eslv_phy) -> Euk_phy

plot_grid(Bbla_phy, Bslv_phy) -> Bac_phy

Euk_phy

Bac_phy

ggplot2::ggsave(file="fig/Eukphy_blast_silva_Abundace.pdf", Euk_phy, width = 7, height = 5, dpi = 600)

ggplot2::ggsave(file="fig/Eukphy_blast_silva_Abundace.png", Euk_phy, width = 7, height = 5, dpi = 600)

ggplot2::ggsave(file="fig/Bacphy_blast_silva_Abundace.pdf", Bac_phy, width = 7, height = 5, dpi = 600)

ggplot2::ggsave(file="fig/Bacphy_blast_silva_Abundace.png", Bac_phy, width = 7, height = 5, dpi = 600)
