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
library(RColorBrewer)
library(microbiomeutilities) # messes up with ggplot2::alpha

## using the devel
#devtools::load_all("/SAN/Susanas_den/MultiAmplicon/")

library(forcats)

source("bin/4_MA_SA_filtering.R")

allTSS

sin18TSS

T.all

## going to plot with ACS
tall <- psmelt(T.all)

head(tall)
f.all.lp

phyTSS <- aggregate_top_taxa2(allTSS, level="phylum", 11)
                                        #genTSS <- aggregate_top_taxa2(allTSS, level="genus", 11)

sin18TSS@tax_table[1,6]

genTSS <- tax_glom(sin18TSS, taxrank="genus")
#genPS <- aggregate_top_taxa2(T.all, level="genus", 11)

gen <- psmelt(genPS)

gen.p <- tax_glom(sin.p, taxrank="genus")

gen.all <- tax_glom(all.p, taxrank="genus")

gen.p  <- transform_sample_counts(gen.p, function(x) x / sum(x))

gen.all  <- transform_sample_counts(gen.all, function(x) x / sum(x))

gen.pp <- psmelt(gen.p)

gen.allp <- psmelt(gen.all)

gen18TSS <- psmelt(gen18TSS)

genTSS <- psmelt(genTSS)

summary(gen.pp$phylum=="Unknown")

sin.p@tax_table[,3]

gen.all

nb.cols=14
mycolor=colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)

composition_dpi_18S <- ggplot(gen.pp, aes(x=Sample, y=Abundance, fill=fct_reorder(genus, Abundance)))+
    geom_bar(position="stack", stat="identity")+
#    scale_x_discrete(labels=gen.pp$EH_ID, breaks=gen.pp$Sample)+
    scale_fill_manual("genus", values=mycolor)+
    facet_grid(~dpi, scales="free")+
    labs(fill="Genus", y="Proportion within total sequence")+
#    rremove("x.text")+
    theme_bw(base_size=14)+
    theme(axis.text.y = element_text(colour = 'black', size = 14),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.key = element_blank(),
#        strip.background = element_rect(colour="black", fill="white"),
        legend.text = element_text(colour = 'black', size = 12, face = 'italic')) 

composition_dpi_all <- ggplot(gen.allp, aes(x=Sample, y=Abundance, fill=fct_reorder(phylum, Abundance)))+
    geom_bar(position="stack", stat="identity")+
#    scale_x_discrete(labels=gen$EH_ID, breaks=gen$Sample)+
#    scale_fill_manual("genus", values=mycolor)+
    facet_grid(~dpi, scales="free")+
    labs(fill="Phylum", y="Proportion within total sequence")+
#    rremove("x.text")+
    theme_bw(base_size=14)+
    theme(axis.text.y = element_text(colour = 'black', size = 14),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.key = element_blank(),
#        strip.background = element_rect(colour="black", fill="white"),
        legend.text = element_text(colour = 'black', size = 12)) 

composition_dpi_all

heat_gen <- ggplot(genTSS, aes(x=Sample, y=genus))+
    geom_tile(aes(fill=Abundance))+
    scale_fill_distiller("Proportion within total sequence", palette = "Spectral")+
    facet_grid(~dpi, scales="free")+
    labs(y="Genus")+
    theme_bw(base_size=14)+
    theme(axis.text.y = element_text(colour = 'black', size = 14, face = 'italic'),
        axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.key = element_blank(),
        strip.background = element_rect(colour="black", fill="white"),
        legend.position="bottom") 

heat_gen

ggsave("fig/Genus_18S_compostion.pdf", composition_dpi_18S, height=6, width=8, dpi=400)
ggsave("fig/Genus_18S_compostion.png", composition_dpi_18S, height=6, width=8, dpi=400)

ggsave("fig/Phylum_18S_16S_compostion.pdf", composition_dpi_18S, height=6, width=8, dpi=400)
ggsave("fig/Phylum_18S_16S_compostion.png", composition_dpi_18S, height=6, width=8, dpi=400)


ggsave("fig/Genus_18Sheat_composition.pdf", heat_gen, height=8, width=8, dpi=400)
ggsave("fig/Genus_18Sheat_composition.png", heat_gen, height=8, width=8, dpi=400)


#ggsave("fig/Genus_heat_composition.png", heat_gen, height=10, width=10, dpi=400)
#ggsave("fig/Genus_heat_compositionTSS.png", heat_genTSS, height=10, width=10, dpi=400)

# 20 most abundant ASV's
dpiPS <- merge_samples(T.all, "dpi")

T.all@sam_data$fakeC <- "little_cheat"

meanPS <- merge_samples(T.all, "fakeC")

meanPS

IDPS <- merge_samples(T.all, "EH_ID")
#phyDPI <- aggregate_taxa(dpiPS, level = "phylum")

m.melt <- psmelt(meanPS)

id.melt <- psmelt(IDPS)
dpi.melt <- psmelt(dpiPS)

abundTaxa <- m.melt$OTU[1:50] # 50 most abundant ASV'S

tall50 <- tall[tall$OTU%in%abundTaxa,]

# fixing genus name
for (i in which(is.na(tall50$genus))){
    tall50$genus[i] <- paste("f", tall50$family[i], sep="__")
}

for (i in which(is.na(tall50$family))){
    tall50$genus[i] <- paste("p", tall50$phylum[i], sep="__")
}

top50asv <- ggplot(tall50, aes(x=Sample, y=genus))+
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

top50asv


ggsave("fig/Genus_top50ASV_heat_composition.pdf", top50asv, height=10, width=25, dpi=600)


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

############# for Emanuel, not reproducable

myt <- data.frame(sdtc$logTSS_Eim, sdtc$logGC, sdtc$dpi)
myt <- na.omit(myt)

names(myt) <- c("logTSS", "logGC", "dpi")

head(myt)

RA <-ggplot(myt, aes(x=logTSS, y=logGC))+
    geom_jitter(shape=21, position=position_jitter(0.002),
                size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    scale_fill_brewer(palette="Spectral")+
    geom_smooth(method = "lm", se=TRUE, na.rm=TRUE, colour="black") +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                                          parse = TRUE) +  
    ylab("Genome copies/ngDNA (log)")+
    xlab("Proportion of Eimeria sequence (log)")+
#    ggtitle("Total sum scaling")+
#        labs(tag= "c)")+
#        annotate(geom="text", x=13, y=0.5, label="Spearman rho=0.94, p<0.001")+
    theme_bw()+
    theme(text = element_text(size=16),
#          axis.title.x = element_blank(),
          legend.position= "bottom")+
    guides(fill=guide_legend(nrow=2, byrow=TRUE))



ggplot2::ggsave(file="fig/EMANUEL_SA.pdf", RA, width = 8, height = 8, dpi = 400)
ggplot2::ggsave(file="fig/EMANUEL_SA.png", RA, width = 8, height = 8, dpi = 400)
