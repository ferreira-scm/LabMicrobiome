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
library(RColorBrewer)

## using the devel
#devtools::load_all("/SAN/Susanas_den/MultiAmplicon/")

set.seed(500)
source("bin/PlottingCor.R")

## All runs pooled
all.PS <- readRDS("tmp/PhyloSeqData_All.Rds")
all.PS.l <- readRDS("tmp/PhyloSeqList_All.Rds")
# single amplicon from MA run
all.PSwang <- all.PS.l[[37]]

## Single amplicon 18S
sin.PS18S <- readRDS("tmp/PS_18Swang.Rds")
sin.PS18S.slv <- readRDS("tmp/PS_18Swang_SILVA.Rds")

## Single amplicon "pooled"
sin.PS <- readRDS("tmp/PhyloSeqData18S.Rds")
sin.PS.slv <- readRDS("tmp/PhyloSeqData18S_SILVA.Rds")

source("bin/PlottingCor.R")
# let's filter
f.sin18 <- fil(sin.PS18S)
f.all <- fil(all.PS)
f.allwang <- fil(all.PSwang)
f.sin18.slv <- fil(sin.PS18S.slv)
## Single amplicon "pooled"
f.sin <- fil(sin.PS)
f.sin.slv <- fil(sin.PS.slv)
# and transform
T.sin18 <- f.sin18
T.all <- f.all
T.allwang <-f.allwang
T.sin18.slv <- f.sin18.slv
T.sin <-f.sin
T.sin.slv <-f.sin.slv
otu_table(T.sin18) <- otu_table(T.sin18)*sample_data(T.sin18)$DNA_g_feces
otu_table(T.all) <- otu_table(T.all)*sample_data(T.all)$DNA_g_feces
otu_table(T.allwang) <- otu_table(T.allwang)*sample_data(T.allwang)$DNA_g_feces
otu_table(T.sin18.slv) <- otu_table(T.sin18.slv)*sample_data(T.sin18.slv)$DNA_g_feces
otu_table(T.sin) <- otu_table(T.sin)*sample_data(T.sin)$DNA_g_feces
otu_table(T.sin.slv) <- otu_table(T.sin.slv)*sample_data(T.sin.slv)$DNA_g_feces

#all.TSS <- transform_sample_counts(f.all, function(x) x / sum(x))

## OK, now we want all the Eimeria sequences
Eim <- subset_taxa(T.all, family%in%"Eimeriidae")
Eim2 <- subset_taxa(T.sin18, family%in%"Eimeriidae")
Eim.slv <- subset_taxa(T.sin18.slv, Family%in%"Eimeriorina")

# load silva taxonomic annotation
#PSslv <- readRDS("tmp/PhyloSeqData18S_SILVA.Rds")
#PS18slv <- readRDS("tmp/PS_18Swang_SILVA.Rds")

# now plotting
Plotting_cor(ps=all.PS, "MA", dir="fig/MA/")

Plotting_cor(ps=all.PSwang, "MA_wang", dir="fig/MA/")

Plotting_cor(ps=sin.PS18S, "SA", dir="fig/SA/")

## now plotting SA silva
#Plotting_cor(ps=PSslv, "SA_slv", dir="fig/SA/")

#Plotting_cor(ps=PSwang, name="wang1141_13_F.Nem_0425_6_3_R", dir="fig/MA/")

for (i in 1:48) {
    nm <- names(all.PS.l)[i]
    ps <- all.PS.l[[i]]
    print(nm)
    try(Plotting_cor(ps, name=nm, dir="fig/MA/"))
}

for (i in 1:48) {
    nm <- names(all.PS.l)[i]
    ps <- all.PS.l[[i]]
    try(NoFilPlotting_cor(ps, name=nm, dir="fig/MA/NoFil/"))
}

#how many primers amplify Apicomplexa and which families?
for (i in 1:48) {
#    print(names(PS.l)[i])
    try(p <- subset_taxa(all.PS.l[[i]],phylum=="Apicomplexa"), silent=TRUE)
    try(get_taxa_unique(p, "family"), silent=TRUE)
    if (exists("p")) {
        a <- get_taxa_unique(p, "family")
        print(paste(i, "- ", names(all.PS.l[i]), ": ", length(a), sep=""))
        print(a)
}
    rm(p)
}


for (i in 1:48) {
#    print(names(all.PS.l)[i])
    try(p <- subset_taxa(all.PS.l[[i]],family=="Eimeriidae"), silent=TRUE)
    try(get_taxa_unique(p, "family"), silent=TRUE)
    if (exists("p")) {
        a <- get_taxa_unique(p, "family")
        print(paste(i, "- ", names(all.PS.l[i]), ": ", length(a), sep=""))
        print(a)
}
    rm(p)
}

###################################
all.PS.l[14]
all.PS.l[34]
all.PS.l[37]
all.PS.l[45]

#######################
# how does host DNA change with infection
# no Mus reads in SA
PSclr = microbiome::transform(f.all, transform="clr")
#PSeimf <- subset_taxa(PSclr, family%in%"Eimeriidae")

Mus <- subset_taxa(T.all, genus%in%"Mus")
Mus <-aggregate_taxa(Mus, level="genus")
Mus.clr <- subset_taxa(PSclr, genus%in%"Mus")
Mus.clr <-aggregate_taxa(Mus.clr, level="genus")

nb <- length(levels(as.factor(m$EH_ID)))

coul <- colorRampPalette(brewer.pal(12, "Accent"))(nb)

mycol <- coul[as.numeric(as.factor(m$EH_ID))]

m <- psmelt(Mus)
m.clr <- psmelt(Mus.clr)


mp.acs <- ggplot(m, aes(dpi, Abundance, color=EH_ID))+
    geom_point(aes(fill=EH_ID), shape=21, size=4, position=position_jitter(0.1), alpha=0.6)+
    scale_color_manual(values=coul)+
        scale_fill_manual(values=coul)+           
    geom_line(aes(group=EH_ID))+
    labs(x="Days post infection", y="Host reads (per g of faeces)")+
    labs(tag= "a)")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
           legend.position = "none",
          axis.line = element_line(colour = "black"))

mp.clr <- ggplot(m.clr, aes(dpi, Abundance, color=EH_ID))+
    geom_point(aes(fill=EH_ID), shape=21, size=4, position=position_jitter(0.1), alpha=0.6)+
    scale_color_manual(values=coul)+
        scale_fill_manual(values=coul)+           
    geom_line(aes(group=EH_ID))+
    labs(x="Days post infection", y="Host reads (per g of faeces)")+
    labs(tag= "a)")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
           legend.position = "none",
          axis.line = element_line(colour = "black"))


# save plots of what we have so far
plot_grid(mp.acs,mp.clr) -> fCor

ggplot2::ggsave(file="fig/Mus_abundance.pdf", fCor, width = 8, height = 5, dpi = 600)
ggplot2::ggsave(file="fig/Mus_abundance.png", fCor, width = 8, height = 5, dpi = 600)


# Are host abundance and Eimeria associated
cor.test(m$Abundance, m$Genome_copies_gFaeces, method="pearson")
cor.test(m$Abundance, m$OPG, method="pearson")

ggplot(m, aes(Abundance, OPG))+
    geom_point()

ggplot(m, aes(Abundance, Genome_copies_gFaeces))+
    geom_point()

mp

m$logA <- log10(1+m$Abundance)

m$logGC <- log10(1+m$Genome_copies_gFaeces)

m$logDNA <- log10(1+m$DNA_g_feces)

#mmus <- lm(Abundance~Genome_copies_gFaeces+ DNA_g_feces+dpi, data=m.clr)

m$OPG

mmus <- lm(Abundance~Genome_copies_gFaeces+ DNA_g_feces+m$OPG + dpi, data=m.clr)
summary(mmus)

mwl <- glm(weightloss~Abundance+Genome_copies_gFaeces+ DNA_g_feces+dpi, data=m.clr)

mwl <- glm(weightloss~logA+logGC+ DNA_g_feces+dpi, data=m)

summary(mmus)

summary(mwl)

library(lme4)

library(MASS)

m2 <- lmer(weightloss~logA*logGC+ logDNA +(1|dpi), data=m)

m2 <- lmer(weightloss~logA+logGC +EH_ID+(1|dpi), data=m)

m2 <- glmm(weightloss~logA+logGC +EH_ID+(1|dpi), data=m)

summary(m2)

Sdem <- with(m, data.frame(
                             weightloss=weightloss- ave(weightloss, EH_ID),
                             logGC=logGC-ave(logGC, EH_ID, FUN=function(x) mean(x, na.rm=T)),
                             logA=logA-ave(logA, EH_ID, FUN=function(x) mean(x, na.rm=T)),
                             dpi=dpi,
                             EH_ID=EH_ID))


STdem <- with(Sdem, data.frame(
                            weightloss=weightloss- ave(weightloss, dpi),
                            logGC=logGC-ave(logGC, dpi, FUN=function(x) mean(x, na.rm=T)),
                            logA=logA-ave(logA, dpi, FUN=function(x) mean(x, na.rm=T)),
                            dpi=dpi,
                            EH_ID=EH_ID))

head(m)

## analysing at peak
# it isn't working so well here...
library(dplyr)

MaxDNA <- m %>% dplyr::group_by(EH_ID) %>%
    dplyr::filter(logGC==max(logGC, na.rm=T))%>%
    dplyr::select(EH_ID, logDNA, logA, logGC)%>% data.frame()

Maxwl <- m %>% dplyr::group_by(EH_ID) %>%
    dplyr::filter(weightloss==max(weightloss, na.rm=T))%>%
    dplyr::select(EH_ID, dpi, weightloss)%>% data.frame()

Maxdf <- merge(Maxwl, MaxDNA, by="EH_ID")

mdem <- lm(weightloss~logA*logGC, data=m)

summary(mdem)

maxm <- lm(weightloss~logA*logGC, data=Maxdf)

summary(maxm)

calc.relimp(mdem)

calc.relimp(maxm)

ggplot(m, aes(Abundance, DNA_g_feces))+
          geom_point()

ggplot(m, aes(Abundance, Genome_copies_gFaeces))+
          geom_point()

# how does DNA g/faeces change with infection
dnap <- ggplot(m, aes(y=DNA_g_feces, x=dpi, color=EH_ID))+
    geom_point()+
    geom_line(aes(group=EH_ID))+
    labs(x="Days post infection", y="DNA (per g of faeces)")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none")


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


##################################################################
### OK, so let's remove food
plant <- subset_taxa(fPS, !phylum%in%"Streptophyta")
plant <- subPS(plant)

Mus <- subset_taxa(fPS, !phylum%in%"Chordata")
Mus <- subPS(Mus)

worms <- subset_taxa(fPS, !phylum%in%"Nematoda")
worms <- subPS(worms)

plant18 <- subset_taxa(fPS18S, !phylum%in%"Streptophyta")
plant18 <- subPS(plant18)

#no Mus in fPS18S
worms18 <- subset_taxa(fPS18S, !phylum%in%"Nematoda")
worms18 <- subPS(worms18)

TSS <- subPS(fPS)
TSS18 <- subPS(fPS18S)

PlantMusWorms <-  subset_taxa(fPS, !(phylum%in%"Streptophyta"| phylum%in%"Nematoda" | phylum%in%"Chordata"))
PlantMusWorms <- subPS(PlantMusWorms)

PlantMus <-  subset_taxa(fPS, !(phylum%in%"Streptophyta"|phylum%in%"Chordata"))
PlantMus <- subPS(PlantMus)

PlantMusWorms18 <-  subset_taxa(fPS18S, !(phylum%in%"Streptophyta"| phylum%in%"Nematoda" | phylum%in%"Chordata"))
PlantMusWorms18 <- subPS(PlantMusWorms18)

plantw <- subset_taxa(fPSwang, !phylum%in%"Streptophyta")
plantw <- subPS(plantw)

wormsw <- subset_taxa(fPSwang, !phylum%in%"Nematoda")
wormsw <- subPS(wormsw)

TSSwang <- subPS(fPSwang)

PlantWormsw <-  subset_taxa(fPSwang, !(phylum%in%"Streptophyta"| phylum%in%"Nematoda"))
PlantWormsw <- subPS(PlantWormsw)


cor.test(TSS$logGC, TSS$logA, method="pearson")

cor.test(plant$logGC, plant$logA, method="pearson")

cor.test(worms$logGC, worms$logA, method="pearson")

cor.test(Mus$logGC, Mus$logA, method="pearson")

cor.test(PlantMusWorms$logGC, PlantMusWorms$logA, method="pearson")

cor.test(PlantMus$logGC, PlantMus$logA, method="pearson")


cor.test(TSS18$logGC, TSS18$logA, method="pearson")

cor.test(plant18$logGC, plant18$logA, method="pearson")

cor.test(worms18$logGC, worms18$logA, method="pearson")

cor.test(PlantMusWorms18$logGC, PlantMusWorms18$logA, method="pearson")


cor.test(TSSwang$logGC, TSSwang$logA, method="pearson")

cor.test(plantw$logGC, plantw$logA, method="pearson")

cor.test(wormsw$logGC, wormsw$logA, method="pearson")

cor.test(PlantWormsw$logGC, PlantWormsw$logA, method="pearson")



# plotting TSS correlation
a <- p_tss(TSS, "a)", "MA")
b <- p_tss(plant, "b)", "MA no plant")
c <- p_tss(Mus, "c)", "MA no host")
d <- p_tss(worms, "d)", "MA no nematodes")
e <- p_tss(PlantMusWorms, "e)", "MA no plants, host or nematodes")
f <- p_tss(PlantMus, "f)", "MA no plants or host")

plot_grid(a,b,c,d,e, f) -> p_cor

ggplot2::ggsave(file="fig/MA/Biological_rem_MA.pdf", p_cor, width = 15, height = 8, dpi = 600)
ggplot2::ggsave(file="fig/MA/Biological_rem_MA.png", p_cor, width = 15, height = 8, dpi = 600)

a1 <- p_tss(TSS18, "a)", "SA")
b1 <- p_tss(plant18, "b)", "SA no plant")
#c1 <- p_tss(Mus18, "c)", "SA no host")
d1 <- p_tss(worms18, "c)", "SA no nematodes")
e1 <- p_tss(PlantMusWorms18, "d)", "SA no plants or nematodes")


plot_grid(a1,b1,d1,e1) -> p_cor1
ggplot2::ggsave(file="fig/SA/Biological_rem_SA.pdf", p_cor1, width = 15, height = 10, dpi = 600)
ggplot2::ggsave(file="fig/SA/Biological_rem_SA.png", p_cor1, width = 15, height = 10, dpi = 600)

a2 <- p_tss(TSSwang, "a)", "MA wang")
b2 <- p_tss(plantw, "b)", "MA wang no plant")
#c1 <- p_tss(Mus18, "c)", "SA no host")
d2 <- p_tss(wormsw, "c)", "MA wang no nematodes")
e2 <- p_tss(PlantWormsw, "d)", "MA wang no plants or nematodes")

plot_grid(a2,b2,d2,e2) -> p_cor2
ggplot2::ggsave(file="fig/MA/Biological_rem_MA_wang.pdf", p_cor2, width = 15, height = 10, dpi = 600)
ggplot2::ggsave(file="fig/MA/Biological_rem_MA_wang.png", p_cor2, width = 15, height = 10, dpi = 600)

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


silva_gen <- ggplot(Tslv, aes(dpi, Abundance, fill=fct_reorder(Genus, Abundance)))+
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


