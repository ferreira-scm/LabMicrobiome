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
library(relaimpo)
library(lme4)
library(MASS)
library(DescTools)
library(cocor)
library(reshape2)
library(Hmisc)

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

p <- c(0.6113, 0.00001, 0.00001, 0.2989, 0.0225)
p.adjust(p, method="BH")

Plotting_cor(ps=all.PSwang, "MA_wang", dir="fig/MA/")

p <- c(0.0187, 0.00001, 0.00001, 0.8687, 0.0028)
p.adjust(p, method="BH")

Plotting_cor(ps=sin.PS18S, "SA", dir="fig/SA/")

p <- c(0.4525, 0.00001, 0.0142, 0.0005, 0.3511)
p.adjust(p, method="BH")

##################################################################
### OK, so let's remove food
plant <- subset_taxa(f.all, !phylum%in%"Streptophyta")
plant <- subPS(plant)

Mus <- subset_taxa(f.all, !phylum%in%"Chordata")
Mus <- subPS(Mus)

worms <- subset_taxa(f.all, !phylum%in%"Nematoda")
worms <- subPS(worms)

plant18 <- subset_taxa(f.sin18, !phylum%in%"Streptophyta")
plant18 <- subPS(plant18)

#no Mus in fPS18S
worms18 <- subset_taxa(f.sin18, !phylum%in%"Nematoda")
worms18 <- subPS(worms18)

TSS <- subPS(f.all)
TSS18 <- subPS(f.sin18)

PlantMusWorms <-  subset_taxa(f.all, !(phylum%in%"Streptophyta"| phylum%in%"Nematoda" | phylum%in%"Chordata"))
PlantMusWorms <- subPS(PlantMusWorms)

PlantMus <-  subset_taxa(f.all, !(phylum%in%"Streptophyta"|phylum%in%"Chordata"))
PlantMus <- subPS(PlantMus)

PlantMusWorms18 <-  subset_taxa(f.sin18, !(phylum%in%"Streptophyta"| phylum%in%"Nematoda" | phylum%in%"Chordata"))
PlantMusWorms18 <- subPS(PlantMusWorms18)

plantw <- subset_taxa(f.allwang, !phylum%in%"Streptophyta")
plantw <- subPS(plantw)

wormsw <- subset_taxa(f.allwang, !phylum%in%"Nematoda")
wormsw <- subPS(wormsw)

TSSwang <- subPS(f.allwang)

PlantWormsw <-  subset_taxa(f.allwang, !(phylum%in%"Streptophyta"| phylum%in%"Nematoda"))
PlantWormsw <- subPS(PlantWormsw)

## person tests
cor.test(TSS$logGC, TSS$logA, method="pearson")
cor.test(plant$logGC, plant$logA, method="pearson")
cor.test(worms$logGC, worms$logA, method="pearson")
cor.test(Mus$logGC, Mus$logA, method="pearson")
#cor.test(PlantMusWorms$logGC, PlantMusWorms$logA, method="pearson")
cor.test(PlantMus$logGC, PlantMus$logA, method="pearson")

plant <- plant[,c("labels", "logA")]
names(plant) <- c("labels", "logA_plant")
tss.df <- merge(TSS, plant, by="labels")

worms <- worms[,c("labels", "logA")]
names(worms) <- c("labels", "logA_worms")
tss.df <- merge(tss.df, worms, by="labels")

Mus <- Mus[,c("labels", "logA")]
names(Mus) <- c("labels", "logA_Mus")
tss.df <- merge(tss.df, Mus, by="labels")

PlantMus <- PlantMus[,c("labels", "logA")]
names(PlantMus) <- c("labels", "logA_PlantMus")
tss.df <- merge(tss.df, PlantMus, by="labels")

cocor(~logGC + logA | logGC + logA_plant, data = tss.df,
            test = c("hittner2003", "zou2007"))

cocor(~logGC + logA | logGC + logA_Mus, data = tss.df,
            test = c("hittner2003", "zou2007"))

cocor(~logGC + logA | logGC + logA_worms, data = tss.df,
            test = c("hittner2003", "zou2007"))

cocor(~logGC + logA | logGC + logA_PlantMus, data = tss.df,
            test = c("hittner2003", "zou2007"))

p.adj <- c(0.0008, 0.2735, 0.1159, 0.0002)

cor.test(TSS18$logGC, TSS18$logA, method="pearson")
cor.test(plant18$logGC, plant18$logA, method="pearson")
cor.test(worms18$logGC, worms18$logA, method="pearson")
#cor.test(PlantMusWorms18$logGC, PlantMusWorms18$logA, method="pearson")

plant18 <- plant18[,c("labels", "logA")]
names(plant18) <- c("labels", "logA_plant")
tss.df18 <- merge(TSS18, plant18, by="labels")

worms18 <- worms18[,c("labels", "logA")]
names(worms18) <- c("labels", "logA_worms")
tss.df18 <- merge(tss.df18, worms18, by="labels")

cocor(~logGC + logA | logGC + logA_plant, data = tss.df18,
            test = c("hittner2003", "zou2007"))

cocor(~logGC + logA | logGC + logA_worms, data = tss.df18,
            test = c("hittner2003", "zou2007"))

p <- c(0.0345, 0.6375)
p.adjust(p, method="BH")

cor.test(TSSwang$logGC, TSSwang$logA, method="pearson")
cor.test(plantw$logGC, plantw$logA, method="pearson")
cor.test(wormsw$logGC, wormsw$logA, method="pearson")
#cor.test(PlantWormsw$logGC, PlantWormsw$logA, method="pearson")

plantw <- plantw[,c("labels", "logA")]
names(plantw) <- c("labels", "logA_plant")
tss.dfw <- merge(TSSwang, plantw, by="labels")

wormsw <- wormsw[,c("labels", "logA")]
names(wormsw) <- c("labels", "logA_worms")
tss.dfw <- merge(tss.dfw, wormsw, by="labels")

cocor(~logGC + logA | logGC + logA_plant, data = tss.dfw,
            test = c("hittner2003", "zou2007"))

cocor(~logGC + logA | logGC + logA_worms, data = tss.dfw,
            test = c("hittner2003", "zou2007"))

p <- c(0.0045, 0.3608)
p.adjust(p, method="BH")


### let's expand table 2
PSeimf <-subset_taxa(all.PS, family%in%"Eimeriidae")
PSeimf18 <- subset_taxa(sin.PS18S, family%in%"Eimeriidae")
PSeimfw <- subset_taxa(all.PSwang, family%in%"Eimeriidae")

PSeimf <-subset_taxa(f.all, family%in%"Eimeriidae")

#create total sums and Eimeria sums data frame
df <- data.frame(sample_sums(otu_table(PSeimf)))

df

df$labels <- rownames(df)

head(df)

names(df) <- c("Eimeriidae_f", "labels")

df$logEimeriidae_f <- log(1+df$Eimeriidae_f)

tss.df <- merge(tss.df, df, by="labels")

head(tss.df)

cocor(~logGC + logEimeriidae | logGC + logA, data = tss.df,
            test = c("hittner2003", "zou2007"))

cor.test(tss.df$logGC, tss.df$logEimeriidae, method="pearson")

cor.test(tss.df$logGC, tss.df$logEimeriidae_f, method="pearson")

eimf <-as.data.frame(sample_sums(PSeimf))
eimf$labels <- rownames(eimf)
names(eimf) <- c("EimeriaSums", "labels")
names(df) <- c("TotalSums", "labels")


cocor(~logGC + logA | logGC + logA_plant, data = tss.dfw,
            test = c("hittner2003", "zou2007"))
p <- c(0.6113, 0.00001, 0.00001, 0.2989, 0.0225)
p.adjust(p, method="BH")


ps=all.PS
name="MA"
dir="fig/MA/"

Plotting_cor(ps=all.PSwang, "MA_wang", dir="fig/MA/")

p <- c(0.0187, 0.00001, 0.00001, 0.8687, 0.0028)
p.adjust(p, method="BH")

Plotting_cor(ps=sin.PS18S, "SA", dir="fig/SA/")

p <- c(0.4525, 0.00001, 0.0142, 0.0005, 0.3511)
p.adjust(p, method="BH")


# plotting TSS correlation
a <- p_tss(TSS, "a)", "MA")
b <- p_tss(plant, "b)", "MA no plant")
c <- p_tss(Mus, "c)", "MA no host")
d <- p_tss(worms, "d)", "MA no nematodes")
#e <- p_tss(PlantMusWorms, "e)", "MA no plants, host or nematodes")
e <- p_tss(PlantMus, "e)", "MA no plants or host")

plot_grid(a,b,c,d,e) -> p_cor

MA.a <- lm(data=TSS,logGC~logA)
MA.b <- lm(data=plant,logGC~logA)
MA.c <- lm(data=Mus,logGC~logA)
MA.d <- lm(data=worms,logGC~logA)
MA.e <- lm(data=PlantMus,logGC~logA)

print(cocor(~logGC + logFilEimeriaSums | logGC + logEimeriaSums, data = cor.df,
            test = c("hittner2003", "zou2007")))


summary(MA.e)


ggplot2::ggsave(file="fig/MA/Biological_rem_MA.pdf", p_cor, width = 15, height = 8, dpi = 600)
ggplot2::ggsave(file="fig/MA/Biological_rem_MA.png", p_cor, width = 15, height = 8, dpi = 600)

a1 <- p_tss(TSS18, "a)", "SA")
b1 <- p_tss(plant18, "b)", "SA no plant")
#c1 <- p_tss(Mus18, "c)", "SA no host")
d1 <- p_tss(worms18, "c)", "SA no nematodes")
#e1 <- p_tss(PlantMusWorms18, "d)", "SA no plants or nematodes")

SA.a <- lm(data=TSS18,logGC~logA)
SA.b <- lm(data=plant18,logGC~logA)
SA.d <- lm(data=worms18,logGC~logA)

summary(SA.d)

plot_grid(a1,b1,d1,e1) -> p_cor1
ggplot2::ggsave(file="fig/SA/Biological_rem_SA.pdf", p_cor1, width = 15, height = 10, dpi = 600)
ggplot2::ggsave(file="fig/SA/Biological_rem_SA.png", p_cor1, width = 15, height = 10, dpi = 600)

a2 <- p_tss(TSSwang, "a)", "MA wang")
b2 <- p_tss(plantw, "b)", "MA wang no plant")
#c1 <- p_tss(Mus18, "c)", "SA no host")
d2 <- p_tss(wormsw, "c)", "MA wang no nematodes")
e2 <- p_tss(PlantWormsw, "d)", "MA wang no plants or nematodes")

MAw.a <- lm(data=TSSwang,logGC~logA)
MAw.b <- lm(data=plantw,logGC~logA)
MAw.d <- lm(data=wormsw,logGC~logA)

summary(MAw.d)

plot_grid(a2,b2,d2,e2) -> p_cor2
ggplot2::ggsave(file="fig/MA/Biological_rem_MA_wang.pdf", p_cor2, width = 15, height = 10, dpi = 600)
ggplot2::ggsave(file="fig/MA/Biological_rem_MA_wang.png", p_cor2, width = 15, height = 10, dpi = 600)



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


f.all.l <- list()

for (i in 1:48) {
    try(f.all.l[[i]] <- fil(all.PS.l[[i]]), silent=TRUE)
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

################################### Mus sequences
for (i in 1:48) {
    print(names(all.PS.l)[i])
    try(p <- subset_taxa(f.all.l[[i]],genus=="Mus"), silent=TRUE)
    if (exists("p")) {
        try(a <- rownames(p@tax_table), silent=TRUE)
        print(paste(i, "- ", names(f.all.l[[i]]), ": ", length(a), sep=""))
        print(a)
}
    rm(p)
}

mus.s <- c(rownames(subset_taxa(f.all.l[[12]]@tax_table, genus=="Mus")), rownames(subset_taxa(f.all.l[[33]]@tax_table, genus=="Mus")))


library(ShortRead)
writeFasta(DNAStringSet(mus.s), "tmp/Mus_ASV.fasta")

#######################
# how does host DNA change with infection
# no Mus reads in SA
PSclr = microbiome::transform(f.all, transform="clr")
#PSeimf <- subset_taxa(PSclr, family%in%"Eimeriidae")

Mus <- subset_taxa(T.all, genus%in%"Mus")

Mus@tax_table

Mus <-aggregate_taxa(Mus, level="genus")
Mus.clr <- subset_taxa(PSclr, genus%in%"Mus")
Mus.clr <-aggregate_taxa(Mus.clr, level="genus")

nb <- length(levels(as.factor(m$EH_ID)))

coul <- colorRampPalette(brewer.pal(8, "Accent"))(nb)

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


mp.clr

# Are host abundance and Eimeria associated
cor.test(log(1+m$Abundance), log(1+m$Genome_copies_gFaeces), method="spearman")

cor.test(log(1+m$Abundance), log(1+m$OPG), method="spearman")

cor.test(m.clr$Abundance, log(1+m.clr$Genome_copies_gFaeces), method="spearman")
#cor.test(m.clr$Abundance, log(1+m.clr$OPG), method="spearman")

m.host <- ggplot(m, aes(log(1+Abundance), log(1+Genome_copies_gFaeces)))+
    geom_point(aes(fill=dpi), shape=21, size=4, alpha=0.6)+
    scale_fill_brewer(palette="Spectral")+           
    labs(x="Host reads per g of Faeces, log(+1)", y="Genome copies per g of faeces log(1+)")+
    labs(tag= "b)")+
#    geom_smooth(method = "lm", se=FALSE, na.rm=TRUE) +
#    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE,  label.x=0.1, label.y=0.9) +  
    annotate(geom="text", x=5, y=19, label="Spearman rho=0.41, p<0.001")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          legend.position = "top",
          axis.line = element_line(colour = "black"))+
    guides(fill=guide_legend(nrow=1, byrow=TRUE))


m.clr.host <- ggplot(m.clr, aes(Abundance, log(1+Genome_copies_gFaeces)))+
    geom_point(aes(fill=dpi), shape=21, size=4, alpha=0.6)+
#    scale_color_brewer(palette="Accent")+
    scale_fill_brewer(palette="Spectral")+
#    geom_smooth(method = "lm", se=FALSE, na.rm=TRUE) +
#    stat_poly_eq(aes(label = paste(..eq.label..,
#                                   ..rr.label..,
#                                   sep = "~~~")),
#                 parse = TRUE, label.x=0.85, label.y=0.23) + 
    annotate(geom="text", x=6, y=3, label="Spearman rho=0.44, p<0.001")+
#    geom_line(aes(group=EH_ID), alpha=0.3)+
    labs(x="Host reads per g of Faeces, clr(+1)", y="Genome copies per g of faeces log(1+)")+
    labs(tag= "b)")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          legend.position = "top",
          axis.line = element_line(colour = "black"))+
    guides(fill=guide_legend(nrow=1, byrow=TRUE))

# save plots of what we have so far
plot_grid(mp.acs, m.host, nrow=2) -> m.acsCor
plot_grid(mp.clr,m.clr.host, nrow=2) -> m.clrCor

ggplot2::ggsave(file="fig/Mus_Eimeria.clr.pdf", m.clrCor, width = 8, height = 8, dpi = 600)
ggplot2::ggsave(file="fig/Mus_Eimeria.clr.png", m.clrCor, width = 8, height = 8, dpi = 600)

ggplot2::ggsave(file="fig/Mus_Eimeria.acs.pdf", m.acsCor, width = 8, height = 8, dpi = 600)
ggplot2::ggsave(file="fig/Mus_Eimeria.acs.png", m.acsCor, width = 8, height = 8, dpi = 600)


# ordination
# host reads, total dna, opg, Genome_copies
# can't have zeros in all columns

m.ord <- m[,c("Genome_copies_gFaeces", "Abundance", "OPG")][rowSums(m[,c("Genome_copies_gFaeces", "Abundance", "OPG")])>0,]

m.clr.ord <- m.clr[,c("Genome_copies_gFaeces", "Abundance", "OPG", "DNA_g_feces")][rowSums(m[,c("Genome_copies_gFaeces", "Abundance", "OPG", "DNA_g_feces")])>0,]


library(vegan)

m.mds <- metaMDS(m.ord)
m.clr.euc <- vegdist(m.clr.ord, method="euclidean")
m.bray <- vegdist(m.ord, "bray")
#sanity check
rownames(mAxis)==rownames(m.ord)

dpi <- as.factor(m$dpi[as.numeric(rownames(m.ord))])
dpi.clr <- as.factor(m.clr$dpi[as.numeric(rownames(m.clr.ord))])

weighloss <- m$weightloss[as.numeric(rownames(m.ord))]
weighloss.clr <- m.clr$weightloss[as.numeric(rownames(m.clr.ord))]

adonis2(m.bray~dpi+weighloss, method="bray", permutations=10000)
adonis2(m.clr.euc~dpi+weighloss.clr, method="euclidean", permutations=10000)

#m.anosim <- anosim(m.bray,dpi, distance="bray", permutations=1000)
#summary(m.anosim)

#mAxis <- m.mds$points[,1:2]
#mAxis <- as.data.frame(mAxis)
#ggplot(mAxis, aes(MDS1, MDS2))+
#    geom_point(size=2, alpha=0.8, aes(color=dpi))

ord.fit <- envfit(m.mds~dpi+weighloss)
ord.fit

ord.fit <- envfit(m.bray~dpi+weighloss)
ord.fit

#mmus <- lm(Abundance~Genome_copies_gFaeces+ DNA_g_feces+dpi, data=m.clr)

m.clr$logGC <- log(1+m.clr$Genome_copies_gFaeces)

m$logA <- log(1+m$Abundance)
m$logGC <- log(1+m$Genome_copies_gFaeces)

##CLR
mmus <- lm(Abundance~logGC+ DNA_g_feces, data=m.clr)
summary(mmus)
mwl <- lm(weightloss~Abundance+Genome_copies_gFaeces+ DNA_g_feces, data=m.clr)
summary(mwl)
anova(mwl)

cm1 <- lmer(weightloss~Abundance+Genome_copies_gFaeces +dpi+(1|EH_ID), data=m.clr)
cm2 <- lmer(weightloss~Abundance+Genome_copies_gFaeces +(1|EH_ID), data=m.clr)
cm3 <- lmer(weightloss~Abundance+dpi+(1|EH_ID), data=m.clr)
cm4 <- lmer(weightloss~Genome_copies_gFaeces +dpi+(1|EH_ID), data=m.clr)

## ACS
#mwl.acs <- lm(weightloss~logA+logGC, data=m)
#plot(mwl.acs)
mwl.acs <- lm(weightloss~Abundance+Genome_copies_gFaeces, data=m)

calc.relimp(mwl.acs)

mmwl.acs1 <- lmer(weightloss~Abundance+Genome_copies_gFaeces +dpi+(1|EH_ID), data=m)

mmwl.acs2 <- lmer(weightloss~Abundance+Genome_copies_gFaeces+(1|EH_ID), data=m)

mmwl.acs3 <- lmer(weightloss~Abundance +dpi+(1|EH_ID), data=m)

mmwl.acs4 <- lmer(weightloss~Genome_copies_gFaeces +dpi+(1|EH_ID), data=m)

sink("fig/weightLoss_lmm.txt")
summary(mmwl.acs1)
sink()

mmwl.acs <- lmer(weightloss~Abundance+Genome_copies_gFaeces +(1|EH_ID)+(1|dpi), data=m)

ranova(mmwl.acs1)

anova(mmwl.acs1, mmwl.acs2)
anova(mmwl.acs1,mmwl.acs3)
anova(mmwl.acs1,mmwl.acs4)


anova(mmwl.acs)

mmmus <- lmer(Abundance~Genome_copies_gFaeces+(1|EH_ID)+(1|dpi), data=m)

mmmus0 <- lmer(Abundance~1+(1|EH_ID), data=m)

summary(mmmus)

library(lmerTest)

ranova(mmmus)

anova(mmmus,mmmus0)


# how does DNA g/faeces change with infection
cor.test(m$logA, log(1+m$weightloss), method="spearman")


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


######################################################

cor.p <- p.adjust(mat.cor$P, method="BH")
print(cor.p)



mat.cor <- rcorr(as.matrix(cor.df), type="pearson") 
mat.cor$r <- abs(mat.cor$r)
upper_tri <- get_upper_tri(mat.cor$r)
upper_tri
m.mat <- melt(upper_tri)

library(ggplot2)
ggplot(data = m.mat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "#006633", mid = "white",
                               midpoint = 0, limit = c(0,1), space = "Lab",
                               name="Pearson\nCorrelation") +
        theme_minimal()+
        theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                                size = 12, hjust = 1))+
             coord_fixed()

