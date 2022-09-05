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
## filtering MA by amplicon
f.all.l <- list()
for (i in 1:48) {
    try(f.all.l[[i]] <- fil(all.PS.l[[i]]), silent=TRUE)
}

f.all.lp <- f.all.l[[1]]
for (i in 2:47){
    f.all.lp <- try(merge_phyloseq(f.all.lp,f.all.l[[i]]))
    print(f.all.lp)}


# sanity check
#merge_phyloseq(f.all.l[[45]], f.all.l[[44]], f.all.l[[47]], f.all.l[[37]],f.all.l[[1]], f.all.l[[7]],
#               f.all.l[[12]], f.all.l[[13]], f.all.l[[13]], f.all.l[[22]], f.all.l[[23]], f.all.l[[28]],
#               f.all.l[[33]], f.all.l[[35]])
## Single amplicon "pooled"
f.sin <- fil(sin.PS)
f.sin.slv <- fil(sin.PS.slv)
# and transform
T.sin18 <- f.sin18
T.all <- f.all
T.all.l <- f.all.lp
T.allwang <-f.allwang
T.sin18.slv <- f.sin18.slv
T.sin <-f.sin
T.sin.slv <-f.sin.slv
otu_table(T.sin18) <- otu_table(T.sin18)*sample_data(T.sin18)$DNA_g_feces
otu_table(T.all) <- otu_table(T.all)*sample_data(T.all)$DNA_g_feces
otu_table(T.all.l) <- otu_table(T.all.l)*sample_data(T.all.l)$DNA_g_feces
otu_table(T.allwang) <- otu_table(T.allwang)*sample_data(T.allwang)$DNA_g_feces
otu_table(T.sin18.slv) <- otu_table(T.sin18.slv)*sample_data(T.sin18.slv)$DNA_g_feces
otu_table(T.sin) <- otu_table(T.sin)*sample_data(T.sin)$DNA_g_feces
otu_table(T.sin.slv) <- otu_table(T.sin.slv)*sample_data(T.sin.slv)$DNA_g_feces

############################################
##############################################
#### exploring MA
for (i in 1:48) {
    nm <- names(all.PS.l)[i]
    ps <- all.PS.l[[i]]
    print(nm)
    try(Plotting_cor(ps, name=nm, dir="fig/MA/"))
}

#for (i in 1:48) {
#    nm <- names(all.PS.l)[i]
#    ps <- all.PS.l[[i]]
#    try(NoFilPlotting_cor(ps, name=nm, dir="fig/MA/NoFil/"))
#}


#########################################################
#how many primers amplify Apicomplexa and which families?
for (i in 1:48) {
#    print(names(all.PS.l)[[i]])
    try(p <- subset_taxa(all.PS.l[[i]],phylum=="Apicomplexa"), silent=TRUE)
    try(get_taxa_unique(p, "family"), silent=TRUE)
    if (exists("p")) {
        a <- get_taxa_unique(p, "family")
        print(paste(i, "- ", names(all.PS.l[i]), ": ", length(a), sep=""))
        print(a)
}
    rm(p)
}



##### and how many amplicons have eimeria?
for (i in 1:48) {
#    print(names(all.PS.l)[i])
    try(p <- subset_taxa(all.PS.l[[i]],genus=="Eimeria"), silent=TRUE)
    try(get_taxa_unique(p, "genus"), silent=TRUE)
    if (exists("p")) {
        a <- get_taxa_unique(p, "genus")
        print(paste(i, "- ", names(all.PS.l[i]), ": ", nrow(p@tax_table), sep=""))
        print(a)
}
    rm(p)
}

### and what happens when we filter?
for (i in 1:48) {
#    print(names(all.PS.l)[i])
    try(p <- subset_taxa(f.all.l[[i]],genus=="Eimeria"), silent=TRUE)
    try(get_taxa_unique(p, "genus"), silent=TRUE)
    if (exists("p")) {
        a <- get_taxa_unique(p, "genus")
        print(paste(i, "- ", names(f.all.l[i]), ": ", nrow(p@tax_table), sep=""))
        print(a)
}
    rm(p)
}


## how many eimeria ASV reads do we have in single amplicon before and after filtering
subset_taxa(sin.PS, genus=="Eimeria")
# sanity check
subset_taxa(sin.PS18S, genus=="Eimeria")
# after filtering
subset_taxa(f.sin, genus=="Eimeria")
#sanity check
subset_taxa(f.sin18, genus=="Eimeria")

## OK, now we want all the Eimeria sequences
Eim <- subset_taxa(T.all.l, genus%in%"Eimeria")
Eim2 <- subset_taxa(T.sin18, genus%in%"Eimeria")
#Eim.slv <- subset_taxa(T.sin18.slv, Family%in%"Eimeriorina")
Eim_nf <- subset_taxa(all.PS, genus%in%"Eimeria")
Eim2_nf <- subset_taxa(sin.PS18S, genus%in%"Eimeria")
#### sensitivity, specificity, positive predicted value and negative predictive value

#### for MA
# GC- Eim -
summary(sample_sums(Eim_nf@otu_table)==0&Eim_nf@sam_data$Genome_copies_gFaeces==0)
# GC+ Eim -
summary(sample_sums(Eim_nf@otu_table)==0&Eim_nf@sam_data$Genome_copies_gFaeces>0)
# GC+ Eim+
summary(sample_sums(Eim_nf@otu_table)>0&Eim_nf@sam_data$Genome_copies_gFaeces>0)
# GC- Eim+
summary(sample_sums(Eim_nf@otu_table)>0&Eim_nf@sam_data$Genome_copies_gFaeces==0)
# GC+ GC-
summary(Eim_nf@sam_data$Genome_copies_gFaeces==0)
# Eim + Eim -
summary(sample_sums(Eim_nf@otu_table)==0)

e.s

e.s <- matrix(c(53, 3, 50, 129), ncol=2, byrow=TRUE)
colnames(e.s) <- c("Eim-","Eim+")
rownames(e.s) <- c("GC-", "GC+")
margin1 <- margin.table(e.s, margin=1)
margin2 <- margin.table(e.s, margin=2)

#sensitivity
e.s[2,2]/margin1[2]*100
#specificity
e.s[1,1]/margin1[1]*100
#ppv
e.s[2,2]/margin2[2]*100
#npv
e.s[1,1]/margin2[1]*100

#### for SA
# GC- Eim -
summary(sample_sums(Eim2@otu_table)==0&Eim2@sam_data$Genome_copies_gFaeces==0)
# GC+ Eim -
summary(sample_sums(Eim2@otu_table)==0&Eim2@sam_data$Genome_copies_gFaeces>0)
# GC+ Eim+
summary(sample_sums(Eim2@otu_table)>0&Eim2@sam_data$Genome_copies_gFaeces>0)
# GC- Eim+
summary(sample_sums(Eim2@otu_table)>0&Eim2@sam_data$Genome_copies_gFaeces==0)
# GC+ GC-
summary(Eim2@sam_data$Genome_copies_gFaeces==0)
# Eim + Eim -
summary(sample_sums(Eim2@otu_table)==0)

e.sa <- matrix(c(23, 14, 10 152), ncol=2, byrow=TRUE)
colnames(e.sa) <- c("Eim-","Eim+")
rownames(e.sa) <- c("GC-", "GC+")
margin1 <- margin.table(e.sa, margin=1)
margin2 <- margin.table(e.sa, margin=2)

#sensitivity
e.sa[2,2]/margin1[2]*100
#specificity
e.sa[1,1]/margin1[1]*100
#ppv
e.sa[2,2]/margin2[2]*100
#npv
e.sa[1,1]/margin2[1]*100

# now plotting and doing correlation analysis and comparisons
## for "pooled" MA
Plotting_cor(ps=all.PS, "MA", dir="fig/MA/")
# Seqf, TSS, RLE, CLR, ACS
p <- c(0.0208, 0.0134, 0.1037, 0.00001, 0.00001)
p.adjust(p, method="BH")

#only for wang
Plotting_cor(ps=all.PSwang, "MA_wang", dir="fig/MA/")
# seq-f,tss,rle,clr,acs
p <- c(0.00001, 0.0138, 0.0880, 0.0030, 0.00001)
p.adjust(p, method="BH")

## for single amplicon
Plotting_cor(ps=sin.PS18S, "SA", dir="fig/SA/")
# seq-f,tss,rle,clr,acs
p <- c(0.0117, 0.1045, 0.00001, 0.00001, 0.0024)
p.adjust(p, method="BH")

# for MA but individually filtered
Plotting_cor_MA.l(ps=all.PS, ps.f=f.all.lp, "MA_individually_filtered", dir="fig/MA/")

# seq-f,tss,rle,clr,acs
p <- c(0.00001, 0.0146, 0.1003, 0.00001, 0.00001)
p.adjust(p, method="BH")

##################################################################
### OK, so let's remove food
# First we prepare the datasets for MA
plant <- subset_taxa(f.all.lp, !phylum%in%"Streptophyta")
plant <- subPS(plant)
Mus <- subset_taxa(f.all.lp, !phylum%in%"Chordata")
Mus <- subPS(Mus)
worms <- subset_taxa(f.all.lp, !phylum%in%"Nematoda")
worms <- subPS(worms)
TSS <- subPS(f.all.lp)
PlantMus <-  subset_taxa(f.all.lp, !(phylum%in%"Streptophyta"|phylum%in%"Chordata"))
PlantMus <- subPS(PlantMus)

## Now for SA
plant18 <- subset_taxa(f.sin18, !phylum%in%"Streptophyta")
plant18 <- subPS(plant18)
#No Mus in fPS18S
worms18 <- subset_taxa(f.sin18, !phylum%in%"Nematoda")
worms18 <- subPS(worms18)
TSS18 <- subPS(f.sin18)
PlantMusWorms18 <-  subset_taxa(f.sin18, !(phylum%in%"Streptophyta"| phylum%in%"Nematoda" | phylum%in%"Chordata"))
PlantMusWorms18 <- subPS(PlantMusWorms18)

##########################################################
###########################################################
# plotting TSS correlation
a <- p_tss(TSS, "a)", "MA-TSS")
b <- p_tss(plant, "b)", "MA-TSS: -Streptophyta")
c <- p_tss(Mus, "c)", "MA-TSS: -Chordata")
d <- p_tss(worms, "d)", "MA-TSS: -Nematoda")
#e <- p_tss(PlantMusWorms, "e)", "MA no plants, host or nematodes")
e <- p_tss(PlantMus, "e)", "MA no plants or host")

plot_grid(a,b,c,d,e) -> p_cor

MA.a <- lm(data=TSS,logGC~logA)
MA.b <- lm(data=plant,logGC~logA)
MA.c <- lm(data=Mus,logGC~logA)
MA.d <- lm(data=worms,logGC~logA)
MA.e <- lm(data=PlantMus,logGC~logA)

ggplot2::ggsave(file="fig/MA/Biological_rem_MA.pdf", p_cor, width = 15, height = 8, dpi = 600)
ggplot2::ggsave(file="fig/MA/Biological_rem_MA.png", p_cor, width = 15, height = 8, dpi = 600)

a1 <- p_tss(TSS18, "a)", "SA-TSS")
b1 <- p_tss(plant18, "b)", "SA-TSS: -Streptophyta")
#c1 <- p_tss(Mus18, "c)", "SA no host")
c1 <- p_tss(worms18, "c)", "SA-TSS: -Nematoda")
#e1 <- p_tss(PlantMusWorms18, "d)", "SA no plants or nematodes")

SA.a <- lm(data=TSS18,logGC~logA)
SA.b <- lm(data=plant18,logGC~logA)
SA.d <- lm(data=worms18,logGC~logA)

plot_grid(a1,b1,c1) -> p_cor1
ggplot2::ggsave(file="fig/SA/Biological_rem_SA.pdf", p_cor1, width = 15, height = 10, dpi = 600)
ggplot2::ggsave(file="fig/SA/Biological_rem_SA.png", p_cor1, width = 15, height = 10, dpi = 600)

##################################################################
######### pearson tests
# for MA
cor.test(TSS$logGC, TSS$logA, method="pearson")
cor.test(plant$logGC, plant$logA, method="pearson")
cor.test(worms$logGC, worms$logA, method="pearson")
cor.test(Mus$logGC, Mus$logA, method="pearson")
cor.test(PlantMus$logGC, PlantMus$logA, method="pearson")
# for SA
cor.test(TSS18$logGC, TSS18$logA, method="pearson")
cor.test(plant18$logGC, plant18$logA, method="pearson")
cor.test(worms18$logGC, worms18$logA, method="pearson")
#cor.test(PlantMusWorms18$logGC, PlantMusWorms18$logA, method="pearson")

## Now we need to pool into the same dataset for the cocor function
# First for MA
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
# We need to add Eimeria sums from Seq together with the above data frame
PSeimf <-subset_taxa(all.PS, genus%in%"Eimeria")
df <- data.frame(sample_sums(otu_table(PSeimf)))
df$labels <- rownames(df)
names(df) <- c("Eimeriidae", "labels")
df$logEimeriidae <- log(df$Eimeriidae)
df$logEimeriidae
tss.df <- tss.df[,c("labels", "logGC", "logA", "logA_plant", "logA_worms", "logA_Mus", "logA_PlantMus")]
df <- merge(df, tss.df, by="labels")

# cool, now we test differences of correlations with the default Seq
#cocor(~logGC + logEimeriidae | logGC + logA, data = df,
#            test = c("hittner2003", "zou2007"))
cocor(~logGC +  logEimeriidae| logGC + logA_Mus, data = df,
            test = c("hittner2003", "zou2007"))

cocor(~logGC +  logEimeriidae| logGC + logA_plant, data = df,
            test = c("hittner2003", "zou2007"))

cocor(~logGC +  logEimeriidae| logGC + logA_worms, data = df,
      test = c("hittner2003", "zou2007"))

cocor(~logGC +  logEimeriidae| logGC + logA_PlantMus, data = df,
            test = c("hittner2003", "zou2007"))

# Adjust for multiple testing
## those are the corresponding p values from the other correlations above plus the sub TSS).
# seq-f,tss,rle,clr,acs,tss-mus, tss-plant, tss-worms, tss-plant-mus
p <- c(0.00001, 0.0146, 0.1003, 0.00001, 0.00001, 0.0055, 0.1010, 0.0173, 0.461)
round(p.adjust(p, method="BH"), 3)

#### comparison between TSS and TSS-removal
cocor(~logGC +  logA| logGC + logA_plant, data = df,
            test = c("hittner2003", "zou2007"))

cocor(~logGC +  logA| logGC + logA_Mus, data = df,
            test = c("hittner2003", "zou2007"))

cocor(~logGC +  logA| logGC + logA_worms, data = df,
      test = c("hittner2003", "zou2007"))

cocor(~logGC +  logA| logGC + logA_PlantMus, data = df,
            test = c("hittner2003", "zou2007"))

p <- c(0.0012, 0.0163, 0.7727, 0.0529)
round(p.adjust(p, method="BH"), 3)

######################################################
# Then for SA
plant18 <- plant18[,c("labels", "logA")]
names(plant18) <- c("labels", "logA_plant")
tss.df18 <- merge(TSS18, plant18, by="labels")
worms18 <- worms18[,c("labels", "logA")]
names(worms18) <- c("labels", "logA_worms")
tss.df18 <- merge(tss.df18, worms18, by="labels")

# We need to add Eimeria sums from Seq together with the above data frame
PSeimf18 <-subset_taxa(sin.PS18S, genus%in%"Eimeria")
df18 <- data.frame(sample_sums(otu_table(PSeimf18)))
df18$labels <- rownames(df18)
names(df18) <- c("Eimeriidae", "labels")
df18$logEimeriidae <- log(df18$Eimeriidae)
tss.df18 <- tss.df18[,c("labels", "logGC", "logA", "logA_plant", "logA_worms")]
df18 <- merge(df18, tss.df18, by="labels", all=TRUE)

################## comparing to TSS
cocor(~logGC + logA | logGC + logA_plant, data = df18,
            test = c("hittner2003", "zou2007"))

cocor(~logGC + logA | logGC + logA_worms, data = df18,
            test = c("hittner2003", "zou2007"))

p <- c(0.00001, 0.0225)
round(p.adjust(p, method="BH"),3)

## Comparing to Seq
cocor(~logGC + logEimeriidae | logGC + logA_plant, data = df18,
            test = c("hittner2003", "zou2007"))

cocor(~logGC +  logEimeriidae | logGC + logA_worms, data = df18,
            test = c("hittner2003", "zou2007"))

# seq-f,tss,rle,clr,acs, tss-plant, tss-worms
p <- c(0.0117, 0.1045, 0.00001, 0.00001, 0.0024, 0.0597, 0.0289)
round(p.adjust(p, method="BH"),3)

##################################################



##################### how many amplicons have mus?
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

### let's save the mus sequences
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



