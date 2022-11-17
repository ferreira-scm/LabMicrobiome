
library(relaimpo)
library(lme4)
library(MASS)
library(DescTools)
library(cocor)
library(reshape2)
library(microbiome)
library(phyloseq)

source("bin/4_MA_SA_filtering.R")

############################################
##############################################

#### sensitivity, specificity, positive predicted value and negative predictive value
############# sensitivity and specificity

sensit <- function(Eim_nf){
# GC- Eim -
tneg <- summary(sample_sums(Eim_nf@otu_table)==0&Eim_nf@sam_data$Genome_copies_gFaeces==0)[[3]]

# GC+ Eim -
fneg <- summary(sample_sums(Eim_nf@otu_table)==0&Eim_nf@sam_data$Genome_copies_gFaeces>0)[[3]]

# GC+ Eim+
tpos <- summary(sample_sums(Eim_nf@otu_table)>0&Eim_nf@sam_data$Genome_copies_gFaeces>0)[[3]]

# GC- Eim+
fpos <- summary(sample_sums(Eim_nf@otu_table)>0&Eim_nf@sam_data$Genome_copies_gFaeces==0)[[3]]

# GC+ GC-
gcpos <- summary(Eim_nf@sam_data$Genome_copies_gFaeces==0)[[2]]
gcneg <- summary(Eim_nf@sam_data$Genome_copies_gFaeces==0)[[3]]

# Eim + Eim -
epos <- summary(sample_sums(Eim_nf@otu_table)==0)[[2]]
eneg <- summary(sample_sums(Eim_nf@otu_table)==0)[[3]]

e.s <- matrix(as.numeric(c(tneg, fpos, fneg, tpos)), ncol=2, byrow=TRUE)
colnames(e.s) <- c("Eim-","Eim+")
rownames(e.s) <- c("GC-", "GC+")
margin1 <- margin.table(e.s, margin=1)
margin2 <- margin.table(e.s, margin=2)

#sensitivity
print("Sensitivity:")
print(e.s[2,2]/margin1[2]*100)
#specificity
print("Specificity:")
print(e.s[1,1]/margin1[1]*100)
#ppv
print("positive predictive value")
print(e.s[2,2]/margin2[2]*100)
#npv
print("negative predictive value")
print(e.s[1,1]/margin2[1]*100)
}

#### for MA
sensit(Eim_nf)
sensit(Eim)

MA.asv1 <- prune_taxa(rownames(tax_table(Eim2))[1], Eim)
MA.asv12 <- prune_taxa(rownames(tax_table(Eim2))[c(1,2)], Eim)
MA.asv13 <- prune_taxa(rownames(tax_table(Eim2))[c(1,3)], Eim)
MA.asv2 <- prune_taxa(rownames(tax_table(Eim2))[2], Eim)
MA.asv3 <- prune_taxa(rownames(tax_table(Eim2))[3], Eim)


sensit(MA.asv1)

sensit(MA.asv13)

sensit(MA.asv12)

sensit(MA.asv2)

sensit(MA.asv3)

#### for SA
sensit(Eim2_nf)

sensit(Eim2)

SA.asv1 <- prune_taxa(rownames(tax_table(Eim2))[1], Eim2)
SA.asv12 <- prune_taxa(rownames(tax_table(Eim2))[c(1,2)], Eim2)
SA.asv123 <- prune_taxa(rownames(tax_table(Eim2))[c(1,2,3)], Eim2)
SA.asv13 <- prune_taxa(rownames(tax_table(Eim2))[c(1,3)], Eim2)
SA.asv2 <- prune_taxa(rownames(tax_table(Eim2))[2], Eim2)
SA.asv3 <- prune_taxa(rownames(tax_table(Eim2))[3], Eim2)
SA.asv4 <- prune_taxa(rownames(tax_table(Eim2))[4], Eim2)
SA.asv5 <- prune_taxa(rownames(tax_table(Eim2))[5], Eim2)

sensit(SA.asv1)

sensit(SA.asv2)

sensit(SA.asv12)

sensit(SA.asv13)

sensit(SA.asv123)

sensit(SA.asv3)

sensit(SA.asv4)

sensit(SA.asv5)

# now plotting and doing correlation analysis and comparisons
## for "pooled" MA
#Plotting_cor(ps=all.PS, "MA", dir="fig/MA/")
# Seqf, TSS, RLE, CLR, ACS
#p <- c(0.0208, 0.0134, 0.1037, 0.00001, 0.00001)
#p.adjust(p, method="BH")

## for single amplicon
Plotting_cor_MA.l(ps=sin.PS18S, f.sin18, "SA", dir="fig/SA/")

# seq-f,tss,rle,clr,acs
p <- c(0.0406, 0.1725, 0.00001, 0.00001, 0.6218, 0.9807)
p.adjust(p, method="BH")

# for MA but individually filtered
Plotting_cor_MA.l(ps=all.PS, ps.f=f.all.lp, "MA_individually_filtered", dir="fig/MA/")

# seq-f,tss,rle,clr,acs

p <- c(0.0001, 0.0052, 0.0037, 0.00001, 0.0006, 0.1941)
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
#cor.test(PlantMus$logGC, PlantMus$logA, method="pearson")
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

p <- c(0.0029, 0.0393, 0.0077)
round(p.adjust(p, method="BH"), 3)


# Adjust for multiple testing
## those are the corresponding p values from the other correlations above plus the sub TSS).
# seq-f,tss,rle,clr,acs,tss-mus, tss-plant, tss-worms
p <- c(0.0001, 0.0052, 0.0037, 0.00001, 0.0006, 0.0023, 0.0337, 0.0064)
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

p <- c(0.0070, 0.1331, 0.781)
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

p <- c(0.0001, 0.0107)
round(p.adjust(p, method="BH"),3)

## Comparing to Seq
cocor(~logGC + logEimeriidae | logGC + logA_plant, data = df18,
            test = c("hittner2003", "zou2007"))

cocor(~logGC +  logEimeriidae | logGC + logA_worms, data = df18,
            test = c("hittner2003", "zou2007"))

# seq-f,tss,rle,clr,acs, tss-plant, tss-worms
p <- c(0.0406, 0.1725, 0.00001, 0.00001, 0.6218, 0.0597, 0.0289)
round(p.adjust(p, method="BH"),3)

##################################################

##################### how many amplicons have mus?
################################### Mus sequences
for (i in 1:48) {
#    print(names(all.PS.l)[i])
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
Mus <- subset_taxa(T.all, genus%in%"Mus")
MusTSS <- subset_taxa(allTSS, genus%in%"Mus")

Mus <-aggregate_taxa(Mus, level="genus")
MusTSS <-aggregate_taxa(MusTSS, level="genus")

m <- psmelt(Mus)
mTSS <- psmelt(MusTSS)


nb <- length(levels(as.factor(m$EH_ID)))

coul <- colorRampPalette(brewer.pal(8, "Accent"))(nb)

mycol <- coul[as.numeric(as.factor(m$EH_ID))]

mp.acs <- ggplot(m, aes(dpi, Abundance, color=EH_ID))+
    geom_point(aes(fill=EH_ID), shape=21, size=4, position=position_jitter(0.1), alpha=0.6)+
    scale_color_manual(values=coul)+
        scale_fill_manual(values=coul)+           
    geom_line(aes(group=EH_ID))+
    labs(x="Days post infection", y="Host reads (per ngDNA)")+
    labs(tag= "a)")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
           legend.position = "none",
          axis.line = element_line(colour = "black"))

mp.tss <- ggplot(mTSS, aes(dpi, Abundance, color=EH_ID))+
    geom_point(aes(fill=EH_ID), shape=21, size=4, position=position_jitter(0.1), alpha=0.6)+
    scale_color_manual(values=coul)+
        scale_fill_manual(values=coul)+           
    geom_line(aes(group=EH_ID))+
    labs(x="Days post infection", y="Host reads (relative abundance)")+
    labs(tag= "a)")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
           legend.position = "none",
          axis.line = element_line(colour = "black"))


mp.tss

mp.acs

# Are host abundance and Eimeria associated
cor.test(m$Abundance, m$Genome_copies_ngDNA, method="spearman")

cor.test(m$Abundance, m$OPG, method="spearman")

cor.test(m$Abundance, m$weightloss, method="spearman")

cor.test(mTSS$Abundance, mTSS$Genome_copies_ngDNA, method="spearman")
#cor.test(m.clr$Abundance, log(1+m.clr$OPG), method="spearman")

m.host <- ggplot(m, aes(log(1+Abundance), log(1+Genome_copies_ngDNA)))+
    geom_point(aes(fill=dpi), shape=21, size=4, alpha=0.6)+
    scale_fill_brewer(palette="Spectral")+           
    labs(x="Host reads/ngDNA, log(+1)", y="Genome copies/ngDNA log(1+)")+
    labs(tag= "b)")+
#    geom_smooth(method = "lm", se=FALSE, na.rm=TRUE) +
#    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE,  label.x=0.1, label.y=0.9) +  
    annotate(geom="text", x=0.7, y=5, label="Spearman rho=0.36, p<0.001")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          legend.position = "top",
          axis.line = element_line(colour = "black"))+
    guides(fill=guide_legend(nrow=1, byrow=TRUE))

m.tss.host <- ggplot(mTSS, aes(Abundance, log(1+Genome_copies_ngDNA)))+
    geom_point(aes(fill=dpi), shape=21, size=4, alpha=0.6)+
#    scale_color_brewer(palette="Accent")+
    scale_fill_brewer(palette="Spectral")+
#    geom_smooth(method = "lm", se=FALSE, na.rm=TRUE) +
#    stat_poly_eq(aes(label = paste(..eq.label..,
#                                   ..rr.label..,
#                                   sep = "~~~")),
#                 parse = TRUE, label.x=0.85, label.y=0.23) + 
    annotate(geom="text", x=0.4, y=5, label="Spearman rho=0.37, p<0.001")+
#    geom_line(aes(group=EH_ID), alpha=0.3)+
    labs(x="Host reads relative abundace", y="Genome copies/ngDNAg log(1+)")+
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
plot_grid(mp.tss,m.tss.host, nrow=2) -> m.tssCor

ggplot2::ggsave(file="fig/Mus_Eimeria.tss.pdf", m.tssCor, width = 8, height = 8, dpi = 600)
ggplot2::ggsave(file="fig/Mus_Eimeria.tss.png", m.tssCor, width = 8, height = 8, dpi = 600)

ggplot2::ggsave(file="fig/Mus_Eimeria.acs.pdf", m.acsCor, width = 8, height = 8, dpi = 600)
ggplot2::ggsave(file="fig/Mus_Eimeria.acs.png", m.acsCor, width = 8, height = 8, dpi = 600)


# ordination
# host reads, total dna, opg, Genome_copies
# can't have zeros in all columns

m.ord <- m[,c("Genome_copies_gFaeces", "Abundance", "OPG")][rowSums(m[,c("Genome_copies_gFaeces", "Abundance", "OPG")])>0,]

m.tss.ord <- mTSS[,c("Genome_copies_gFaeces", "Abundance", "OPG", "DNA_g_feces")][rowSums(m[,c("Genome_copies_gFaeces", "Abundance", "OPG", "DNA_g_feces")])>0,]

library(vegan)

m.mds <- metaMDS(m.ord)

m.tss.euc <- vegdist(m.tss.ord, method="euclidean")

m.bray <- vegdist(m.ord, "bray")

m.tss.bray <- vegdist(m.tss.ord, "bray")


dpi <- as.factor(m$dpi[as.numeric(rownames(m.ord))])
dpi.tss <- as.factor(mTSS$dpi[as.numeric(rownames(m.tss.ord))])

weighloss <- m$weightloss[as.numeric(rownames(m.ord))]

weighloss.tss <- mTSS$weightloss[as.numeric(rownames(m.tss.ord))]

adonis2(m.bray~dpi+weighloss, method="bray", permutations=10000)

adonis2(m.tss.bray~dpi.tss+weighloss.tss, method="euclidean", permutations=10000)

#sanity check
rownames(m)==rownames(m.ord)


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

mTSS$logGC <- log(1+mTSS$Genome_copies_ngDNA)

m$logA <- log(1+m$Abundance)
m$logGC <- log(1+m$Genome_copies_ngDNA)

##TSS
mmus <- lm(Abundance~logGC + Total_DNA, data=mTSS)
summary(mmus)

mwl <- lm(weightloss~Abundance+Genome_copies_ngDNA+ Total_DNA, data=mTSS)
summary(mwl)
anova(mwl)

histogram(m$weightloss)

cm1 <- lmer(weightloss~Abundance+Genome_copies_gFaeces +dpi+(1|EH_ID), data=mTSS)
cm2 <- lmer(weightloss~Abundance+Genome_copies_gFaeces +(1|EH_ID), data=mTSS)
cm3 <- lmer(weightloss~Abundance+dpi+(1|EH_ID), data=mTSS)
cm4 <- lmer(weightloss~Genome_copies_gFaeces +dpi+(1|EH_ID), data=mTSS)

## ACS
#mwl.acs <- lm(weightloss~logA+logGC, data=m)
#plot(mwl.acs)




mwl.acs <- lm(weightloss~Abundance+Genome_copies_ngDNA, data=m)

calc.relimp(mwl.acs)

mmwl.acs1 <- lmer(weightloss~Abundance+Genome_copies_ngDNA +dpi+(1|EH_ID), data=m)

mmwl.acs2 <- lmer(weightloss~Abundance+Genome_copies_ngDNA+(1|EH_ID), data=m)

mmwl.acs3 <- lmer(weightloss~Abundance +dpi+(1|EH_ID), data=m)

mmwl.acs4 <- lmer(weightloss~Genome_copies_ngDNA +dpi+(1|EH_ID), data=m)

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


# how does Mus DNA change with infection
cor.test(m$logA, log(1+m$weightloss), method="spearman")



