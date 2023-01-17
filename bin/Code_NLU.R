# this was taken from 4_MA_Correlations.R



#######################################################
##################################################################
### OK, so let's remove food
# First we prepare the datasets for MA

get_taxa_unique(f.all.lp.slv, "Phylum")

plant <- subset_taxa(f.all.lp.slv, !Phylum%in%c("p_Anthophyta", "c__Monocotyledonae", "p__Charophyta", "p__Phragmoplastophyta"))
plant <- subPS(plant)
Mus <- subset_taxa(f.all.lp.slv, !Phylum%in%"p__Vertebrata")
Mus <- subPS(Mus)
worms <- subset_taxa(f.all.lp.slv, !Phylum%in%"p__Nematozoa")
worms <- subPS(worms)
TSS <- subPS(f.all.lp.slv)
PlantMus <-  subset_taxa(f.all.lp.slv, !Phylum%in%c("p_Anthophyta", "c__Monocotyledonae", "p__Charophyta", "p__Phragmoplastophyta", "p__Vertebrata"))
PlantMus <- subPS(PlantMus)

## Now for SA
plant18 <- subset_taxa(f.sin18.slv, !Phylum%in%c("p_Anthophyta", "c__Monocotyledonae", "p__Charophyta", "p__Phragmoplastophyta"))
plant18 <- subPS(plant18)
Mus18 <- subset_taxa(f.sin18.slv, !Phylum%in%"p__Vertebrata")
Mus18 <- subPS(Mus)
#No Mus in fPS18S
worms18 <- subset_taxa(f.sin18.slv, !Phylum%in%"p__Nematozoa")
worms18 <- subPS(worms18)
TSS18 <- subPS(f.sin18.slv)
PlantMus18 <-  subset_taxa(f.sin18.slv, !Phylum%in%c("p_Anthophyta", "c__Monocotyledonae", "p__Charophyta", "p__Phragmoplastophyta", "p__Vertebrata"))
PlantMus18 <- subPS(PlantMus18)

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

ggplot2::ggsave(file="fig/MA/Biological_rem_MA.pdf", p_cor, width = 15, height = 8, dpi = 600)
ggplot2::ggsave(file="fig/MA/Biological_rem_MA.png", p_cor, width = 15, height = 8, dpi = 600)

a1 <- p_tss(TSS18, "a)", "SA-TSS")
b1 <- p_tss(plant18, "b)", "SA-TSS: -diet")
c1 <- p_tss(Mus18, "c)", "SA no host")
d1 <- p_tss(worms18, "d)", "SA-TSS: -nematodes")
e1 <- p_tss(PlantMus18, "d)", "SA no diet or host")

plot_grid(a1,b1,c1, d1, e1 ) -> p_cor1
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

# cool, now we test differences of correlations with the default Seq
#cocor(~logGC + logEimeriidae | logGC + logA, data = df,
#            test = c("hittner2003", "zou2007"))
cocor(~logGC +  logA| logGC + logA_Mus, data = tss.df,
            test = c("hittner2003", "zou2007"))

cocor(~logGC +  logA| logGC + logA_plant, data = tss.df,
            test = c("hittner2003", "zou2007"))

cocor(~logGC +  logA| logGC + logA_worms, data = tss.df,
      test = c("hittner2003", "zou2007"))

## human, plant, worms
p <- c(0.00001, 0.00001, 0.5833)
round(p.adjust(p, method="BH"), 3)


######################################################
# Then for SA
plant18 <- plant18[,c("labels", "logA")]
names(plant18) <- c("labels", "logA_plant")
tss.df18 <- merge(TSS18, plant18, by="labels")
worms18 <- worms18[,c("labels", "logA")]
names(worms18) <- c("labels", "logA_worms")
tss.df18 <- merge(tss.df18, worms18, by="labels")
#tss.df18 <- tss.df18[,c("labels", "logGC", "logA", "logA_plant", "logA_worms")]

################## comparing to TSS
cocor(~logGC + logA | logGC + logA_plant, data = tss.df18,
            test = c("hittner2003", "zou2007"))

cocor(~logGC + logA | logGC + logA_worms, data = tss.df18,
      test = c("hittner2003", "zou2007"))

cocor(~logGC + logA | logGC + logA_plant, data = tss.df18,
            test = c("hittner2003", "zou2007"))

p <- c(0.0001, 0.05, 0.00001)
round(p.adjust(p, method="BH"),3)

##################################################

get_taxa_unique(f.sin18.slv, "Genus")
get_taxa_unique(f.sin18, "genus")

get_taxa_unique(T.all.slv, "Genus")
get_taxa_unique(T.all, "genus")


##################### how many amplicons have mammal asvs?
################################### Mus sequences
for (i in 1:48) {
#    print(names(all.PS.l)[i])
    try(p <- subset_taxa(f.all.l.slv[[i]],Genus=="g__Mammalia"), silent=TRUE)
    if (exists("p")) {
        try(a <- rownames(p@tax_table), silent=TRUE)
        print(paste(i, "- ", names(f.all.l[[i]]), ": ", length(a), sep=""))
        print(a)
}
    rm(p)
}

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


f.all.l

f.all.l.slv[[12]]@tax_table

f.all.l[[12]]@tax_table

f.all.l.slv[[33]]@tax_table

f.all.l[[33]]@tax_table

#same MUS sequence in MA and SA
#a==rownames(subset_taxa(f.sin18.slv,Genus=="g__Mammalia")@tax_table)

### let's save the mus sequences
#library(ShortRead)
#writeFasta(DNAStringSet(a), "tmp/Mus_ASV.fasta")

################### stop here

#######################
# how does host DNA change with infection
# no Mus reads in SA
Mus <- subset_taxa(T.all.slv, Genus%in%"g__Mammalia")
MusTSS <- subset_taxa(allTSS.slv, Genus%in%"g__Mammalia")
m <- psmelt(Mus)
mTSS <- psmelt(MusTSS)

Mus18 <- subset_taxa(T.sin18.slv, Genus%in%"g__Mammalia")
MusTSS18 <- subset_taxa(sin18TSS.slv, Genus%in%"g__Mammalia")
m18 <- psmelt(Mus18)
mTSS18 <- psmelt(MusTSS18)


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

mp18.acs <- ggplot(m18, aes(dpi, Abundance, color=EH_ID))+
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


m.tss <- ggplot(mTSS, aes(dpi, Abundance, color=EH_ID))+
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

m18.tss <- ggplot(mTSS18, aes(dpi, Abundance, color=EH_ID))+
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


ggplot(m, aes(dpi, weightloss, color=EH_ID))+
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


mp.tss

mp.acs

m18.tss
mp18.acs

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
#    annotate(geom="text", x=0.7, y=5, label="Spearman rho=0.36, p<0.001")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          legend.position = "top",
          axis.line = element_line(colour = "black"))+
    guides(fill=guide_legend(nrow=1, byrow=TRUE))

m.host18 <- ggplot(m18, aes(Abundance, log(1+Genome_copies_ngDNA)))+
    geom_point(aes(fill=dpi), shape=21, size=4, alpha=0.6)+
#    scale_color_brewer(palette="Accent")+
    scale_fill_brewer(palette="Spectral")+
#    geom_smooth(method = "lm", se=FALSE, na.rm=TRUE) +
#    stat_poly_eq(aes(label = paste(..eq.label..,
#                                   ..rr.label..,
#                                   sep = "~~~")),
#                 parse = TRUE, label.x=0.85, label.y=0.23) + 
#    annotate(geom="text", x=0.4, y=5, label="Spearman rho=0.37, p<0.001")+
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



