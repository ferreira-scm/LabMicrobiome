#!/usr/bin/Rscript

#library(reshape)
library(DECIPHER)
library(phangorn)
library(ShortRead)
library(RColorBrewer)
library("seqinr")
library(treeio)
library(ggtree)
library(magrittr)
library(ape)
library(BiocManager)
library(lmerTest)
library(lme4)

source("bin/4_MA_SA_filtering.R")

# Do ASV match beween MA and SA?
#unfiltered
rownames(Eim2@tax_table) %in% rownames(Eim@tax_table)
#filtered
colnames(Eim@otu_table)[1] == colnames(Eim2@otu_table)[1]
colnames(Eim@otu_table)[2] == colnames(Eim2@otu_table)[2]

#plotting individual ASV's
MA.e <- psmelt(Eim)
SA.e <- psmelt(Eim2)

Eim2.g <- tax_glom(Eim2, taxrank="Genus")
SA.e.g <- psmelt(Eim2.g)
Eim.g <- tax_glom(Eim, taxrank="Genus")
MA.e.g <- psmelt(Eim.g)

SA.e$ASV <- "ASV"
SA.e$ASV[which(SA.e$OTU==colnames(Eim2@otu_table)[5])] <- "ASV5"
SA.e$ASV[which(SA.e$OTU==colnames(Eim2@otu_table)[4])] <- "ASV4"
SA.e$ASV[which(SA.e$OTU==colnames(Eim2@otu_table)[3])] <- "ASV3"
SA.e$ASV[which(SA.e$OTU==colnames(Eim2@otu_table)[2])] <- "ASV2"
SA.e$ASV[which(SA.e$OTU==colnames(Eim2@otu_table)[1])] <- "ASV1"

MA.e$ASV <- "ASV"
MA.e$ASV[which(MA.e$OTU==colnames(Eim@otu_table)[3])] <- "ASV3"
MA.e$ASV[which(MA.e$OTU==colnames(Eim@otu_table)[2])] <- "ASV2"
MA.e$ASV[which(MA.e$OTU==colnames(Eim@otu_table)[1])] <- "ASV1"

SA.e5 <- SA.e[which(SA.e$OTU==colnames(Eim2@otu_table)[5]),]
SA.e4 <- SA.e[which(SA.e$OTU==colnames(Eim2@otu_table)[4]),]
SA.e3 <- SA.e[which(SA.e$OTU==colnames(Eim2@otu_table)[3]),]
SA.e2 <- SA.e[which(SA.e$OTU==colnames(Eim2@otu_table)[2]),]
SA.e1 <- SA.e[which(SA.e$OTU==colnames(Eim2@otu_table)[1]),]

## dataset with asv1 and 2

MA.e3 <- MA.e[which(MA.e$OTU==colnames(Eim@otu_table)[3]),]
MA.e2 <- MA.e[which(MA.e$OTU==colnames(Eim@otu_table)[2]),]
MA.e1 <- MA.e[which(MA.e$OTU==colnames(Eim@otu_table)[1]),]

# No GC and ASV abundance zeros
SA.e.0 <- SA.e[SA.e$Genome_copies_ngDNA>0,]
SA.e.0 <- SA.e.0[SA.e.0$Abundance>0,]
SA.e.0 <-SA.e.0[!is.na(SA.e.0$Abundance),]

MA.e.0 <- MA.e[MA.e$Genome_copies_ngDNA>0,]
MA.e.0 <- MA.e.0[MA.e.0$Abundance>0,]
MA.e.0 <-MA.e.0[!is.na(MA.e.0$Abundance),]

## Let's do the correlations first
SA.ASV <- SA.e[SA.e$Genome_copies_ngDNA>0,]
SA.e.g0 <- SA.e.g[SA.e.g$Genome_copies_ngDNA>0,]
SA.e.g0 <- SA.e.g0[SA.e.g0$Abundance>0,]
SA.e.g0 <- SA.e.g0[!is.na(SA.e.g0$Genome_copies_ngDNA),]

MA.e.g0 <- MA.e.g[MA.e.g$Genome_copies_ngDNA>0,]
MA.e.g0 <- MA.e.g0[MA.e.g0$Abundance>0,]
MA.e.g0 <- MA.e.g0[!is.na(MA.e.g0$Genome_copies_ngDNA),]

# sanity test here
cor.test(log(SA.e.g0$Abundance), log(SA.e.g0$Genome_copies_ngDNA))
cor.test(log(MA.e.g0$Abundance), log(MA.e.g0$Genome_copies_ngDNA))

SA.e10 <- SA.e1[SA.e1$Genome_copies_ngDNA>0,]
SA.e10 <- SA.e10[!is.na(SA.e10$Genome_copies_ngDNA),]
SA.e10 <- SA.e10[SA.e10$Abundance>0,]

SA.e20 <- SA.e2[SA.e2$Genome_copies_ngDNA>0,]
SA.e20 <- SA.e20[!is.na(SA.e10$Genome_copies_ngDNA),]
SA.e20 <- SA.e20[SA.e20$Abundance>0,]

SA.e30 <- SA.e3[SA.e3$Genome_copies_ngDNA>0,]
SA.e30 <- SA.e30[!is.na(SA.e30$Genome_copies_ngDNA),]
SA.e30 <- SA.e30[SA.e30$Abundance>0,]

SA.e40 <- SA.e4[SA.e4$Genome_copies_ngDNA>0,]
SA.e40 <- SA.e40[!is.na(SA.e40$Genome_copies_ngDNA),]
SA.e40 <- SA.e40[SA.e40$Abundance>0,]

### they are correlate pretty well, ASV1 being the best
cor.test(log(SA.e10$Genome_copies_ngDNA), log(SA.e10$Abundance))
cor.test(log(SA.e20$Genome_copies_ngDNA), log(SA.e20$Abundance))
cor.test(log(SA.e30$Genome_copies_ngDNA), log(SA.e30$Abundance))
cor.test(log(SA.e40$Genome_copies_ngDNA), log(SA.e40$Abundance))

### now same for MA
MA.e10 <- MA.e1[MA.e1$Genome_copies_ngDNA>0,]
MA.e10 <- MA.e10[!is.na(MA.e10$Genome_copies_ngDNA),]
MA.e10 <- MA.e10[MA.e10$Abundance>0,]

MA.e20 <- MA.e2[MA.e2$Genome_copies_ngDNA>0,]
MA.e20 <- MA.e20[!is.na(MA.e10$Genome_copies_ngDNA),]
MA.e20 <- MA.e20[MA.e20$Abundance>0,]

## ASV1 is also best for MA
cor.test(log(MA.e10$Genome_copies_ngDNA), log(MA.e10$Abundance))
cor.test(log(MA.e20$Genome_copies_ngDNA), log(MA.e20$Abundance))

####
## Let's do the correlations first
SA.opg <- SA.e.g[SA.e.g$OPG>0,]
SA.opg <- SA.opg[SA.opg$Abundance>0,]

### OPG correlation
cor.test(log(SA.opg$OPG), log(SA.opg$Abundance))

col7 <- c("#F1B6DA", "#C51B7D", "#DE77AE","#01665E", "#35978F", "#80CDC1","#C7EAE5")

OPG_Abundance <- ggplot(SA.opg, aes(y=log(OPG), x=log(Abundance), fill=dpi))+
    geom_point(shape=21, size=4, alpha=0.8)+
    scale_fill_manual(values=col7)+
    ylab("Oocysts/g faeces (log)")+
#    guides(fill=guide_legend(ncol=3, byrow=TRUE))+
    xlab("Eimeria ASV abundance (log)")+
    annotate(geom="text", x=min(log(SA.opg$Abundance)), y=max(log(SA.opg$OPG)), hjust=0.05, label="Pearson rho=0.67, p<0.001", size=3)+
    theme_bw(base_size=10)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.line = element_line(colour = "black"))

## now same but for MA
MA.opg <- MA.e.g[MA.e.g$OPG>0,]
MA.opg <- MA.opg[MA.opg$Abundance>0,]
### OPG correlation
cor.test(log(MA.opg$OPG), log(MA.opg$Abundance))

OPG_Abundance_MA <- ggplot(MA.opg, aes(y=log(OPG), x=log(Abundance), fill=dpi))+
    geom_point(shape=21, size=4, alpha=0.8)+
    scale_fill_manual(values=col7)+
    ylab("Oocysts/g faeces (log)")+
#    guides(fill=guide_legend(ncol=3, byrow=TRUE))+
    xlab("Eimeria ASV abundance (log)")+
    annotate(geom="text", x=min(log(MA.opg$Abundance)), y=max(log(MA.opg$OPG)), hjust=0.05, label="Pearson rho=0.59, p<0.001", size=3)+
    theme_bw(base_size=10)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
           legend.position = "none",
          axis.line = element_line(colour = "black"))

OPG_Abundance_MA

legend <- get_legend(OPG_Abundance_MA+
          guides(fill=guide_legend(nrow=1, byrow=TRUE))+
          theme(legend.position="top"))

OPG_ab <- plot_grid(OPG_Abundance, OPG_Abundance_MA, labels="auto")

OPG_ab_panel <- plot_grid(legend, OPG_ab,  nrow=2, rel_heights=c(0.1, 0.8))



ggplot2::ggsave(file="fig/Eimeria_ASV_OPG.pdf", OPG_ab_panel, width = 8, height = 6, dpi = 300)
ggplot2::ggsave(file="fig/Eimeria_ASV_OPG.png", OPG_ab_panel, width = 8, height = 6, dpi = 300)

# saving object so that we can plot them together with wild
saveRDS(OPG_Abundance, "tmp/OPG_Abundance.R")
# saving object so that we can plot them together with wild
saveRDS(OPG_Abundance_MA, "tmp/OPG_Abundance_MA.R")
saveRDS(OPG_ab_panel, "tmp/OPG_Abundance_MA_panel.R")


#### Checking specificity

mean(SA.e.g$Abundance[SA.e.g$Genome_copies_ngDNA==0], na.rm=TRUE)

SA.e.g$Abundance[SA.e.g$Genome_copies_ngDNA==0]

sd(SA.e.g$Abundance[SA.e.g$Genome_copies_ngDNA==0], na.rm=TRUE)

MA.e$Abundance[MA.e$Genome_copies_ngDNA==0]

mean(MA.e.g$Abundance[MA.e.g$Genome_copies_ngDNA==0], na.rm=TRUE)

sd(MA.e.g$Abundance[MA.e.g$Genome_copies_ngDNA==0], na.rm=TRUE)

mean(SA.e.g$Abundance)
sd(SA.e.g$Abundance)

length(SA.e.g$Abundance)

mean(MA.e.g$Abundance)

sd(MA.e.g$Abundance)
length(MA.e.g$Abundance)


############### Which ASV explains Eimeria genome copies?
################ preparing dataset for regressions
# removeing zeros from qPCR
# can't remove zeros from indidual ASV's since then I end up with a handful of samples
Eim2.0 = phyloseq::prune_samples(Eim2@sam_data$Genome_copies_ngDNA>0, Eim2)
Eim.0 = phyloseq::prune_samples(Eim@sam_data$Genome_copies_ngDNA>0, Eim)

## OK, I will now do a linear model with only postive PCR samples and positive sequencing reads (total Eimeria sums) and include the sum of all Eimeria ASV's.

SA.df <- data.frame(Eim2.0@sam_data$Genome_copies_ngDNA, Eim2.0@otu_table[,1], Eim2.0@otu_table[,2],Eim2.0@otu_table[,3],Eim2.0@otu_table[,4],Eim2.0@sam_data$dpi, Eim2.0@sam_data$EH_ID)
names(SA.df) <- c("Genome_copies", "ASV1", "ASV2", "ASV3", "ASV4", "dpi", "EH_ID")

SA.df$TotalE <- SA.df$ASV1+SA.df$ASV2+SA.df$ASV3+SA.df$ASV4
SA.df <- SA.df[SA.df$TotalE>0,]
nrow(SA.df)

lmm.sa <- lmer(data=SA.df, log(Genome_copies)~log(1+ASV1)+log(1+ASV2)+log(1+ASV3)+log(1+ASV4)+ (1|dpi))
lmm.sa0 <- lmer(data=SA.df, log(Genome_copies)~1+ (1|dpi))

ranova(lmm.sa)
summary(lmm.sa)

# calculating r-squared following Schielzeth et al. 2012
# Extraction of fitted value for the alternative model
Fixed <- fixef(lmm.sa)[2]*getME(lmm.sa, "X")[,2]+fixef(lmm.sa)[3]*getME(lmm.sa, "X")[,3]+fixef(lmm.sa)[4]*getME(lmm.sa, "X")[,4]++fixef(lmm.sa)[5]*getME(lmm.sa, "X")[,5]

# Calculation of the variance in fitted values
VarF <- var(Fixed)
# R2GLMM(m) - marginal R2GLMM
# VarCorr() extracts variance components
# attr(VarCorr(lmer.model),'sc')^2 extracts the residual variance
VarF/(VarF + VarCorr(lmm.sa)$dpi[1] + attr(VarCorr(lmm.sa), "sc")^2)

# R2GLMM(c) - conditional R2GLMM for full model
(VarF + VarCorr(lmm.sa)$dpi[1])/(VarF + VarCorr(lmm.sa)$dpi[1] + (attr(VarCorr(lmm.sa), "sc")^2))

#plot(lmm.sa)
anova(lmm.sa, test="LRT")
anova(lmm.sa, lmm.sa0)
# saving model results
sink("fig/Eimeira_asv_SA.txt")
summary(lmm.sa)
sink()

## plotting estimated random effects for each dpi and its interval estimate

library(merTools)

plotREsim(REsim(lmm.sa))

confint(lmm.sa)


## is using negative binomial better than transforming?
gau.m.0.nb <- glmer.nb(data=SA.df, Genome_copies~ASV1+ASV2+ASV3+ASV4+(1|dpi))
summary(gau.m.0.nb)


#### for multiamplicon
MA.df <- data.frame(Eim.0@sam_data$Genome_copies_ngDNA, Eim.0@otu_table[,1], Eim.0@otu_table[,2],Eim.0@sam_data$dpi)
names(MA.df) <- c("Genome_copies", "ASV1", "ASV2", "dpi")
MA.df$TotalE <- MA.df$ASV1+MA.df$ASV2
MA.df <- MA.df[MA.df$TotalE>0,]
nrow(MA.df)

lmm.ma <- lmer(data=MA.df, log(Genome_copies)~log(1+ASV1)+log(1+ASV2)+(1|dpi))
summary(lmm.ma)
lmm.ma0 <- lmer(data=MA.df, log(Genome_copies)~1+(1|dpi))

sink("fig/Eimeira_asv_MA.txt")
print(summary(lmm.ma))
sink()

# LR tests for significance
anova(lmm.ma, test="LRT")
anova(lmm.ma, lmm.ma0)
ranova(lmm.ma)
# calculating r-squared following Schielzeth et al. 2012
# Extraction of fitted value for the alternative model
mFixed <- fixef(lmm.ma)[2]*getME(lmm.ma, "X")[,2]+fixef(lmm.ma)[3]*getME(lmm.ma, "X")[,3]

# Calculation of the variance in fitted values
mVarF <- var(mFixed)
# R2GLMM(m) - marginal R2GLMM
# VarCorr() extracts variance components
# attr(VarCorr(lmer.model),'sc')^2 extracts the residual variance
mVarF/(mVarF + VarCorr(lmm.ma)$dpi[1] + attr(VarCorr(lmm.ma), "sc")^2)

# R2GLMM(c) - conditional R2GLMM for full model
(mVarF + VarCorr(lmm.ma)$dpi[1])/(mVarF + VarCorr(lmm.ma)$dpi[1] + (attr(VarCorr(lmm.ma), "sc")^2))

## same thing using a package, as sanity check

library(MuMIn)

r.squaredGLMM(lmm.sa)
r.squaredGLMM(lmm.ma)

## is using negative binomial better than transforming?
gau.m.0.ma.nb <- glmer.nb(data=MA.df, Genome_copies~ASV1+ASV2+(1|dpi))
summary(gau.m.0.ma.nb)

### testing Victor's hypothesis for oocysts
gau.o <- lm(log(1+Eim2@sam_data$OPG)~log(1+Eim2@otu_table[,1])+log(1+Eim2@otu_table[,2])+log(1+Eim2@otu_table[,3])+log(Eim2@otu_table[,4]+1))
summary(gau.o)
cor.test(log(1+SA.e2$Abundance), log(1+SA.e2$OPG))
cor.test(log(1+SA.e2$Abundance), log(1+SA.e2$Genome_copies_gFaeces))

#ggplot(SA.e1, aes(log(1+Abundance), log(1+OPG)))+
#    geom_jitter()

#ggplot(SA.e1, aes(log(1+Abundance), log(1+Genome_copies_gFaeces)))+
#    geom_jitter()


# ploting ASV ~ GC by dpi
summary(SA.e10$dpi)
## got 10 levels but only 9 dpi's
coul <- c("#F6E8C3","#DFC27D","#BF812D","#F1B6DA", "#C51B7D", "#DE77AE","#01665E", "#35978F", "#80CDC1","#C7EAE5")

coul1 <- coul[-c(1)]
coul2 <- coul1[-1]
coul3 <- coul2[-1]

plot_SA_1 <- ggplot(SA.e10, aes(x=log(Abundance), y=log(Genome_copies_ngDNA), colour=dpi))+
    geom_point(aes(fill=dpi), shape=21, colour="black", position=position_jitter(0.2), size=4, alpha=0.8)+
    scale_color_manual(values=coul2)+
    scale_fill_manual(values=coul1)+
    labs(y="Eimeria genome copies (log)", x="Eimeria ASV1 abundance (log)")+
    ggtitle("Eimeria ASV1")+
    geom_smooth(method=lm, se=FALSE)+
    coord_cartesian(xlim=c(-8.5, 2), ylim=c(-0.5, 12.1))+
    guides(colour="none")+
    theme_bw(base_size=10)+
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=10),
          legend.position = "none",
          axis.line = element_line(colour = "black"))

plot_SA_1

SA.e20 <- SA.e20[!is.na(SA.e20$Abundance),]
plot_SA_2 <- ggplot(SA.e20, aes(x=log(Abundance), y=log(Genome_copies_ngDNA), colour=dpi))+
    geom_point(aes(fill=dpi), shape=21, colour="black", position=position_jitter(0.2), size=4, alpha=0.8)+
    scale_color_manual(values=coul1)+
    scale_fill_manual(values=coul)+
    labs(y="Eimeria genome copies (log)", x="Eimeria ASV2 abundance (log)")+
    ggtitle("Eimeria ASV2")+
    geom_smooth(method=lm, se=FALSE)+
    guides(colour="none")+
    coord_cartesian(xlim=c(-8.5, 2), ylim=c(-0.5, 12.1))+
    theme_bw()+
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=10),
          legend.position = "none",
          axis.line = element_line(colour = "black"))


plot_SA_2

cor.test(log(SA.e.g0$Abundance), log(SA.e.g0$Genome_copies_ngDNA))
cor.test(log(MA.e.g0$Abundance), log(MA.e.g0$Genome_copies_ngDNA))

plot_SA_all <- ggplot(SA.e.g0, aes(x=log(Abundance), y=log(Genome_copies_ngDNA)))+
    geom_point(aes(fill=dpi), shape=21, position=position_jitter(0.2), size=4, alpha=0.8)+
    scale_fill_manual(values=coul)+
    labs(y="Eimeria genome copies (log)", x="Eimeria ASV abundance (log)")+
    ggtitle("Single-amplicon")+
    coord_cartesian(xlim=c(-8.5, 2), ylim=c(-0.5, 12.1))+
    theme_bw(base_size=10)+
    annotate("text", x=min(log(SA.e.g0$Abundance)), y=max(log(SA.e.g0$Genome_copies_ngDNA)), hjust=0.05, label="rho=0.93, p<0.001")+
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=10),
          legend.position = "none",
          axis.line = element_line(colour = "black"))

plot_MA_all <- ggplot(MA.e.g0, aes(x=log(Abundance), y=log(Genome_copies_ngDNA)))+
    geom_point(aes(fill=dpi), shape=21, position=position_jitter(0.2), size=4, alpha=0.8)+
    scale_fill_manual(values=coul1)+
    labs(y="Eimeria genome copies (log)", x="Eimeria ASV abundance (log)")+
    ggtitle("Multi-amplicon")+
    coord_cartesian(xlim=c(-8.5, 2), ylim=c(-0.5, 12.1))+
    theme_bw(base_size=10)+
    annotate("text", x=min(log(MA.e.g0$Abundance)), y=max(log(MA.e.g0$Genome_copies_ngDNA)), hjust=0.05, label="rho=0.91, p<0.001")+
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=10),
          legend.position = "none",
          axis.line = element_line(colour = "black"))

plot_MA_1 <- ggplot(MA.e10, aes(x=log(Abundance), y=log(Genome_copies_ngDNA), colour=dpi))+
    geom_point(aes(fill=dpi), shape=21, colour="black", position=position_jitter(0.2), size=4, alpha=0.8)+
    scale_color_manual(values=coul2)+
    scale_fill_manual(values=coul2)+
    labs(y="Eimeria genome copies (log)", x="Eimeria ASV1 abundance (log)")+
    ggtitle("Eimeria ASV1")+
    geom_smooth(method=lm, se=FALSE)+
    coord_cartesian(xlim=c(-8.5, 2), ylim=c(-0.5, 12.1))+
    guides(colour="none")+
    theme_bw(base_size=10)+
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=10),
          legend.position = "none",
          axis.line = element_line(colour = "black"))

plot_MA_2 <-ggplot(MA.e20, aes(x=log(Abundance), y=log(Genome_copies_ngDNA), colour=dpi))+
    geom_point(aes(fill=dpi), shape=21, colour="black", position=position_jitter(0.2), size=4, alpha=0.8)+
    scale_color_manual(values=coul2)+
    scale_fill_manual(values=coul1)+
    labs(y="Eimeria genome copies (log)", x="Eimeria ASV2 abundance (log)")+
    ggtitle("Eimeria ASV2")+
    geom_smooth(method=lm, se=FALSE)+
    guides(colour="none")+
    coord_cartesian(xlim=c(-8.5, 2), ylim=c(-0.5, 12.1))+
    theme_bw(base_size=10)+
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=10),
          legend.position = "none",
          axis.line = element_line(colour = "black"))

plot_SA_ASV <- plot_grid(plot_SA_all,plot_SA_1, plot_SA_2, nrow=3, labels="auto") +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=2))

plot_SA_ASV


plot_MA_ASV <- plot_grid(plot_MA_all, plot_MA_1, plot_MA_2, nrow=3, labels=c("d", "e", "f"))+
        theme(panel.border = element_rect(colour = "black", fill=NA, size=2))

plot_MA_SA_ASV <- plot_grid(plot_SA_ASV, plot_MA_ASV, nrow=1)

legend <- get_legend(plot_SA_all+
          guides(fill=guide_legend(nrow=1, byrow=TRUE))+
          theme(legend.position="top"))

plot_MA_SA_ASV.l <- plot_grid(legend,plot_MA_SA_ASV, nrow=2, rel_heights=c(0.1, 1))

ggplot2::ggsave(file="fig/Eimeria_SA_MA_GC_ASVs_dpi.pdf", plot_MA_SA_ASV.l, width = 6, height = 10, dpi = 300)
ggplot2::ggsave(file="fig/Eimeria_SA_MA_GC_ASVs_dpi.png", plot_MA_SA_ASV.l, width = 6, height = 10, dpi = 350)

forPrese <-plot_grid(plot_SA_all, plot_MA_all, nrow=1)
forPrese.l <-plot_grid(legend, forPrese, nrow=2, rel_heights=c(0.1, 1))
    
ggplot2::ggsave(file="fig/FOR_PRESENT.Eimeria_SA_MA_GC_ASVs_dpi.png", forPrese.l, width = 10, height = 4, dpi = 300)
ggplot2::ggsave(file="fig/FOR_PRESENT.Eimeria_SA_MA_GC_ASVs_dpi.pdf", forPrese.l, width = 10, height = 4, dpi = 300)

### make glm qPCR~ASV1+ASV2+...
## ASV5+ASV5 useful?
#### Plot by EH_ID, only positive animals

keep <- SA.e4$EH_ID[SA.e4$Abundance>0]
SA.e4 <- SA.e4[which(SA.e4$EH_ID%in%keep),]
SA4 <- ggplot(SA.e4, aes(x=dpi, y=log(1+Abundance), fill=EH_ID))+
    geom_point(shape=21, position=position_jitter(0.2), size=2, alpha=0.7)+
    geom_line(aes(group=EH_ID), alpha=0.2)+
    labs(y="ASV4 abundance log(+1)", x="Days post infection")+
    #    scale_fill_manual(values=coul)+
    theme_bw(base_size=12)+
    ggtitle("SA")+
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          legend.position = "none",
          axis.line = element_line(colour = "black"))


keep <- SA.e3$EH_ID[SA.e3$Abundance>0]
SA.e3 <- SA.e3[which(SA.e3$EH_ID%in%keep),]
SA3 <- ggplot(SA.e3, aes(x=dpi, y=log(1+Abundance), fill=EH_ID))+
    geom_point(shape=21, position=position_jitter(0.2), size=2, alpha=0.7)+
    geom_line(aes(group=EH_ID), alpha=0.2)+
#    scale_fill_manual(values=coul)+
    labs(y="ASV3 abundance log(+1)", x="Days post infection")+
    ggtitle("SA")+
    theme_bw(base_size=12)+
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          legend.position = "none",
          axis.line = element_line(colour = "black"))


keep <- SA.e2$EH_ID[SA.e2$Abundance>0]
SA.e2 <- SA.e2[which(SA.e2$EH_ID%in%keep),]
SA2 <- ggplot(SA.e2, aes(x=dpi, y=log(1+Abundance), fill=EH_ID))+
    geom_point(shape=21, position=position_jitter(0.2), size=2, alpha=0.7)+
    geom_line(aes(group=EH_ID), alpha=0.2)+
#    scale_fill_manual(values=coul)+
    labs(y="ASV2 abundance log(+1)", x="Days post infection")+
    ggtitle("SA")+
    theme_bw(base_size=12)+
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          legend.position = "none",
          axis.line = element_line(colour = "black"))


keep <- SA.e1$EH_ID[SA.e1$Abundance>0]
SA.e1 <- SA.e1[which(SA.e1$EH_ID%in%keep),]
SA1 <- ggplot(SA.e1, aes(x=dpi, y=log(1+Abundance), fill=EH_ID))+
    geom_point(shape=21, position=position_jitter(0.2), size=2, alpha=0.7)+
    geom_line(aes(group=EH_ID), alpha=0.2)+
#    scale_fill_manual(values=coul)+
    labs(y="ASV1 abundance log(+1)", x="Days post infection")+
    ggtitle("SA")+
    theme_bw(base_size=12)+
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          legend.position = "none",
          axis.line = element_line(colour = "black"))


### ASV from MA run plotted
keep <- MA.e2$EH_ID[MA.e2$Abundance>0]
MA.e2 <- MA.e2[which(MA.e2$EH_ID%in%keep),]
MA2 <- ggplot(MA.e2, aes(x=dpi, y=log(1+Abundance), fill=EH_ID))+
    geom_point(shape=21, position=position_jitter(0.2), size=2, alpha=0.7)+
    geom_line(aes(group=EH_ID), alpha=0.2)+
#    scale_fill_manual(values=coul)+
    labs(y="ASV2 abundance log(+1)", x="Days post infection")+
    ggtitle("MA")+
    theme_bw(base_size=12)+
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          legend.position = "none",
          axis.line = element_line(colour = "black"))


keep <- MA.e1$EH_ID[MA.e1$Abundance>0]
MA.e1 <- MA.e1[which(MA.e1$EH_ID%in%keep),]
MA1 <- ggplot(MA.e1, aes(x=dpi, y=log(1+Abundance), fill=EH_ID))+
    geom_point(shape=21, position=position_jitter(0.2), size=2, alpha=0.7)+
    geom_line(aes(group=EH_ID), alpha=0.2)+
#    scale_fill_manual(values=coul)+
    labs(y="ASV1 abundance log(+1)", x="Days post infection")+
    ggtitle("MA")+
    theme_bw(base_size=12)+
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          legend.position = "none",
          axis.line = element_line(colour = "black"))


qpcr <- ggplot(MA.e1, aes(x=dpi, y=log(1+Genome_copies_gFaeces), fill=EH_ID))+
    geom_point(shape=21, position=position_jitter(0.2), size=2, alpha=0.7)+
    geom_line(aes(group=EH_ID), alpha=0.2)+
#    scale_fill_manual(values=coul)+
    ggtitle("qPCR")+
    labs(y="Eimeira Genome copies log(1+)", x="Days post infection")+
    theme_bw(base_size=10)+
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=10),
          legend.position = "none",
          axis.line = element_line(colour = "black"))

MA.all <- ggplot(MA.e.g, aes(x=dpi, y=log(1+Abundance), fill=EH_ID))+
    geom_point(shape=21, position=position_jitter(0.1), size=2, alpha=0.7)+
    geom_line(aes(group=EH_ID), alpha=0.2)+
#    scale_fill_manual(values=coul)+
    ggtitle("Multi-amplicon")+
    labs(y="Eimeria ASV abundance( log(1+)", x="Days post infection")+
        theme_bw(size_base=10)+
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=10),
          legend.position = "none",
          axis.line = element_line(colour = "black"))

SA.all <- ggplot(SA.e.g, aes(x=dpi, y=log(1+Abundance), fill=EH_ID))+
    geom_point(shape=21, position=position_jitter(0.1), size=2, alpha=0.7)+
    geom_line(aes(group=EH_ID), alpha=0.2)+
#    scale_fill_manual(values=coul)+
    ggtitle("Single amplicon")+
    labs(y="Eimeria ASV abundance log(1+)", x="Days post infection")+
    theme_bw(base_size=10)+
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=10),
          legend.position = "none",
          axis.line = element_line(colour = "black"))

MA.all

SA.all

ma.sa_1 <- SA.e[, c("ASV", "Abundance", "EH_ID", "dpi", "labels")]
ma.sa_2 <- MA.e[, c("ASV", "Abundance", "labels")]

ma.sa <- merge(ma.sa_1, ma.sa_2, by=c("ASV", "labels"))
# remove zeros

ma.sa.0 <- ma.sa[ma.sa$Abundance.x>0,]
ma.sa.0 <- ma.sa.0[ma.sa.0$Abundance.y>0,]

head(ma.sa)

ma.sa

cor.test(ma.sa$Abundance.x, ma.sa$Abundance.y, method="pearson")
cor.test(ma.sa.0$Abundance.x, ma.sa.0$Abundance.y, method="pearson")

ASV.c <- ggplot(ma.sa, aes(x=Abundance.x, y=Abundance.y, fill=ASV))+
    geom_point(size=4, shape=21, alpha=0.7)+
    scale_fill_manual(values=c("#009E73", "mediumvioletred"), name="")+
    xlab("Eimeria ASV abundance - single amplicon")+
    ylab("Eimeria ASV abundance - multi-amplicon")+
    annotate(geom="text", x=1, y=4, label="Pearson rho=0.99, p<0.001", size=3)+
    theme_bw(base_size=10)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=10),
          legend.position = "top",
          axis.line = element_line(colour = "black"))

ASV.c

sa <- SA.e[, c("ASV", "Abundance", "EH_ID", "dpi", "labels")]
sa$amp <- "sa"
ma <- MA.e[, c("ASV", "Abundance", "EH_ID", "dpi", "labels")]
ma$amp <- "ma"
sama <- rbind(sa, ma)

ASV_sama <- ggplot(sama, aes(x=(Abundance), y=ASV, fill=amp))+
    geom_point(size=4, shape=21, alpha=0.7, position=position_jitterdodge(dodge.width=0.75, jitter.width=0.1))+
    geom_boxplot(alpha=0.3, colour="black", outlier.shape = NA)+
    scale_fill_manual(values=c("#CC6677", "#DDCC77"), name="", labels=c("multi-amplicon", "single-amplicon"))+
    xlab("Eimeria ASV abundance")+
    ylab("")+
    guides(fill = guide_legend(override.aes = list(linetype = 0)),
           color = guide_legend(override.aes = list(linetype = 0)))+
    theme_bw(base_size=10)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=10),
          legend.position = "top",
          axis.line = element_line(colour = "black"))

ASV_sa

Fig1 <- plot_grid(ASV.c, ASV_sama, labels="auto", nrow=1)

ASV.c

ggplot2::ggsave(file="fig/Figure1-ASV_concordance.pdf", Fig1, width=8, height=4, dpi=300)
ggplot2::ggsave(file="fig/Figure1-ASV_concordance.png", Fig1, width=8, height=4, dpi=300)


ASV_sama

head(ma.sa)

sama_1 <- ma.sa[ma.sa$ASV=="ASV1",]

sama_2 <- ma.sa[ma.sa$ASV=="ASV2",]

t.test(log(1+sama_1$Abundance.x), log(1+sama_1$Abundance.y))
wilcox.test(sama_1$Abundance.x, sama_1$Abundance.y)

t.test(log(1+sama_2$Abundance.x), log(1+sama_2$Abundance.y))
wilcox.test(sama_2$Abundance.x, sama_2$Abundance.y)

ma.sa[ma.sa$ASV=="ASV1",]

SA_Eimeria.ASVs <- ggplot(SA.e.0, aes(y=log(Genome_copies_ngDNA), x=log(Abundance), fill=ASV))+
    geom_point(shape=21, size=4, alpha=0.7)+
    scale_fill_manual(values=c("#009E73", "#F0E442", "#0072B2",
                               "#D55E00", "#E63C9A"))+
    geom_smooth(method=lm, colour="black", aes(colour="ASV"))+
    ylab("Genome copies/ng DNA (log)")+
#    guides(fill=guide_legend(ncol=3, byrow=TRUE))+
    xlab("ASV abundance (log)")+
#    annotate(geom="text", x=12, y=19, label="Pearson rho=0.92, p<0.001", size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          legend.position = "top",
          axis.line = element_line(colour = "black"))

SA_Eimeria.ASVs

MA_Eimeria.ASVs <- ggplot(MA.e.0, aes(x=log(Genome_copies_ngDNA), y=log(Abundance), fill=ASV))+
    geom_point(shape=21, size=4, alpha=0.7)+
    scale_fill_manual(values=c("#009E73", "#F0E442", "#0072B2"))+
    geom_smooth(method=lm, colour="black", aes(colour="ASV"))+
    ylab("Genome copies/ng DNA (log)")+
    xlab("ASV abundance (log)")+
#    annotate(geom="text", x=12, y=19, label="Pearson rho=0.92, p<0.001", size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=11),
          legend.position = "top",
          axis.line = element_line(colour = "black"))

SA_Eimeria.ASVs

ASV.SA.MA <- plot_grid(SA_Eimeria.ASVs, MA_Eimeria.ASVs, ASV.c, nrow=1, labels="auto")

plot_grid(SA1,SA2,SA3,SA4) -> SA.asv
plot_grid(MA1,MA2, nrow=2) -> MA.asv

row1 <- plot_grid(SA1,MA1, nrow=1, labels=c("b", "c"))
row2 <- plot_grid(SA2,MA2, nrow=1, labels=c("d", "e"))

SAMA.asv <- plot_grid(qpcr, row1, row2, labels=c("a", "", ""), ncol=1)

SAMA.all <- plot_grid(qpcr, SA.all, MA.all, labels="auto", nrow=1)

SAMA.asv

SAMA.all

ggplot2::ggsave(file="fig/Eimeria_ASVs_dpi.pdf", SAMA.asv, width = 10, height = 15, dpi = 300)
ggplot2::ggsave(file="fig/Eimeria_ASVs_dpi.png", SAMA.asv, width = 10, height = 15, dpi = 300)

ggplot2::ggsave(file="fig/Sup1_Eimeria_SA_MA_all_dpi.pdf", SAMA.all, width = 10, height = 4, dpi = 300)

ggplot2::ggsave(file="fig/Sup1_Eimeria_SA_MA_all_dpi.png", SAMA.all, width = 10, height = 4, dpi = 300)

ggplot2::ggsave(file="fig/MA_SA_Eimeria_ASVs.pdf", ASV.SA.MA, width = 13, height = 5, dpi = 300)
ggplot2::ggsave(file="fig/MA_SA_Eimeria_ASVs.png", ASV.SA.MA, width = 13, height = 5, dpi = 300)


ggplot2::ggsave(file="fig/ASV_abundace_SAMA.pdf", ASV_sama, width=4, height=4, dpi=300)
ggplot2::ggsave(file="fig/ASV_abundace_SAMA.png", ASV_sama, width=4, height=4, dpi=300)



#################### plotting individuals by ASV
library(cowplot) # to plot a list of plots

cl <- colorRampPalette(brewer.pal(8, "Accent"))(6)

length(levels(as.factor(SA.e$EH_ID)))

SA.e$EH_ID

p.ID <- function(i){
    SA.e%>%
    dplyr::filter(EH_ID%in%SA.e$EH_ID[i])%>%
    dplyr::select(EH_ID, dpi, ASV, Genome_copies_ngDNA, Abundance)%>%
    ggplot(aes(x=dpi, y=log(Abundance+1), fill=ASV))+
    geom_point(shape=21, position=position_jitter(0.2), size=4, alpha=0.8, aes(fill=ASV), color="black")+
    geom_line(aes(group=ASV), color="gray", alpha=0.5)+
    scale_fill_manual(values=c("#009E73", "mediumvioletred","#F0E442","#0072B2"), name="")+
    ylab("Single-amplicon ASV abundance /ng DNA (log(1+)")+
    theme_bw()+
    theme(text = element_text(size=16), axis.title.x = element_blank(), legend.position = "top") -> nm
    nm
}

library(dplyr)
p.EH <- list()
for (i in 1:22){
    p <- p.ID(i)
    p.EH[[i]] <- p
}

pdf("fig/SA/Eimeria_ASVs_ID.pdf")
for (i in 1:22) {
    print(p.EH[[i]])
}
dev.off()

p.ID.MA <- function(i){
    MA.e%>%
    dplyr::filter(EH_ID%in%MA.e$EH_ID[i])%>%
    dplyr::select(EH_ID, dpi, ASV, Genome_copies_ngDNA, Abundance)%>%
    ggplot(aes(x=dpi, Abundance+1, fill=ASV))+
    geom_point(shape=21, position=position_jitter(0.2), size=4, alpha=0.7, aes(fill=ASV), color="black")+
    ylab("Multi-amplicon ASV abundance /ng DNA (log(1+)")+
    geom_line(aes(group=ASV), color="gray", alpha=0.5)+
    scale_fill_manual(values=c("#009E73", "mediumvioletred"), name="")+
    theme_bw()+
    theme(text = element_text(size=16), axis.title.x = element_blank(), legend.position = "top") -> nm
    nm
}

p.EH.M <- list()
for (i in 1:22){
    p <- p.ID.MA(i)
    p.EH.M[[i]] <- p
}

p.EH.M[1]

pdf("fig/MA/Eimeria_ASVs_ID.pdf")
for (i in 1:22) {
    print(p.EH.M[[i]])
}
dev.off()

################################



ggplot(SA.e, aes(x=Abundance, y=ASV))+
    geom_point(aes(fill=ASV), size=4, shape=21, apha=0.5, position=position_jitter(height=0.2))+
#      geom_density(aes(color = ASV), size = 1) +
    scale_fill_manual(values=c("#009E73", "#F0E442", "#0072B2", "#D55E00", "#E63C9A"))+
#        scale_colour_manual(values=c("#009E73", "#F0E442", "#0072B2", "#D55E00", "#E63C9A"))+
#    ylab("SA - ASV abundance (log1+)")+
#    xlab("MA - ASV abundance (log1+)")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
           legend.position = "top",
          axis.line = element_line(colour = "black"))


ASV1 <- ma.sa[ma.sa$ASV=="ASV1",]
ASV2 <- ma.sa[ma.sa$ASV=="ASV2",]
ASV3 <- ma.sa[ma.sa$ASV=="ASV3",]

cor.test(ASV1$Abundance.x, ASV1$Abundance.y)
cor.test(ASV2$Abundance.x, ASV2$Abundance.y)
cor.test(log(1+ASV1$Abundance.x), log(1+ASV1$Abundance.y))
cor.test(log(1+ASV2$Abundance.x), log(1+ASV2$Abundance.y))


############################
#######################################################################################
seqs <- DNAStringSet(getSequences(colnames(Eim@otu_table)))
seqs2 <- DNAStringSet(getSequences(colnames(Eim2@otu_table)))

##aligments
#first better name reads
names(seqs2) <- c("ASV1", "ASV2", "ASV3", "ASV4")

## let's get all the sequences used in Jarquín-Díaz et al. 2019


access.l <- c("JQ993669", "JQ993670","JQ993671","JQ993665", "JQ993659", "KU174470", "KU174461", "KU174464", "KU174480", "KU174483", "KU174462", "KU174463", "KU174468", "KU174479", "KU174449", "KU174450", "KU174465", "KU174467", "KU174474", "KU174451", "KU174459", "KU174485", "KU174487", "KU174472", "JQ993666","JQ993649", "JQ993650", "AF080614", "MH751998", "KT360995", "JF304148", "U40263","JQ993651", "AF311643","AF311644", "JQ993654", "JQ993655", "JQ993656", "JQ993657", "JQ993658", "JQ993660", "KU174475", "JQ993661", "JQ993662", "JQ993663", "JQ993664", "JQ993652", "JQ993667", "KU174454", "KU174469", "KU174481", "KU174456","KU174484", "KU174478", "KU174473", "KU174471",  "KU174455", "KU174457", "KU174476","KU174466", "KU174486", "KU174453","KU174458", "KU174460", "KU174452","KU174482", "AF246717", "KT184355","JQ993653", "AF307880", "AF339489","AF307878", "AF307876", "AF324214","AF339490", "AF339491", "AF307879", "AF339492", "AF307877", "AF311642", "KU192965", "KU192958", "KU192936", "KU192961", "KU192956", "KU192931", "KU192938", "KU192916")
eim.db <- read.GenBank(access.l)
eim.db.ID <- paste(names(eim.db), attr(eim.db, "species"), sep="_")
# extra 18S falciformis
Eimf <- read.GenBank("KT184339.1")
Eimf.ID <- paste(names(Eimf), attr(Eimf, "species"), sep="_")



# convert to DNAStringset
eim.DB <- eim.db %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
names(eim.DB) <- eim.db.ID

Eim.f <- Eimf %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
names(Eim.f) <- Eimf.ID

names(eim.DB)

refEim <- readDNAStringSet("/SAN/Susanas_den/AmpMarkers/wildEimeria18S/Eim_ref.fa")
names(refEim) <- gsub("(\\s)", "_", names(refEim))

#gotta correct here
which(names(refEim)=="MH751946.1_Eimeria_vermiformis")
names(refEim)[22] <- "MH751946.1_Eimeria_ferrisi"

names(seqs2) <- c("Lab_single-multi-amplicon_ASV1", "Lab_single-multi-amplicon_ASV2", "Lab_single-amplicon_ASV3", "Lab_single-amplicon_ASV4")

allSeqs <- c(eim.DB, Eim.f, refEim, seqs2)

Ref_Eims <- c(refEim, eim.DB)

alignment <- AlignSeqs(seqs2, anchor=NA, verbose=FALSE, iterations=10, refinements=10, processors=90)
Allal <- AlignSeqs(allSeqs, anchor=NA, verbose=FALSE, iterations=10, refinements=10, processors=90)
Allal <- AdjustAlignment(Allal)

Ref_Align <- AlignSeqs(Ref_Eims, anchor=NA, verbose=FALSE, iterations=20, refinements=20, processors=90)

writeFasta(Ref_Align, "tmp/Eimeria_references.fa")
writeFasta(Eim.f, "tmp/Eimeria_falKT184339.1.fa")
writeFasta(seqs2, "tmp/Eimeria_lab_ASV.fa")
writeFasta(allSeqs, "tmp/Eimeria_seqs.fa")
writeFasta(alignment, "tmp/Eimeria_reads_alignment.fa")

#BrowseSeqs(Allal, highlight=0)

# read our tree constructed with iqtree2
#~/iqtree-2.2.0-Linux/bin/iqtree2 -s tmp/Eimeria_alignment.fa -m TN+F+R2 -B 5000 -redo
